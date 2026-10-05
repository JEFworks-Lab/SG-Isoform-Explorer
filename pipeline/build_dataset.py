#!/usr/bin/env python3
"""Build browser-ready isoform assets from a declarative JSON manifest.

The browser uses 1-based inclusive genomic coordinates. BAM queries are
converted from pysam's 0-based half-open convention only at the I/O boundary.
STAR junctions are converted to exon-edge coordinates exactly as:

    donor = STAR column 2 - 1
    acceptor = STAR column 3 + 1

Inputs are configured in ``config/browser_build.json`` so adding a gene or
replacing a sample does not require editing this script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import shutil
from collections import defaultdict
from pathlib import Path

try:
    import pysam
except ImportError as exc:  # pragma: no cover
    raise SystemExit("pysam is required: python -m pip install -r requirements.txt") from exc


def dump(path: Path, value, *, compact: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compact:
        path.write_text(json.dumps(value, separators=(",", ":")) + "\n")
    else:
        path.write_text(json.dumps(value, indent=2, sort_keys=False) + "\n")


def resolve_path(config_dir: Path, value: str | None) -> Path | None:
    if value is None:
        return None
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (config_dir / path).resolve()


def load_config(path: Path) -> dict:
    config = json.loads(path.read_text())
    for key in ("dataset", "conditions", "comparison", "genes", "samples", "coverage"):
        if key not in config:
            raise ValueError(f"Configuration is missing required key: {key}")
    condition_ids = {row["id"] for row in config["conditions"]}
    if len(condition_ids) != len(config["conditions"]):
        raise ValueError("Condition IDs must be unique")
    comparison = config["comparison"]
    if {comparison["numerator"], comparison["denominator"]} - condition_ids:
        raise ValueError("Comparison conditions must be declared in conditions")

    config_dir = path.parent
    gene_ids = set()
    for gene in config["genes"]:
        if gene["gene_id"] in gene_ids:
            raise ValueError(f"Duplicate gene_id: {gene['gene_id']}")
        gene_ids.add(gene["gene_id"])
        gene["gtf"] = str(resolve_path(config_dir, gene["gtf"]))
    sample_ids = set()
    for sample in config["samples"]:
        if sample["sample_id"] in sample_ids:
            raise ValueError(f"Duplicate sample_id: {sample['sample_id']}")
        sample_ids.add(sample["sample_id"])
        if sample["condition"] not in condition_ids:
            raise ValueError(f"Unknown condition for {sample['sample_id']}: {sample['condition']}")
        for key in ("bam", "bai", "sj", "star_log", "quant"):
            sample[key] = str(resolve_path(config_dir, sample.get(key)))
    statistics = config.get("statistics")
    if statistics and statistics.get("metrics_csv"):
        statistics["metrics_csv"] = str(resolve_path(config_dir, statistics["metrics_csv"]))
    return config


def gtf_attributes(raw: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'(\S+) "([^"]+)"', raw)}


def gtf_tags(raw: str) -> list[str]:
    return re.findall(r'tag "([^"]+)"', raw)


def versioned(base: str, version: str | None) -> str:
    return f"{base}.{version}" if version else base


def interval(start: int, end: int, **extra) -> dict:
    return {"start": start, "end": end, **extra}


def parse_subset_models(gene_configs: list[dict]) -> tuple[list[dict], list[dict]]:
    genes = []
    transcripts = []
    for meta in gene_configs:
        gene_id = meta["gene_id"]
        path = Path(meta["gtf"])
        if not path.exists():
            raise FileNotFoundError(path)
        tx_by_id: dict[str, dict] = {}
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"{path}:{line_number}: expected 9 GTF columns")
                chrom, source, feature, start, end, _score, strand, frame, raw_attrs = fields
                attrs = gtf_attributes(raw_attrs)
                if attrs.get("gene_id") != gene_id:
                    raise ValueError(f"{path}:{line_number}: unexpected gene_id {attrs.get('gene_id')}")
                transcript_base = attrs.get("transcript_id")
                if not transcript_base:
                    continue
                transcript_id = versioned(transcript_base, attrs.get("transcript_version"))
                start_i, end_i = int(start), int(end)
                tx = tx_by_id.setdefault(
                    transcript_id,
                    {
                        "transcript_id": transcript_id,
                        "transcript_id_base": transcript_base,
                        "transcript_name": attrs.get("transcript_name", transcript_id),
                        "gene_id": gene_id,
                        "gene_id_version": meta.get("gene_id_version", gene_id),
                        "gene_name": meta["gene_name"],
                        "chrom": chrom,
                        "start": start_i,
                        "end": end_i,
                        "strand": strand,
                        "biotype": attrs.get("transcript_biotype"),
                        "source": source,
                        "tags": gtf_tags(raw_attrs),
                        "exons": [],
                        "cds": [],
                        "five_prime_utr": [],
                        "three_prime_utr": [],
                    },
                )
                if chrom != meta["chrom"] or strand != meta["strand"]:
                    raise ValueError(f"{path}:{line_number}: locus chromosome/strand does not match config")
                tx["start"] = min(tx["start"], start_i)
                tx["end"] = max(tx["end"], end_i)
                if feature == "exon":
                    tx["exons"].append(interval(start_i, end_i, exon_number=int(attrs.get("exon_number", 0))))
                elif feature == "CDS":
                    tx["cds"].append(interval(start_i, end_i, phase=None if frame == "." else int(frame)))
                elif feature in ("five_prime_utr", "three_prime_utr"):
                    tx[feature].append(interval(start_i, end_i))

        for tx in tx_by_id.values():
            if not tx["exons"]:
                raise ValueError(f"{path}: {tx['transcript_id']} has no exons")
            for key in ("exons", "cds", "five_prime_utr", "three_prime_utr"):
                tx[key].sort(key=lambda item: (item["start"], item["end"]))
            tx["is_mane_select"] = "MANE_Select" in tx["tags"]
            tx["is_ensembl_canonical"] = "Ensembl_canonical" in tx["tags"]
        ordered = sorted(
            tx_by_id.values(),
            key=lambda tx: (not tx["is_mane_select"], not tx["is_ensembl_canonical"], tx["transcript_name"]),
        )
        expected = meta.get("expected_transcripts")
        if expected is not None and len(ordered) != expected:
            raise AssertionError(f"{meta['gene_name']}: {len(ordered)} transcripts != configured {expected}")
        transcripts.extend(ordered)
        genes.append(
            {
                "gene_id": gene_id,
                "gene_id_version": meta.get("gene_id_version", gene_id),
                "gene_name": meta["gene_name"],
                "chrom": meta["chrom"],
                "start": int(meta["start"]),
                "end": int(meta["end"]),
                "strand": meta["strand"],
                "transcript_ids": [tx["transcript_id"] for tx in ordered],
                # Browser JSON is deployable/public output. Record the source
                # filename without leaking a maintainer's absolute filesystem
                # path; the full resolved path remains in the local build
                # configuration/report used for auditing.
                "display_model_source": path.name,
                "coordinate_system": "1-based inclusive",
            }
        )
    return genes, transcripts


def parse_star_unique_reads(log_path: Path) -> int:
    for line in log_path.read_text().splitlines():
        if "Uniquely mapped reads number" in line:
            return int(line.split("|")[1].strip())
    raise ValueError(f"Missing unique-read count in {log_path}")


def build_samples(sample_configs: list[dict], coverage_config: dict) -> list[dict]:
    rows = []
    for configured in sample_configs:
        sample = dict(configured)
        required = [Path(sample[key]) for key in ("bam", "bai", "sj", "quant")]
        if coverage_config.get("normalization") == "star_unique_mapped_reads":
            required.append(Path(sample["star_log"]))
        for path in required:
            if not path.exists():
                raise FileNotFoundError(path)
        if coverage_config.get("normalization") == "star_unique_mapped_reads":
            sample["normalization_denominator"] = parse_star_unique_reads(Path(sample["star_log"]))
            sample["normalization_denominator_name"] = "STAR uniquely mapped reads"
        elif coverage_config.get("normalization") == "explicit":
            sample["normalization_denominator"] = int(sample["normalization_denominator"])
            sample["normalization_denominator_name"] = sample.get("normalization_denominator_name", "configured alignments")
        else:
            raise ValueError("coverage.normalization must be star_unique_mapped_reads or explicit")
        rows.append(sample)
    return rows


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "t", "yes", "y"}


def load_metric_rows(statistics: dict | None, wanted: set[str]) -> dict[str, dict]:
    if not statistics or not statistics.get("metrics_csv"):
        return {}
    fields = statistics.get("fields", {})
    result = {}
    with Path(statistics["metrics_csv"]).open(newline="") as handle:
        for row in csv.DictReader(handle):
            transcript_id = row[statistics["transcript_id_column"]]
            if transcript_id not in wanted:
                continue
            parsed = {}
            for output_name, source_name in fields.items():
                raw = row.get(source_name, "")
                if output_name == "filter_keep":
                    parsed[output_name] = parse_bool(raw)
                else:
                    parsed[output_name] = None if raw == "" else float(raw)
            result[transcript_id] = parsed
    missing = wanted - result.keys()
    if missing:
        raise ValueError(f"Displayed transcripts missing metrics rows: {sorted(missing)}")
    return result


def build_expression(samples, transcripts, statistics, conditions, comparison) -> tuple[dict, dict]:
    wanted = {tx["transcript_id"] for tx in transcripts}
    metrics = load_metric_rows(statistics, wanted)
    significance = (statistics or {}).get("significance", {})
    effect_field = significance.get("effect_field", "swish_log2fc")
    qvalue_field = significance.get("qvalue_field", "qvalue")
    absolute_log2fc_gt = float(significance.get("absolute_log2fc_gt", 1.0))
    qvalue_lt = float(significance.get("qvalue_lt", 0.01))
    quant_by_sample = {}
    for sample in samples:
        values = {}
        with Path(sample["quant"]).open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["Name"] in wanted:
                    values[row["Name"]] = {"tpm": float(row["TPM"]), "num_reads": float(row["NumReads"])}
        if set(values) != wanted:
            raise ValueError(f"{sample['sample_id']} is missing transcripts: {sorted(wanted - values.keys())}")
        quant_by_sample[sample["sample_id"]] = values

    condition_ids = [condition["id"] for condition in conditions]
    numerator = comparison["numerator"]
    denominator = comparison["denominator"]
    expression_rows = []
    enrichment_rows = []
    for tx in transcripts:
        transcript_id = tx["transcript_id"]
        sample_values = {
            sample["sample_id"]: {
                "condition": sample["condition"],
                **quant_by_sample[sample["sample_id"]][transcript_id],
            }
            for sample in samples
        }
        condition_values = {}
        for condition in condition_ids:
            ids = [sample["sample_id"] for sample in samples if sample["condition"] == condition]
            if not ids:
                raise ValueError(f"No samples found for condition {condition}")
            values = [sample_values[sample_id]["tpm"] for sample_id in ids]
            condition_values[condition] = {"mean_tpm": sum(values) / len(values), "replicate_tpm": values}
        metric = metrics.get(transcript_id, {})
        effect = metric.get(effect_field)
        qvalue = metric.get(qvalue_field)
        metrics_available = (
            effect is not None and qvalue is not None
            and math.isfinite(effect) and math.isfinite(qvalue)
        )
        significant = bool(
            metrics_available
            and abs(effect) > absolute_log2fc_gt
            and qvalue < qvalue_lt
        )
        call = (
            "enriched" if significant and effect > 0
            else "depleted" if significant and effect < 0
            else "not_significant" if metrics_available
            else "not_tested"
        )
        metric_with_call = {
            **metric,
            "metrics_available": metrics_available,
            "significant": significant,
            "call": call,
        }
        expression_rows.append(
            {
                "gene_id": tx["gene_id"],
                "transcript_id": transcript_id,
                "transcript_name": tx["transcript_name"],
                "samples": sample_values,
                "conditions": condition_values,
                "statistics": metric_with_call,
            }
        )
        numerator_mean = condition_values[numerator]["mean_tpm"]
        denominator_mean = condition_values[denominator]["mean_tpm"]
        ratio = (numerator_mean + 0.01) / (denominator_mean + 0.01)
        enrichment_rows.append(
            {
                "gene_id": tx["gene_id"],
                "transcript_id": transcript_id,
                "transcript_name": tx["transcript_name"],
                "numerator_condition": numerator,
                "denominator_condition": denominator,
                "numerator_mean_tpm": numerator_mean,
                "denominator_mean_tpm": denominator_mean,
                "ratio": ratio,
                "log2_ratio": math.log2(ratio),
                **metric_with_call,
            }
        )

    gene_rows = []
    for gene_id in dict.fromkeys(tx["gene_id"] for tx in transcripts):
        rows = [row for row in expression_rows if row["gene_id"] == gene_id]
        numerator_mean = sum(row["conditions"][numerator]["mean_tpm"] for row in rows)
        denominator_mean = sum(row["conditions"][denominator]["mean_tpm"] for row in rows)
        ratio = (numerator_mean + 0.01) / (denominator_mean + 0.01)
        gene_rows.append(
            {
                "gene_id": gene_id,
                "scope": "sum of displayed subset transcripts",
                "numerator_condition": numerator,
                "denominator_condition": denominator,
                "numerator_mean_tpm": numerator_mean,
                "denominator_mean_tpm": denominator_mean,
                "ratio": ratio,
                "log2_ratio": math.log2(ratio),
            }
        )
    return (
        {"unit": "TPM", "conditions": condition_ids, "transcripts": expression_rows},
        {
            "comparison": comparison,
            "definition": "Transcript statistics are supplied analysis results; gene summary is descriptive across displayed subset transcripts.",
            "call_definition": {
                "effect_field": effect_field,
                "qvalue_field": qvalue_field,
                "absolute_log2fc_gt": absolute_log2fc_gt,
                "qvalue_lt": qvalue_lt,
                "rule": f"abs({effect_field}) > {absolute_log2fc_gt:g} and {qvalue_field} < {qvalue_lt:g}",
                "positive_call": f"{comparison['numerator']} enriched",
                "negative_call": f"{comparison['numerator']} depleted",
            },
            "genes": gene_rows,
            "transcripts": enrichment_rows,
        },
    )


def include_read(read, excluded_flags: int, unique_only: bool) -> bool:
    if read.flag & excluded_flags:
        return False
    return not unique_only or (read.has_tag("NH") and read.get_tag("NH") == 1)


def read_segments_and_junctions(read) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    ref = read.reference_start
    segment_start = None
    segments = []
    junctions = []
    for operation, length in read.cigartuples or []:
        if operation in (0, 7, 8):
            if segment_start is None:
                segment_start = ref
            ref += length
        elif operation == 2:
            if segment_start is not None:
                segments.append((segment_start + 1, ref))
                segment_start = None
            ref += length
        elif operation == 3:
            if segment_start is not None:
                segments.append((segment_start + 1, ref))
                segment_start = None
            donor = ref
            ref += length
            junctions.append((donor, ref + 1))
    if segment_start is not None:
        segments.append((segment_start + 1, ref))
    return segments, junctions


def overlaps(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return max(start_a, start_b) <= min(end_a, end_b)


def reservoir_add(items: list, value, seen: int, rng: random.Random, limit: int) -> None:
    if len(items) < limit:
        items.append(value)
    else:
        slot = rng.randrange(seen)
        if slot < limit:
            items[slot] = value


def extract_bam_assets(bam_path: Path, sample: dict, gene: dict, coverage_config: dict) -> tuple[dict, dict]:
    if coverage_config.get("source") != "bam":
        raise ValueError("This builder currently supports coverage.source=bam; BigWig is documented as a future adapter")
    if int(coverage_config.get("source_bin_size", 1)) != 1:
        raise ValueError("v3 requires exact 1-bp source coverage; set coverage.source_bin_size to 1")
    excluded_flags = int(str(coverage_config.get("exclude_flags", "0xF04")), 0)
    unique_only = bool(coverage_config.get("unique_only", True))
    max_reads = int(coverage_config.get("max_sampled_reads", 100))
    predicate = lambda read: include_read(read, excluded_flags, unique_only)
    start0 = gene["start"] - 1
    end0 = gene["end"]
    with pysam.AlignmentFile(bam_path, "rb", index_filename=sample["bai"]) as bam:
        if not bam.has_index():
            raise ValueError(f"Missing usable index for {bam_path}")
        if bam.header.get("HD", {}).get("SO") != "coordinate":
            raise ValueError(f"BAM is not coordinate-declared: {bam_path}")
        if gene["chrom"] not in bam.references:
            raise ValueError(f"{gene['chrom']} missing from {bam_path}")
        arrays = bam.count_coverage(gene["chrom"], start0, end0, quality_threshold=0, read_callback=predicate)
        raw_values = [sum(base_counts) for base_counts in zip(*arrays)]
        denominator = sample["normalization_denominator"]
        normalized_values = [value * 1_000_000 / denominator for value in raw_values]
        positions = list(range(gene["start"], gene["end"] + 1))

        seed = int.from_bytes(
            hashlib.sha256(f"{sample['sample_id']}:{gene['gene_id']}".encode()).digest()[:8], "big"
        )
        rng = random.Random(seed)
        reads = []
        spliced_reads = []
        total_reads = 0
        total_spliced = 0
        for read in bam.fetch(gene["chrom"], start0, end0):
            if not predicate(read):
                continue
            segments, junctions = read_segments_and_junctions(read)
            clipped = [
                {"start": max(start, gene["start"]), "end": min(end, gene["end"])}
                for start, end in segments
                if overlaps(start, end, gene["start"], gene["end"])
            ]
            if not clipped:
                continue
            within_junctions = [
                {"donor": donor, "acceptor": acceptor}
                for donor, acceptor in junctions
                if gene["start"] <= donor < acceptor <= gene["end"]
            ]
            total_reads += 1
            item = {
                "read_id": read.query_name,
                "mapq": read.mapping_quality,
                "strand": "-" if read.is_reverse else "+",
                "segments": clipped,
                "junctions": within_junctions,
            }
            reservoir_add(reads, item, total_reads, rng, max_reads)
            if within_junctions:
                total_spliced += 1
                reservoir_add(spliced_reads, item, total_spliced, rng, max_reads)

    filter_text = f"{'NH=1 and ' if unique_only else ''}exclude SAM flags {hex(excluded_flags)}"
    coverage = {
        "gene_id": gene["gene_id"],
        "sample_id": sample["sample_id"],
        "chrom": gene["chrom"],
        "start": gene["start"],
        "end": gene["end"],
        "coordinate_system": "1-based inclusive",
        "bin_size": 1,
        "positions": positions,
        "values": raw_values,
        "normalized_values": normalized_values,
        "unit": "unique-read depth",
        "normalized_unit": "unique-read depth per million STAR uniquely mapped reads",
        "normalization_denominator": denominator,
        "normalization_denominator_name": sample["normalization_denominator_name"],
        "read_filter": filter_text,
    }
    read_json = {
        "gene_id": gene["gene_id"],
        "sample_id": sample["sample_id"],
        "chrom": gene["chrom"],
        "start": gene["start"],
        "end": gene["end"],
        "coordinate_system": "1-based inclusive",
        "read_filter": filter_text,
        "total_reads": total_reads,
        "displayed_reads": len(reads),
        "reads": reads,
        "total_spliced_reads": total_spliced,
        "displayed_spliced_reads": len(spliced_reads),
        "spliced_reads": spliced_reads,
    }
    return coverage, read_json


def subset_junction_index(transcripts: list[dict]) -> dict:
    by_gene = defaultdict(lambda: defaultdict(list))
    for tx in transcripts:
        exons = sorted(tx["exons"], key=lambda exon: exon["start"])
        for left, right in zip(exons, exons[1:]):
            by_gene[tx["gene_id"]][(left["end"], right["start"])].append(tx["transcript_id"])
    return by_gene


def parse_star_junctions(path: Path, sample_id: str, gene: dict, model_junctions: dict) -> dict:
    junctions = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{line_number}: expected 9 columns")
            chrom = fields[0]
            intron_start, intron_end = int(fields[1]), int(fields[2])
            if chrom != gene["chrom"] or intron_start < gene["start"] or intron_end > gene["end"]:
                continue
            donor = intron_start - 1
            acceptor = intron_end + 1
            if not (gene["start"] <= donor < acceptor <= gene["end"]):
                raise AssertionError(f"Invalid converted junction {path}:{line_number}")
            matching = sorted(model_junctions.get((donor, acceptor), []))
            unique_count = int(fields[6])
            multi_mapping_count = int(fields[7])
            junctions.append(
                {
                    "donor": donor,
                    "acceptor": acceptor,
                    # ``count`` remains the compatibility field used by the
                    # frontend; both explicit names preserve STAR semantics.
                    "count": unique_count,
                    "unique_mapping_count": unique_count,
                    "multi_mapping_count": multi_mapping_count,
                    "evidence_class": "unique_supported" if unique_count > 0 else "multi_mapping_only",
                    "maximum_overhang": int(fields[8]),
                    "star_intron_start": intron_start,
                    "star_intron_end": intron_end,
                    "strand_code": int(fields[3]),
                    "motif_code": int(fields[4]),
                    "star_annotated": bool(int(fields[5])),
                    "matches_subset_annotation": bool(matching),
                    "matching_transcript_ids": matching,
                }
            )
    return {
        "gene_id": gene["gene_id"],
        "sample_id": sample_id,
        "coordinate_system": "1-based inclusive exon-edge coordinates",
        "source": "final STAR SJ.out.tab",
        "conversion": "donor = STAR column 2 - 1; acceptor = STAR column 3 + 1",
        "junctions": junctions,
    }


def build_manifest(config: dict) -> dict:
    coverage = config["coverage"]
    unique = bool(coverage.get("unique_only", True))
    flags = coverage.get("exclude_flags", "0xF04")
    return {
        "schema_version": "3.1.0",
        "dataset": config["dataset"],
        "conditions": config["conditions"],
        "comparison": config["comparison"],
        "coordinate_system": "1-based inclusive",
        "junction_convention": "donor = STAR column 2 - 1; acceptor = STAR column 3 + 1",
        "coverage": {
            "source": "indexed coordinate-sorted BAM",
            "source_bin_size": 1,
            "raw_unit": "unique-read depth" if unique else "alignment depth",
            "normalized_unit": "depth per million STAR uniquely mapped reads",
            "filter": f"{'NH=1 and ' if unique else ''}exclude SAM flags {flags}",
            "display_bin_sizes": [1, 5, 25],
        },
        "samples": "samples.json",
        "resources": {
            "genes": "genes.json",
            "transcripts": "transcripts.json",
            "expression": "expression.json",
            "enrichment": "enrichment.json",
            "coverage_pattern": "coverage/{gene_id}__{sample_id}.json",
            "junction_pattern": "junctions/{gene_id}__{sample_id}.json",
            "reads_pattern": "reads/{gene_id}__{sample_id}.json",
        },
    }


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=project_root / "config" / "browser_build.json")
    parser.add_argument("--no-public-sync", action="store_true", help="Build canonical browser-data only")
    args = parser.parse_args()
    config = load_config(args.config.resolve())

    genes, transcripts = parse_subset_models(config["genes"])
    samples = build_samples(config["samples"], config["coverage"])
    expression, enrichment = build_expression(
        samples, transcripts, config.get("statistics"), config["conditions"], config["comparison"]
    )
    model_junctions = subset_junction_index(transcripts)

    output = (project_root / config.get("output", "browser-data")).resolve()
    public = (project_root / config.get("public_output", "frontend/public/data/dataset")).resolve()
    temporary = project_root / ".browser-data.tmp"
    if temporary.exists():
        shutil.rmtree(temporary)
    temporary.mkdir(parents=True)
    for name in ("coverage", "junctions", "reads"):
        (temporary / name).mkdir()

    dump(temporary / "manifest.json", build_manifest(config))
    dump(
        temporary / "samples.json",
        [
            {
                key: sample[key]
                for key in (
                    "sample_id",
                    "condition",
                    "replicate",
                    "display_name",
                    "normalization_denominator",
                    "normalization_denominator_name",
                )
            }
            for sample in samples
        ],
    )
    dump(temporary / "genes.json", genes)
    dump(temporary / "transcripts.json", transcripts)
    dump(temporary / "expression.json", expression)
    dump(temporary / "enrichment.json", enrichment)

    gene_config_by_id = {gene["gene_id"]: gene for gene in config["genes"]}
    for sample in samples:
        for gene in genes:
            print(f"Processing {gene['gene_name']} / {sample['sample_id']} ...", flush=True)
            coverage, reads = extract_bam_assets(Path(sample["bam"]), sample, gene, config["coverage"])
            junctions = parse_star_junctions(
                Path(sample["sj"]), sample["sample_id"], gene, model_junctions[gene["gene_id"]]
            )
            expected = gene_config_by_id[gene["gene_id"]].get("expected_junction_counts", {}).get(sample["sample_id"])
            if expected is not None and len(junctions["junctions"]) != expected:
                raise AssertionError(
                    f"{sample['sample_id']} {gene['gene_name']}: {len(junctions['junctions'])} junctions != configured {expected}"
                )
            stem = f"{gene['gene_id']}__{sample['sample_id']}.json"
            dump(temporary / "coverage" / stem, coverage, compact=True)
            dump(temporary / "junctions" / stem, junctions)
            dump(temporary / "reads" / stem, reads)

    if output.exists():
        shutil.rmtree(output)
    temporary.rename(output)
    if not args.no_public_sync:
        if public.exists():
            shutil.rmtree(public)
        shutil.copytree(output, public)
    print(
        f"Built {len(genes)} genes, {len(transcripts)} transcripts, and {len(genes) * len(samples)} sample/locus tracks in {output}"
    )
    if not args.no_public_sync:
        print(f"Synchronized browser assets to {public}")


if __name__ == "__main__":
    main()
