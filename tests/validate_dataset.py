#!/usr/bin/env python3
"""Independently validate browser JSON against configured GTF/Salmon/BAM/STAR sources."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]


def load(path: Path):
    return json.loads(path.read_text())


def resolved(config_dir: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (config_dir / path).resolve()


def attrs(raw: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'(\S+) "([^"]+)"', raw)}


def near(actual: float, expected: float, label: str, tolerance: float = 1e-10) -> None:
    assert math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance), (
        f"{label}: got {actual}, expected {expected}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "browser_build.json")
    args = parser.parse_args()
    config_path = args.config.resolve()
    config_dir = config_path.parent
    config = load(config_path)
    output = (ROOT / config.get("output", "browser-data")).resolve()
    public = (ROOT / config.get("public_output", "frontend/public/data/dataset")).resolve()

    manifest = load(output / "manifest.json")
    samples = load(output / "samples.json")
    genes = load(output / "genes.json")
    transcripts = load(output / "transcripts.json")
    expression = load(output / "expression.json")
    enrichment = load(output / "enrichment.json")
    sample_config = {row["sample_id"]: row for row in config["samples"]}
    gene_config = {row["gene_id"]: row for row in config["genes"]}

    assert manifest["schema_version"] == "3.1.0"
    assert manifest["coordinate_system"] == "1-based inclusive"
    assert manifest["junction_convention"] == "donor = STAR column 2 - 1; acceptor = STAR column 3 + 1"
    assert manifest["coverage"]["source_bin_size"] == 1
    assert {gene["gene_id"] for gene in genes} == set(gene_config)
    assert {sample["sample_id"] for sample in samples} == set(sample_config)

    source_features = defaultdict(lambda: defaultdict(list))
    source_meta = {}
    for configured in config["genes"]:
        gtf = resolved(config_dir, configured["gtf"])
        with gtf.open() as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                columns = line.rstrip("\n").split("\t")
                meta = attrs(columns[8])
                transcript_base = meta.get("transcript_id")
                if not transcript_base:
                    continue
                transcript_id = transcript_base + ("." + meta["transcript_version"] if meta.get("transcript_version") else "")
                feature = columns[2]
                source_meta[transcript_id] = (columns[0], columns[6], configured["gene_id"])
                if feature in ("exon", "CDS", "five_prime_utr", "three_prime_utr"):
                    source_features[transcript_id][feature].append((int(columns[3]), int(columns[4])))
    built_transcript_ids = {tx["transcript_id"] for tx in transcripts}
    source_transcript_ids = set(source_meta)
    assert built_transcript_ids == source_transcript_ids, (
        "Transcript set mismatch: generated browser JSON contains "
        f"{len(built_transcript_ids)} models, but the GTFs selected by this config contain "
        f"{len(source_transcript_ids)}. Validate with the same config used to build the JSON. "
        f"Only in JSON: {sorted(built_transcript_ids - source_transcript_ids)[:8]}; "
        f"only in configured GTFs: {sorted(source_transcript_ids - built_transcript_ids)[:8]}"
    )
    for tx in transcripts:
        transcript_id = tx["transcript_id"]
        assert (tx["chrom"], tx["strand"], tx["gene_id"]) == source_meta[transcript_id]
        for source_name, json_name in (
            ("exon", "exons"),
            ("CDS", "cds"),
            ("five_prime_utr", "five_prime_utr"),
            ("three_prime_utr", "three_prime_utr"),
        ):
            actual = sorted((item["start"], item["end"]) for item in tx[json_name])
            expected = sorted(source_features[transcript_id][source_name])
            assert actual == expected, f"{transcript_id} {source_name}"
    for gene in genes:
        expected_count = gene_config[gene["gene_id"]].get("expected_transcripts")
        if expected_count is not None:
            assert len(gene["transcript_ids"]) == expected_count

    transcript_ids = {tx["transcript_id"] for tx in transcripts}
    quant = {}
    for sample in samples:
        source = sample_config[sample["sample_id"]]
        with resolved(config_dir, source["quant"]).open(newline="") as handle:
            quant[sample["sample_id"]] = {
                row["Name"]: float(row["TPM"])
                for row in csv.DictReader(handle, delimiter="\t")
                if row["Name"] in transcript_ids
            }
        assert set(quant[sample["sample_id"]]) == transcript_ids

    expression_by_id = {row["transcript_id"]: row for row in expression["transcripts"]}
    for transcript_id in transcript_ids:
        row = expression_by_id[transcript_id]
        for sample in samples:
            sample_id = sample["sample_id"]
            near(row["samples"][sample_id]["tpm"], quant[sample_id][transcript_id], f"TPM/{sample_id}/{transcript_id}")
        for condition in (item["id"] for item in config["conditions"]):
            values = [
                quant[sample["sample_id"]][transcript_id]
                for sample in samples
                if sample["condition"] == condition
            ]
            near(row["conditions"][condition]["mean_tpm"], sum(values) / len(values), f"mean/{condition}/{transcript_id}")

    statistics = config.get("statistics")
    if statistics and statistics.get("metrics_csv"):
        significance = statistics.get("significance", {})
        effect_field = significance.get("effect_field", "swish_log2fc")
        qvalue_field = significance.get("qvalue_field", "qvalue")
        absolute_log2fc_gt = float(significance.get("absolute_log2fc_gt", 1.0))
        qvalue_lt = float(significance.get("qvalue_lt", 0.01))
        assert enrichment["call_definition"] == {
            "effect_field": effect_field,
            "qvalue_field": qvalue_field,
            "absolute_log2fc_gt": absolute_log2fc_gt,
            "qvalue_lt": qvalue_lt,
            "rule": f"abs({effect_field}) > {absolute_log2fc_gt:g} and {qvalue_field} < {qvalue_lt:g}",
            "positive_call": f"{config['comparison']['numerator']} enriched",
            "negative_call": f"{config['comparison']['numerator']} depleted",
        }
        id_column = statistics["transcript_id_column"]
        with resolved(config_dir, statistics["metrics_csv"]).open(newline="") as handle:
            metrics = {
                row[id_column]: row
                for row in csv.DictReader(handle)
                if row[id_column] in transcript_ids
            }
        assert set(metrics) == transcript_ids
        enrichment_by_id = {row["transcript_id"]: row for row in enrichment["transcripts"]}
        expression_by_id = {row["transcript_id"]: row for row in expression["transcripts"]}
        for transcript_id, source in metrics.items():
            observed = enrichment_by_id[transcript_id]
            for output_name, source_name in statistics.get("fields", {}).items():
                if output_name == "filter_keep":
                    assert observed[output_name] == (source[source_name].lower() == "true")
                else:
                    near(observed[output_name], float(source[source_name]), f"metric/{output_name}/{transcript_id}")
            effect = observed[effect_field]
            qvalue = observed[qvalue_field]
            expected_significant = abs(effect) > absolute_log2fc_gt and qvalue < qvalue_lt
            expected_call = (
                "enriched" if expected_significant and effect > 0
                else "depleted" if expected_significant and effect < 0
                else "not_significant"
            )
            assert observed["metrics_available"] is True
            assert observed["significant"] is expected_significant
            assert observed["call"] == expected_call
            expression_stats = expression_by_id[transcript_id]["statistics"]
            assert expression_stats["metrics_available"] is True
            assert expression_stats["significant"] is expected_significant
            assert expression_stats["call"] == expected_call

    transcripts_by_gene = defaultdict(list)
    for tx in transcripts:
        transcripts_by_gene[tx["gene_id"]].append(tx)
    excluded_flags = int(str(config["coverage"].get("exclude_flags", "0xF04")), 0)
    unique_only = bool(config["coverage"].get("unique_only", True))
    checked_tracks = 0
    for sample in samples:
        sample_id = sample["sample_id"]
        source = sample_config[sample_id]
        bam_path = resolved(config_dir, source["bam"])
        sj_path = resolved(config_dir, source["sj"])
        with pysam.AlignmentFile(bam_path, "rb", index_filename=str(resolved(config_dir, source["bai"]))) as bam:
            assert bam.has_index()
            assert bam.header.get("HD", {}).get("SO") == "coordinate"
            for gene in genes:
                stem = f"{gene['gene_id']}__{sample_id}.json"
                coverage = load(output / "coverage" / stem)
                junction_json = load(output / "junctions" / stem)
                reads = load(output / "reads" / stem)
                expected_positions = list(range(gene["start"], gene["end"] + 1))
                assert coverage["bin_size"] == 1
                assert coverage["positions"] == expected_positions
                assert len(coverage["values"]) == len(coverage["normalized_values"]) == len(expected_positions)

                arrays = bam.count_coverage(
                    gene["chrom"],
                    gene["start"] - 1,
                    gene["end"],
                    quality_threshold=0,
                    read_callback=lambda read: not (read.flag & excluded_flags)
                    and (not unique_only or (read.has_tag("NH") and read.get_tag("NH") == 1)),
                )
                depth = [sum(base_counts) for base_counts in zip(*arrays)]
                assert coverage["values"] == depth, f"raw coverage/{stem}"
                denominator = coverage["normalization_denominator"]
                assert denominator == sample["normalization_denominator"]
                for index, raw in enumerate(depth):
                    near(
                        coverage["normalized_values"][index],
                        raw * 1_000_000 / denominator,
                        f"normalized coverage/{stem}/{index}",
                    )

                expected_star = []
                with sj_path.open() as handle:
                    for line in handle:
                        fields = line.rstrip("\n").split("\t")
                        if fields[0] == gene["chrom"] and int(fields[1]) >= gene["start"] and int(fields[2]) <= gene["end"]:
                            expected_star.append(
                                (int(fields[1]) - 1, int(fields[2]) + 1, int(fields[6]), int(fields[7]))
                            )
                observed_star = [
                    (j["donor"], j["acceptor"], j["unique_mapping_count"], j["multi_mapping_count"])
                    for j in junction_json["junctions"]
                ]
                assert observed_star == expected_star, f"STAR junction conversion/{stem}"
                expected_count = gene_config[gene["gene_id"]].get("expected_junction_counts", {}).get(sample_id)
                if expected_count is not None:
                    assert len(observed_star) == expected_count

                modeled = set()
                for tx in transcripts_by_gene[gene["gene_id"]]:
                    exons = sorted(tx["exons"], key=lambda exon: exon["start"])
                    modeled.update((left["end"], right["start"]) for left, right in zip(exons, exons[1:]))
                for junction in junction_json["junctions"]:
                    assert junction["count"] == junction["unique_mapping_count"]
                    assert junction["evidence_class"] == (
                        "unique_supported" if junction["unique_mapping_count"] > 0 else "multi_mapping_only"
                    )
                    if junction["evidence_class"] == "multi_mapping_only":
                        assert junction["multi_mapping_count"] > 0
                    key = (junction["donor"], junction["acceptor"])
                    assert junction["matches_subset_annotation"] == (key in modeled)
                    if junction["matches_subset_annotation"]:
                        assert key in modeled, f"arc misses displayed exon edge/{stem}/{key}"
                assert reads["total_reads"] >= reads["displayed_reads"] > 0
                assert reads["total_spliced_reads"] >= reads["displayed_spliced_reads"] > 0
                checked_tracks += 1

    output_files = sorted(output.rglob("*.json"))
    for source in output_files:
        served = public / source.relative_to(output)
        assert served.exists() and served.read_bytes() == source.read_bytes(), f"stale served data: {served}"

    print(
        f"PASS: {len(genes)} genes, {len(transcripts)} subset models, "
        f"{len(expression['transcripts'])} expression rows, {checked_tracks} BAM/STAR tracks, "
        f"{len(output_files)} synchronized JSON files."
    )


if __name__ == "__main__":
    main()
