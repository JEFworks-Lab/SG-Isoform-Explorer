#!/usr/bin/env python3
"""Compile a small, human-edited project file into browser_build.json.

This preparation step resolves gene symbols from a full GTF, writes complete
gene-specific GTFs, discovers per-sample files from filename patterns, and
infers all gene IDs, loci, strands, transcript counts, and versioned IDs.
The generated build manifest is intentionally machine-facing and should not be
edited by hand.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

from gtf_tools import extract_gene_subsets, read_table_transcripts


ROOT = Path(__file__).resolve().parents[1]


def resolve(base: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def relative_to_config(path: Path, config_dir: Path) -> str:
    return os.path.relpath(path, config_dir)


def load_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {"sample_id", "condition", "replicate", "display_name"}
    missing = required - set(rows[0] if rows else [])
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if not rows:
        raise ValueError(f"{path}: sample sheet is empty")
    ids = [row["sample_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError(f"{path}: sample_id values must be unique")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=ROOT / "config" / "project.example.json")
    parser.add_argument("--output", type=Path, help="Generated build manifest; defaults to config/browser_build.generated.json")
    parser.add_argument("--subset-dir", type=Path, help="Generated GTF directory; defaults to Browser/generated/gtf")
    parser.add_argument("--report", type=Path, help="Preparation report; defaults beside the generated manifest")
    args = parser.parse_args()

    project_path = args.project.resolve()
    project_dir = project_path.parent
    project = json.loads(project_path.read_text())
    output = (args.output or project_dir / "browser_build.generated.json").resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output_dir = output.parent
    subset_dir = (args.subset_dir or ROOT / "generated" / "gtf").resolve()
    subset_dir.mkdir(parents=True, exist_ok=True)
    report_path = (args.report or output.with_name(output.stem + ".report.json")).resolve()

    symbols = project["genes"]
    if any(not isinstance(symbol, str) for symbol in symbols):
        raise ValueError("project genes must be a list of gene symbols")
    reference_gtf = resolve(project_dir, project["reference_gtf"])
    if not reference_gtf.exists():
        raise FileNotFoundError(reference_gtf)

    selection_config = project.get("transcript_selection", {"mode": "all"})
    selection = selection_config.get("mode", "all")
    table_ids = None
    table_path = None
    if selection == "table":
        table_path = resolve(project_dir, selection_config["table"])
        table_ids = read_table_transcripts(
            table_path,
            symbols,
            gene_column=selection_config.get("gene_column", "external_gene_name"),
            transcript_column=selection_config.get("transcript_column", "ensembl_transcript_id_version"),
        )
    subsets = extract_gene_subsets(reference_gtf, symbols, selection=selection, table_transcripts=table_ids)

    gene_config = []
    gene_report = []
    for subset in subsets:
        gtf_path = subset_dir / f"{subset.symbol}.gtf"
        gtf_path.write_text("".join(subset.lines))
        gene_config.append(
            {
                "gene_id": subset.gene_id,
                "gene_id_version": subset.gene_id_version,
                "gene_name": subset.symbol,
                "chrom": subset.chrom,
                "start": subset.start,
                "end": subset.end,
                "strand": subset.strand,
                "gtf": relative_to_config(gtf_path, output_dir),
            }
        )
        summary = subset.summary()
        summary["gtf"] = str(gtf_path)
        gene_report.append(summary)

    sheet_path = resolve(project_dir, project["sample_sheet"])
    sample_rows = load_samples(sheet_path)
    condition_ids = {condition["id"] for condition in project["conditions"]}
    roots = {name: resolve(project_dir, value) for name, value in project["data_roots"].items()}
    patterns = project["file_patterns"]
    sample_config = []
    sample_report = []
    for row in sample_rows:
        if row["condition"] not in condition_ids:
            raise ValueError(f"{row['sample_id']}: unknown condition {row['condition']}")
        context = {**row, "replicate": int(row["replicate"])}
        resolved_files = {}
        for key in ("bam", "bai", "sj", "star_log", "quant"):
            if row.get(key):
                path = resolve(sheet_path.parent, row[key])
            else:
                root_name = "quantification" if key == "quant" else "alignment"
                path = (roots[root_name] / patterns[key].format(**context)).resolve()
            if not path.exists():
                raise FileNotFoundError(f"{row['sample_id']} {key}: {path}")
            resolved_files[key] = path
        sample_config.append(
            {
                "sample_id": row["sample_id"],
                "condition": row["condition"],
                "replicate": int(row["replicate"]),
                "display_name": row["display_name"],
                **{key: relative_to_config(path, output_dir) for key, path in resolved_files.items()},
            }
        )
        sample_report.append(
            {
                "sample_id": row["sample_id"],
                "condition": row["condition"],
                "replicate": int(row["replicate"]),
                "files": {key: str(path) for key, path in resolved_files.items()},
            }
        )

    statistics = None
    if project.get("statistics"):
        statistics = dict(project["statistics"])
        metrics_path = resolve(project_dir, statistics.pop("metrics_csv"))
        if not metrics_path.exists():
            raise FileNotFoundError(metrics_path)
        statistics["metrics_csv"] = relative_to_config(metrics_path, output_dir)

    build = {
        "generated": {
            "by": "pipeline/prepare_project.py",
            "from": relative_to_config(project_path, output_dir),
            "note": "Machine-generated; edit the project file and sample sheet instead.",
        },
        "dataset": project["dataset"],
        "conditions": project["conditions"],
        "comparison": project["comparison"],
        "output": project.get("output", "browser-data"),
        "public_output": project.get("public_output", "frontend/public/data/dataset"),
        "coverage": project.get(
            "coverage",
            {
                "source": "bam",
                "source_bin_size": 1,
                "exclude_flags": "0xF04",
                "unique_only": True,
                "normalization": "star_unique_mapped_reads",
                "max_sampled_reads": 100,
            },
        ),
        "genes": gene_config,
        "samples": sample_config,
    }
    if statistics:
        build["statistics"] = statistics
    output.write_text(json.dumps(build, indent=2) + "\n")

    report = {
        "status": "PASS",
        "project": str(project_path),
        "generated_build_manifest": str(output),
        "reference_gtf": str(reference_gtf),
        "transcript_selection": selection,
        "selection_table": str(table_path) if table_path else None,
        "genes": gene_report,
        "samples": sample_report,
        "important": [
            "Gene IDs, versions, loci, strands, and transcript counts were inferred from the GTF.",
            "Expected junction counts are observations produced by build/validation, not user-entered requirements.",
            "Preparation does not modify BAM, BAI, STAR, Salmon, or source statistics files.",
        ],
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Prepared {len(gene_config)} genes and {len(sample_config)} samples")
    print(f"Build manifest -> {output}")
    print(f"Preparation report -> {report_path}")


if __name__ == "__main__":
    main()
