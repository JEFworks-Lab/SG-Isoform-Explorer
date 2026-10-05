#!/usr/bin/env python3
"""Create gene-specific Ensembl GTFs using gene symbols instead of IDs.

Examples:

  python pipeline/subset_gtf_by_gene.py --gtf annotation.gtf \
    --genes CAPNS1 ANXA2 --selection all --output-dir generated/gtf

  python pipeline/subset_gtf_by_gene.py --gtf annotation.gtf.gz \
    --genes CAPNS1 ANXA2 --selection table --table metrics.csv \
    --gene-column external_gene_name \
    --transcript-column ensembl_transcript_id_version \
    --output-dir generated/gtf
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from gtf_tools import extract_gene_subsets, read_table_transcripts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gtf", required=True, type=Path, help="Full Ensembl GTF or GTF.gz")
    parser.add_argument("--genes", required=True, nargs="+", help="Exact gene symbols, for example CAPNS1 ANXA2")
    parser.add_argument("--selection", choices=("all", "table", "canonical"), default="all")
    parser.add_argument("--table", type=Path, help="CSV/TSV containing transcript IDs for --selection table")
    parser.add_argument("--gene-column", default="external_gene_name")
    parser.add_argument("--transcript-column", default="ensembl_transcript_id_version")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--summary", type=Path, help="Defaults to OUTPUT_DIR/gene_subsets.summary.json")
    args = parser.parse_args()

    if not args.gtf.exists():
        raise FileNotFoundError(args.gtf)
    table_ids = None
    if args.selection == "table":
        if args.table is None:
            parser.error("--selection table requires --table")
        table_ids = read_table_transcripts(
            args.table,
            args.genes,
            gene_column=args.gene_column,
            transcript_column=args.transcript_column,
        )
    subsets = extract_gene_subsets(args.gtf, args.genes, selection=args.selection, table_transcripts=table_ids)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "source_gtf": str(args.gtf.resolve()),
        "selection": args.selection,
        "genes": [],
    }
    for subset in subsets:
        output = args.output_dir / f"{subset.symbol}.gtf"
        output.write_text("".join(subset.lines))
        summary = subset.summary()
        summary["gtf"] = str(output.resolve())
        payload["genes"].append(summary)
        print(f"{subset.symbol}: {len(subset.transcript_ids)} transcripts, {len(subset.lines)} records -> {output}")

    summary_path = args.summary or args.output_dir / "gene_subsets.summary.json"
    summary_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Summary -> {summary_path}")


if __name__ == "__main__":
    main()
