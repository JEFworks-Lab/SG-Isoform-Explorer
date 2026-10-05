"""Small, dependency-free helpers for selecting genes from an Ensembl GTF."""

from __future__ import annotations

import csv
import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


ATTRIBUTE_RE = re.compile(r'(\S+) "([^"]+)"')


def attributes(raw: str) -> dict[str, str]:
    return {key: value for key, value in ATTRIBUTE_RE.findall(raw)}


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def versioned(base: str, version: str | None) -> str:
    return f"{base}.{version}" if version else base


@dataclass
class GeneSubset:
    symbol: str
    gene_id: str
    gene_id_version: str
    chrom: str
    start: int
    end: int
    strand: str
    transcript_ids: list[str]
    lines: list[str]
    source_record_count: int

    def summary(self) -> dict:
        return {
            "gene_name": self.symbol,
            "gene_id": self.gene_id,
            "gene_id_version": self.gene_id_version,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "transcript_count": len(self.transcript_ids),
            "transcript_ids": self.transcript_ids,
            "record_count": len(self.lines),
            "source_record_count": self.source_record_count,
        }


def read_table_transcripts(
    table: Path,
    symbols: Iterable[str],
    *,
    gene_column: str,
    transcript_column: str,
) -> dict[str, set[str]]:
    """Return versioned transcript IDs from a CSV/TSV, keyed by requested symbol."""
    requested = {symbol.casefold(): symbol for symbol in symbols}
    delimiter = "\t" if table.suffix.lower() in {".tsv", ".txt"} else ","
    selected = {symbol: set() for symbol in requested.values()}
    with table.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        columns = set(reader.fieldnames or [])
        missing = {gene_column, transcript_column} - columns
        if missing:
            raise ValueError(f"{table}: missing columns {sorted(missing)}")
        for row in reader:
            requested_symbol = requested.get(row[gene_column].strip().casefold())
            transcript_id = row[transcript_column].strip()
            if requested_symbol and transcript_id:
                selected[requested_symbol].add(transcript_id)
    empty = [symbol for symbol, ids in selected.items() if not ids]
    if empty:
        raise ValueError(f"No transcript rows found in {table} for: {', '.join(empty)}")
    return selected


def extract_gene_subsets(
    gtf: Path,
    symbols: Iterable[str],
    *,
    selection: str = "all",
    table_transcripts: dict[str, set[str]] | None = None,
) -> list[GeneSubset]:
    """Scan a GTF once and return records for each requested gene symbol.

    ``selection`` is ``all``, ``table``, or ``canonical``. Canonical retains
    MANE Select and Ensembl canonical models (normally the same transcript).
    Table selection accepts versioned or unversioned IDs and always emits the
    complete set of GTF records for each selected transcript plus the gene row.
    """
    requested = list(dict.fromkeys(symbols))
    if not requested:
        raise ValueError("At least one gene symbol is required")
    if selection not in {"all", "table", "canonical"}:
        raise ValueError("selection must be one of: all, table, canonical")
    if selection == "table" and table_transcripts is None:
        raise ValueError("table selection requires transcript IDs")

    symbol_pattern = re.compile(
        r'gene_name "(' + "|".join(re.escape(symbol) for symbol in requested) + r')";',
        re.IGNORECASE,
    )
    requested_by_fold = {symbol.casefold(): symbol for symbol in requested}
    records: dict[str, list[tuple[list[str], dict[str, str], str]]] = {symbol: [] for symbol in requested}
    gene_ids: dict[str, set[str]] = {symbol: set() for symbol in requested}

    with open_text(gtf) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            match = symbol_pattern.search(line)
            if not match:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{gtf}:{line_number}: expected 9 GTF columns")
            meta = attributes(fields[8])
            symbol = requested_by_fold[meta.get("gene_name", match.group(1)).casefold()]
            gene_id = meta.get("gene_id")
            if not gene_id:
                raise ValueError(f"{gtf}:{line_number}: matched record has no gene_id")
            records[symbol].append((fields, meta, line))
            gene_ids[symbol].add(gene_id)

    missing = [symbol for symbol, rows in records.items() if not rows]
    if missing:
        raise ValueError(f"Gene symbol(s) not found in {gtf}: {', '.join(missing)}")
    ambiguous = {symbol: ids for symbol, ids in gene_ids.items() if len(ids) != 1}
    if ambiguous:
        raise ValueError(f"Ambiguous gene symbol(s): {ambiguous}; use a stable gene ID or curated annotation")

    result = []
    for symbol in requested:
        rows = records[symbol]
        gene_id = next(iter(gene_ids[symbol]))
        gene_rows = [(fields, meta, line) for fields, meta, line in rows if fields[2] == "gene"]
        transcript_rows = [(fields, meta, line) for fields, meta, line in rows if fields[2] == "transcript"]
        annotated_ids = {
            versioned(meta["transcript_id"], meta.get("transcript_version"))
            for _fields, meta, _line in transcript_rows
            if meta.get("transcript_id")
        }

        if selection == "all":
            selected_versioned = annotated_ids
        elif selection == "canonical":
            selected_versioned = {
                versioned(meta["transcript_id"], meta.get("transcript_version"))
                for _fields, meta, line in transcript_rows
                if 'tag "MANE_Select"' in line or 'tag "Ensembl_canonical"' in line
            }
            if not selected_versioned:
                raise ValueError(f"{symbol}: no MANE Select or Ensembl canonical transcript in {gtf}")
        else:
            requested_ids = table_transcripts[symbol]
            requested_bases = {item.split(".", 1)[0] for item in requested_ids}
            selected_versioned = {
                item for item in annotated_ids if item in requested_ids or item.split(".", 1)[0] in requested_bases
            }
            missing_ids = sorted(requested_bases - {item.split(".", 1)[0] for item in selected_versioned})
            if missing_ids:
                preview = ", ".join(missing_ids[:8])
                suffix = " ..." if len(missing_ids) > 8 else ""
                raise ValueError(f"{symbol}: {len(missing_ids)} table transcript(s) absent from GTF: {preview}{suffix}")

        selected_bases = {item.split(".", 1)[0] for item in selected_versioned}
        kept_lines = []
        kept_coordinates = []
        for fields, meta, line in rows:
            transcript_id = meta.get("transcript_id")
            if fields[2] == "gene" or transcript_id in selected_bases:
                kept_lines.append(line if line.endswith("\n") else line + "\n")
                kept_coordinates.append((int(fields[3]), int(fields[4])))

        if not selected_versioned:
            raise ValueError(f"{symbol}: selection produced no transcripts")
        locus_fields, locus_meta, _line = gene_rows[0] if gene_rows else rows[0]
        start = int(locus_fields[3]) if gene_rows else min(start for start, _end in kept_coordinates)
        end = int(locus_fields[4]) if gene_rows else max(end for _start, end in kept_coordinates)
        gene_id_version = versioned(gene_id, locus_meta.get("gene_version"))
        result.append(
            GeneSubset(
                symbol=symbol,
                gene_id=gene_id,
                gene_id_version=gene_id_version,
                chrom=locus_fields[0],
                start=start,
                end=end,
                strand=locus_fields[6],
                transcript_ids=sorted(selected_versioned),
                lines=kept_lines,
                source_record_count=len(rows),
            )
        )
    return result
