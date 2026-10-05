# Pipeline guide: source RNA-seq files to browser JSON

This guide is for the person rebuilding the current dataset, adding genes, or
adapting the project to a new RNA-seq experiment. To run or publish the
already-built website, use the shorter instructions in the project README;
raw sequencing files are not required for that workflow.

## What the pipeline does

```text
gene symbols + sample sheet + source locations
                       │
                       ▼
              prepare_project.py
       resolves genes, transcripts, and files
                       │
                       ▼
        browser_build.generated.json
                       │
                       ▼
               build_dataset.py
  GTF + Salmon + Swish + BAM + STAR → browser JSON
                       │
                       ▼
             validate_dataset.py
   independently re-reads sources and checks outputs
                       │
                       ▼
       frontend/public/data/khong2017/
```

The frontend never parses BAM, BAI, GTF, Salmon, or STAR files. All expensive
processing happens before deployment.

## Prerequisites

- Python 3
- Node.js `^20.19.0` or `>=22.12.0`
- enough local storage and access to every configured source file
- indexed, coordinate-sorted BAMs

Create the Python environment once, from the repository root:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
```

Only `pysam` is currently required by the Python pipeline.

## The two files a new user edits

### 1. `config/project.example.json`

This file records scientific intent and file-discovery rules. For a new
project, copy it to a clearly named file and edit the copy:

```bash
cp config/project.example.json config/project.my_dataset.json
```

The important fields are:

| Field | Meaning |
|---|---|
| `dataset` | Public labels and reference description |
| `conditions` | Condition IDs and display labels |
| `comparison` | Numerator and denominator for the biological contrast |
| `reference_gtf` | Full annotation GTF matching the alignment reference |
| `genes` | Gene symbols to include; no Ensembl IDs are required |
| `transcript_selection` | Which annotated transcripts enter the browser |
| `sample_sheet` | CSV described below |
| `data_roots` | Directories containing alignments and quantification results |
| `file_patterns` | How a sample ID maps to its source files |
| `statistics` | Metrics CSV, column mapping, and significance rule |
| `output` | Regenerable pipeline output directory |
| `public_output` | JSON directory served by the frontend |

Paths are resolved relative to the project JSON file, not relative to the
shell's current directory.

### 2. `config/samples.csv`

One row represents one biological sample. Required columns are:

```csv
sample_id,condition,replicate,display_name
SRR5605163,SG,1,Stress Granule RNA · replicate 1
```

- `sample_id` must match the filename/folder patterns.
- `condition` must equal one of the IDs under `conditions`.
- `replicate` identifies the biological replicate within that condition.
- `display_name` is visitor-facing text.

Do not add source paths to every row when a consistent filename pattern can
derive them.

## Complete source-file contract

For every sample, the current adapter expects:

| Source | Current default pattern | Used for |
|---|---|---|
| BAM | `{sample_id}_Aligned.sortedByCoord.out.bam` | exact filtered read depth and sampled read evidence |
| BAI | `{sample_id}_Aligned.sortedByCoord.out.bam.bai` | indexed locus queries |
| STAR junctions | `{sample_id}_SJ.out.tab` | observed junction coordinates and counts |
| STAR log | `{sample_id}_Log.final.out` | uniquely mapped-read normalization denominator |
| Salmon quantification | `{sample_id}/quant.sf` | transcript TPM and estimated counts |

Dataset-level sources are:

- a full GTF using the same genome assembly and chromosome convention as the
  BAMs;
- optionally, a Fishpond/Swish CSV containing transcript IDs, effect sizes,
  q-values, and any selection fields.

Change `file_patterns` if collaborators use different filenames. Do not rename
large source files solely to satisfy the example configuration.

## Select transcripts deliberately

`transcript_selection.mode` controls the browser's transcript universe:

- `all`: every transcript annotated for each gene in the reference GTF;
- `table`: only transcripts for each gene that occur in the configured
  CSV/TSV, such as the transcripts tested by Swish;
- `canonical`: MANE Select and/or Ensembl canonical transcript models.

For the current Ensembl 115 sources:

| Universe | CAPNS1 | ANXA2 | Interpretation |
|---|---:|---:|---|
| Full GTF | 25 | 83 | all annotated transcripts |
| Current Swish table | 11 | 44 | tested/table-selected transcripts |
| Original curated subset | 4 | 7 | hand-selected display models |

The current Browser-v2 generated build uses the table-selected 11/44 universe.
That is not the complete Ensembl annotation, and the README/UI should not imply
that it is.

### Create gene-specific GTFs directly

All annotated models for gene symbols:

```bash
.venv/bin/python pipeline/subset_gtf_by_gene.py \
  --gtf /path/to/Homo_sapiens.GRCh38.115.gtf \
  --genes CAPNS1 ANXA2 \
  --selection all \
  --output-dir generated/gtf
```

Only models listed in an analysis table:

```bash
.venv/bin/python pipeline/subset_gtf_by_gene.py \
  --gtf /path/to/Homo_sapiens.GRCh38.115.gtf \
  --genes CAPNS1 ANXA2 \
  --selection table \
  --table /path/to/swish_results.csv \
  --gene-column external_gene_name \
  --transcript-column ensembl_transcript_id_version \
  --output-dir generated/gtf
```

The utility scans the full GTF, resolves each symbol, keeps complete feature
records, and writes `gene_subsets.summary.json` containing the inferred gene
IDs, locus, strand, transcript IDs, and record counts. It fails on missing or
ambiguous symbols rather than silently selecting the wrong gene.

## Run the pipeline

From the repository root:

```bash
.venv/bin/python pipeline/prepare_project.py \
  --project config/project.example.json

.venv/bin/python pipeline/build_dataset.py \
  --config config/browser_build.generated.json

.venv/bin/python tests/validate_dataset.py \
  --config config/browser_build.generated.json

cd frontend
npm ci
npm run build
```

For a copied project file, pass that file to `prepare_project.py`; it still
emits the generated build manifest named by the preparation script.

### What each command guarantees

1. **Prepare** resolves gene symbols, selects transcript models, generates
   gene-specific GTFs, discovers files from patterns, and checks that required
   inputs exist.
2. **Build** parses the selected GTF, Salmon, Swish, BAM, STAR junction, and
   STAR log inputs and atomically creates browser-ready JSON.
3. **Validate** independently re-reads the configured sources, recalculates
   important values, checks joins and coordinates, and verifies that the
   public copy matches the canonical output.
4. **Frontend build** confirms the React application compiles with its locked
   dependencies and deployable asset paths.

Preparation generates `config/browser_build.generated.json`. Treat that file
as a machine-facing record: edit the project JSON or sample sheet and rerun
preparation instead of hand-editing it.

## Do not mix build configurations

The build and validation commands are a matched pair. If you build with:

```text
config/browser_build.generated.json
```

then validate with that same file. The older `config/browser_build.json`
describes the frozen 4-CAPNS1/7-ANXA2 view. Validating a table-selected 11/44
build against the frozen manifest correctly fails at the transcript-set check.

## Source-to-JSON mapping

| Source | Browser output |
|---|---|
| GTF | `genes.json`, `transcripts.json` |
| Salmon `quant.sf` | `expression.json` |
| Swish/Fishpond metrics CSV | `enrichment.json` |
| BAM + BAI | `coverage/*.json`, `reads/*.json` |
| STAR `SJ.out.tab` | `junctions/*.json` |
| STAR `Log.final.out` | normalization fields in `samples.json` and coverage metadata |
| Dataset/build metadata | `manifest.json` |

`manifest.json` is the runtime index. The frontend first loads it, then uses
its resource names and patterns to find all remaining JSON files.

## Coverage: reads, positions, and display bins

Coverage is calculated at exact 1-bp resolution from the indexed BAM.
For each retained alignment, every aligned reference base contributes one unit
of depth at that position. A spliced read contributes to the aligned exon
blocks on both sides of its skipped intron; it does not fill the intron with
coverage.

The current filter keeps primary, uniquely mapping alignments (`NH=1`) and
excludes unmapped, secondary, QC-failed, duplicate, and supplementary records
with SAM mask `0xF04`.

The 5-bp and 25-bp controls are display summaries of the immutable 1-bp
coverage. Each displayed value is the mean depth of the source positions in
that window. They make a long/dense locus faster and easier to read; they do
not decide whether a read is counted and they do not change junction counts.

Both values are stored:

- raw depth: retained reads covering a position;
- normalized depth: raw depth per million STAR uniquely mapped reads.

## Junctions and the STAR coordinate conversion

Junction counts come from the final STAR `SJ.out.tab`, not from the coverage
bins. For this browser's 1-based inclusive exon-edge convention:

```text
donor = STAR column 2 - 1
acceptor = STAR column 3 + 1
```

This conversion places each arc on the adjacent exon bases. It does not alter
coverage or count reads. A one-base error would make a correct-count arc miss
the transcript's exon edge.

The browser supports two junction scopes:

- **subset-matched arcs**: final STAR junctions whose converted endpoints
  match adjacent exon edges in the displayed transcript subset;
- **all final STAR arcs**: all final STAR junction rows within the gene locus,
  including events not represented by the displayed annotation subset.

Neither option changes the underlying STAR counts.

### Zero uniquely mapping reads

STAR column 7 is uniquely mapping junction support and column 8 is
multimapping support. A retained row can therefore contain zero unique reads
but positive multimapping support. The pipeline preserves it as
`multi_mapping_only` evidence. The browser defaults to at least one unique
read; “Include 0 unique” shows these rows distinctly and never treats them as
positive unique-read support.

## Expression and differential statistics

Salmon TPM and Swish evidence answer different questions:

- `expression.json` stores per-sample transcript TPM and descriptive means;
- `enrichment.json` stores Swish effect sizes, p/q-values, availability, and
  the configured display call.

The current call is reproducibly encoded as:

```json
"significance": {
  "effect_field": "swish_log2fc",
  "qvalue_field": "qvalue",
  "absolute_log2fc_gt": 1.0,
  "qvalue_lt": 0.01
}
```

Both inequalities are strict. A transcript is:

- `enriched` when significant and Swish log2FC is positive;
- `depleted` when significant and Swish log2FC is negative;
- `not_significant` when metrics exist but thresholds are not both met;
- `not_tested` when required metrics are unavailable.

Positive/negative direction is relative to the configured numerator (currently
SG) and denominator (currently Total RNA).

## Frontend view contract

No additional source file or JSON schema is needed for the linked browser
controls. At runtime the frontend joins `transcripts.json`, `expression.json`,
and `enrichment.json` by exact versioned `transcript_id`.

- One shared filter, sort order, and focused transcript drives abundance,
  differential evidence, and the full isoform-model panel.
- Vertical abundance bars show the selected mean/replicate; the labeled points
  always come from the three per-sample TPM entries already stored in
  `expression.json`.
- The one-model strip under each sashimi sample uses exon coordinates from
  `transcripts.json` and the same locus scale as coverage and junctions.
- Selecting that strip is presentation-only. It must never filter BAM depth,
  STAR arcs, or junction counts, because those are locus-level evidence rather
  than transcript assignments.

The existing validator's exact transcript-set checks protect the joins. A
frontend change still requires `npm run build` and visual checks of CAPNS1,
ANXA2, intron compression, linked filtering/sorting, and independent sample
isoform selectors.

The descriptive SG/Total ratio is calculated from summed mean TPM with a small
pseudocount:

```text
(sum SG mean TPM + 0.01) / (sum Total mean TPM + 0.01)
```

It is not a gene-level significance test.

## Multiple metrics exports exist

The historical analysis folders contain two similarly named Swish/Fishpond
exports. Their focus-gene transcript membership and TPM/effect values agree,
but inferential values such as q-values differ. The current project file uses
the copy under `Khong2017/salmon_quasimapping/`.

Before public scientific release, choose the canonical export with the data
owner, document how it was produced, and retain a checksum. Do not silently
change `statistics.metrics_csv` to the other copy. A CSV containing only
significant/sequence-enriched transcripts must not be used as the full tested
transcript universe.

## Output and deployment

The pipeline writes the same dataset to two locations:

- `browser-data/`: regenerable staging/canonical pipeline output;
- `frontend/public/data/khong2017/`: served copy committed for the static site.

Only the second location is needed by GitHub Pages. It must remain committed,
because the Pages workflow does not have access to private BAMs or analysis
directories. The `.gitignore` excludes `browser-data/` to avoid committing two
identical copies.

The frontend uses Vite's base URL when fetching JSON and the favicon. The
included workflow sets that base to `/<repository>/` for a GitHub project site
and `/` for an `<owner>.github.io` repository.

## Validation checklist

Before accepting a rebuilt dataset:

- preparation reports every required source as present;
- the validator exits successfully with the same build config;
- transcript IDs in GTF, Salmon, and statistics sources join as expected;
- BAM reference names/build agree with the GTF;
- BAI queries succeed for every configured locus;
- junction endpoint conversion passes exon-edge checks;
- coverage and STAR values are independently reproduced;
- `browser-data/` and the public copy are byte-identical;
- `npm run build` succeeds;
- CAPNS1 and ANXA2 are manually spot-checked in the UI at 1-bp resolution.

An htslib warning that a BAI timestamp is older than its BAM is not by itself
proof of a bad index. If indexed locus queries and independent validation pass,
do not rewrite source data merely to remove the warning. Re-index when the BAM
contents change or when queries fail.

## Adapting the project for collaborators

For a new gene list or dataset:

1. copy the example project JSON;
2. enter gene symbols and reference/source locations;
3. update the sample sheet;
4. update filename patterns and statistics column mappings if necessary;
5. select the transcript universe explicitly;
6. prepare, build, and validate;
7. review the generated summary/report and UI;
8. commit code, configuration, small GTF subsets, and public JSON—not raw
   sequencing data or local absolute-path reports.

The public repository should explain where authorized maintainers obtain the
large inputs, but it should not expose private storage paths, credentials, or
unpublished samples.
