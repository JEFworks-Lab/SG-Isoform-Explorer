# Project context and handoff

This document is the public-safe handoff record for Browser-v2. It explains
what the browser represents, which decisions were made, what is reproducible,
and what a future maintainer should verify before extending or publishing it.
Machine-specific paths and private storage inventories are intentionally
omitted.

## Project purpose

The Stress Granule Isoform Explorer was built to examine isoform-aware evidence
for stress-granule RNA versus Total RNA in arsenite-stressed U-2 OS cells from
Khong *et al.* 2017. The current proof of concept is limited to CAPNS1 and
ANXA2, with the intention that the same pipeline can later support additional
genes and compatible RNA-seq datasets.

The browser is an evidence viewer. It places transcript annotation, Salmon
abundance, Fishpond/Swish statistics, read depth, and observed STAR splice
junctions in one interface. It does not perform alignment, transcript
quantification, or differential testing itself.

## Current data scope

| Item | Current value |
|---|---|
| Study | Khong *et al.* 2017, GEO GSE99304 |
| Organism | Human |
| Cell model | Arsenite-stressed U-2 OS cells |
| Contrast | Stress Granule RNA / Total RNA |
| Samples | Three SG and three Total RNA replicates |
| Genome | GRCh38.p14 |
| Annotation | Ensembl release 115 |
| Focus genes | CAPNS1 and ANXA2 |
| Transcript universe | Models present in the configured Swish metrics table |
| Current model counts | 11 CAPNS1 and 44 ANXA2 transcripts |

The 11/44 models are the table-selected analysis universe, not every transcript
in Ensembl 115. The original hand-curated browser contained only 4 CAPNS1 and
7 ANXA2 models. Both are subsets of the full Ensembl annotation, but they serve
different purposes and must not be described interchangeably.

## Evidence shown in the browser

### Transcript structures

Gene-specific GTF subsets provide transcript, exon, CDS, and UTR geometry.
Users choose genes by symbol; preparation resolves Ensembl IDs, coordinates,
strand, and transcript membership.

### Abundance

Salmon `quant.sf` supplies per-sample transcript TPM and estimated counts.
The abundance panel reports SG and Total means/replicates. These values are
model-based transcript estimates rather than direct counts of reads uniquely
assigned to a displayed isoform.

### Differential evidence

Fishpond/Swish supplies transcript-level effect sizes and inferential values.
The current display call is significant when both conditions hold:

```text
absolute Swish log2FC > 1
q-value < 0.01
```

Significant positive effects are labeled enriched and significant negative
effects depleted, relative to SG / Total. The thresholds are written into the
configuration and output so the frontend does not invent them.

### Read depth

Coverage is extracted from indexed, coordinate-sorted BAMs at exact 1-bp
resolution. The current filter requires `NH=1` and excludes unmapped,
secondary, QC-failed, duplicate, and supplementary alignments (`0xF04`). Raw
depth and depth per million STAR uniquely mapped reads are stored.

The 5-bp and 25-bp controls average the stored 1-bp positions for display.
They do not change which alignments were retained and do not affect junction
counts.

### Splice junctions

Junction arcs use final STAR `SJ.out.tab` rows. Browser exon-edge coordinates
are calculated as:

```text
donor = STAR column 2 - 1
acceptor = STAR column 3 + 1
```

Subset-matched arcs have endpoints matching adjacent exons in the displayed
transcript subset. All-final-STAR arcs include every final STAR event in the
locus, whether or not the selected GTF models explain it.

STAR rows with zero uniquely mapping reads and positive multimapping support
are retained as `multi_mapping_only`. They are hidden by the default minimum
unique-read filter and are visually distinguished when explicitly included.

Coverage and junctions are locus-level alignment evidence. They must not be
described as proving that a particular read came from one displayed isoform.

## Reproducible data flow

```text
Full Ensembl GTF ────────────────┐
Sample sheet ────────────────────┤
BAM + BAI ───────────────────────┤
STAR SJ.out.tab + Log.final.out ─┼─> preparation/build ─> browser JSON
Salmon quant.sf ─────────────────┤                         │
Swish/Fishpond CSV ──────────────┘                         ▼
                                                  independent validation
                                                           │
                                                           ▼
                                                  React/Vite static site
```

The human-edited inputs are `config/project.example.json` and
`config/samples.csv`. Preparation generates gene-specific GTFs and the detailed
`config/browser_build.generated.json`. The latter is an executable audit record
and should be regenerated rather than manually edited.

The browser-ready data live in `frontend/public/data/khong2017/`. That directory
is intentionally committed so a public Pages build is independent of private
raw-data storage.

## Important historical decisions

1. **Reference consistency:** annotation and coordinates are standardized on
   GRCh38 / Ensembl 115.
2. **Gene-symbol workflow:** collaborators should not have to know Ensembl IDs
   or type chromosome/locus metadata by hand.
3. **Explicit transcript selection:** full annotation, metrics-table models,
   and hand-curated models are distinct selectable universes.
4. **STAR is authoritative for junction counts:** BAM-derived depth is not used
   to invent arc counts.
5. **Exact coverage is preserved:** display binning happens after 1-bp source
   coverage is built.
6. **Descriptive and inferential values remain separate:** Salmon TPM/ratios do
   not replace Swish effect sizes and q-values.
7. **Validation re-reads sources:** it does not merely validate JSON syntax.
8. **Linked transcript navigation:** one filter, sort order, and focused model
   now coordinate abundance, Swish evidence, and transcript architecture.
9. **Sashimi annotation stays non-inferential:** each sample can show one
   independently selected isoform beneath its plot, but the overlay never
   assigns locus-level coverage or junction support to that transcript.

The abundance redesign uses compact vertical SG/Total bars with all three
replicate TPMs retained as labeled points. Mean/replicate switching changes the
bar value while preserving replicate context. These are frontend presentation
changes only; no browser JSON values or scientific thresholds were altered.

## Known caveats and unresolved release decisions

### Canonical statistics export

Two historical directories contain similarly named Swish/Fishpond exports.
Their focus-gene transcript membership and TPM/effect values agree, but their
inferential values, including q-values, differ. The current build uses the
export configured under the Salmon quasimapping analysis directory. Before a
scientific release, the data owner should identify the canonical export,
explain its provenance, and retain a checksum.

### Data redistribution

The static site commits derived JSON, including processed coverage and
junction evidence. Before public release, confirm that redistribution is
consistent with the original study's terms, relevant annotation licenses,
collaborator expectations, and institutional policy.

### Software and content licenses

No license has been chosen by this handoff. The owner should choose a software
license and separately consider whether the included or derived scientific
data need attribution or additional terms.

### Scope of interpretation

This is a focused browser, not a genome-wide differential expression analysis.
The descriptive summed SG/Total TPM ratio is not a gene-level hypothesis test.
Junctions and coverage support locus-level splice/read patterns but do not by
themselves identify transcript-of-origin.

## Handoff: routine tasks

### Run the current browser

```bash
cd frontend
npm ci
npm run dev
```

### Rebuild the current dataset

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python pipeline/prepare_project.py --project config/project.example.json
.venv/bin/python pipeline/build_dataset.py --config config/browser_build.generated.json
.venv/bin/python tests/validate_dataset.py --config config/browser_build.generated.json
cd frontend
npm ci
npm run build
```

### Add genes or adapt a dataset

1. copy and edit the project JSON;
2. update the sample sheet;
3. select a transcript-universe policy;
4. map nonstandard filenames and statistics columns;
5. prepare, build, validate, and visually review;
6. document changed scientific decisions in this file.

Detailed field descriptions are in `docs/PIPELINE_GUIDE.md`.

## Public repository boundary

Appropriate to commit:

- pipeline and frontend source code;
- schemas and tests;
- human-edited configuration examples;
- small generated GTF subsets;
- the reviewed browser-ready JSON;
- public-safe documentation and the deployment workflow.

Keep out of the public repository:

- BAM, BAI, FASTQ, CRAM, and private/full raw-data mirrors;
- `node_modules`, virtual environments, caches, and build directories;
- credentials or access tokens;
- internal absolute filesystem paths and machine-specific audit inventories;
- unpublished or restricted datasets.

The root `.gitignore` encodes these defaults. It does not replace a manual
review of `git status` and the staged diff before publishing.

## Release checklist

- [ ] Confirm the canonical metrics/Swish export with the data owner.
- [ ] Confirm permission to publish the derived JSON and annotation subset.
- [ ] Choose and add an appropriate license or licenses.
- [ ] Remove private paths, credentials, and unpublished information.
- [ ] Run source preparation, build, and independent validation.
- [ ] Run `npm ci` and `npm run build` in `frontend/`.
- [ ] Spot-check both genes, all sample tracks, and both junction scopes.
- [ ] Verify q-values/effects against the canonical source table.
- [ ] Set GitHub Pages Source to GitHub Actions.
- [ ] Verify the deployed site's favicon, JSON requests, and direct URL.
- [ ] Record the release tag, source checksums, and validation date.

## Where to look next

- `README.md`: shortest route to run and publish the browser.
- `docs/PIPELINE_GUIDE.md`: input contract and complete rebuild instructions.
- `config/project.example.json`: human-edited project intent.
- `config/samples.csv`: sample metadata.
- `config/browser_build.generated.json`: resolved machine-facing build record.
- `pipeline/`: preparation and source conversion code.
- `tests/validate_dataset.py`: independent data checks.
- `frontend/src/lib/data.js`: runtime JSON entry point.
- `frontend/public/data/khong2017/manifest.json`: deployed dataset index.
