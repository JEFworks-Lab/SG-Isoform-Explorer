# Stress Granule Isoform Explorer

An interactive, isoform-aware RNA-seq browser for the stress-granule (SG) and
Total RNA samples from Khong *et al.* 2017. The included dataset focuses on
**CAPNS1** and **ANXA2**, aligned to GRCh38 and annotated with Ensembl release
115.

The browser combines four related but distinct forms of evidence:

- transcript structures from a gene-specific GTF;
- transcript abundance from Salmon TPM;
- transcript-level differential evidence from Fishpond/Swish;
- locus-level read depth and splice junctions from BAM and STAR outputs.

The deployed website is fully static. It reads the prepared JSON files in
`frontend/public/data/khong2017/`; it does **not** open BAM, GTF, Salmon, or
STAR files in a visitor's browser.

Project background and scientific decisions are recorded in
[`docs/PROJECT_CONTEXT.md`](docs/PROJECT_CONTEXT.md). To build a new dataset or
add genes, use [`docs/PIPELINE_GUIDE.md`](docs/PIPELINE_GUIDE.md). The linked
filters, vertical abundance plots, and per-sample isoform overlays are explained
in [`docs/UI_GUIDE.md`](docs/UI_GUIDE.md).

## Choose the workflow you need

| Goal | What you need | What to run |
|---|---|---|
| View the included browser locally | Node.js and this repository | `npm ci`, then `npm run dev` in `frontend/` |
| Publish the included browser | A GitHub repository with Pages enabled | Push `main`; the included workflow builds and deploys it |
| Rebuild CAPNS1/ANXA2 from source | Python, `pysam`, GTF, BAM/BAI, STAR, Salmon, and Swish files | Prepare, build, and validate the dataset |
| Adapt the project | The same source types plus a sample sheet and gene symbols | Edit `config/project.example.json` and `config/samples.csv`, then run the pipeline |

The first two workflows do not require access to the raw sequencing files.

## Run the included browser locally

Vite 8 requires Node.js `^20.19.0` or `>=22.12.0`. From the repository root:

```bash
cd frontend
npm ci
npm run dev
```

Open the local URL printed by Vite. If port 5173 is already occupied, Vite
will select another port; that is normal.

To verify the production build:

```bash
cd frontend
npm ci
npm run build
npm run preview
```

`npm ci` uses the committed lockfile, so it is preferred over `npm install`
for a reproducible checkout.

## Publish with GitHub Pages

This repository includes `.github/workflows/deploy-pages.yml`. The workflow:

1. installs the locked frontend dependencies;
2. selects the correct Vite base path for either a project site or a
   `<owner>.github.io` site;
3. builds `frontend/dist/`;
4. uploads and deploys the static artifact to GitHub Pages.

After creating a GitHub repository:

```bash
git init
git add .
git commit -m "Initial public browser"
git branch -M main
git remote add origin https://github.com/<owner>/<repository>.git
git push -u origin main
```

Replace `<owner>` and `<repository>` with the actual GitHub values; do not
copy the angle-bracket placeholders literally. Then, on GitHub:

1. open **Settings → Pages**;
2. set **Source** to **GitHub Actions**;
3. open the **Actions** tab and wait for “Deploy browser to GitHub Pages.”

For a normal project repository, the site URL is usually:

```text
https://<owner>.github.io/<repository>/
```

No source BAMs or private filesystem paths are needed by the Pages build.

## Files required to deploy the existing site

Commit these paths:

| Path | Purpose |
|---|---|
| `frontend/src/` | React components, data loading, and styles |
| `frontend/public/data/khong2017/` | Browser-ready dataset served to visitors |
| `frontend/public/dna-favicon.png` | Site icon |
| `frontend/index.html` | HTML entry point |
| `frontend/package.json` and `package-lock.json` | JavaScript dependencies and commands |
| `frontend/vite.config.js` | Local and GitHub Pages asset base path |
| `.github/workflows/deploy-pages.yml` | Automated Pages build and deployment |

The JSON directory must include `manifest.json`, `samples.json`,
`genes.json`, `transcripts.json`, `expression.json`, `enrichment.json`, and
the per-gene/per-sample files below `coverage/`, `junctions/`, and `reads/`.
`manifest.json` is the index that tells the frontend where to find them.

## Rebuild the included biological data

The pipeline is optional for website visitors but required when the data or
gene set changes. From the repository root:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt

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

Always validate with the same configuration used to build. The current
gene-symbol workflow generates `config/browser_build.generated.json`; do not
validate that output with the older frozen `config/browser_build.json`.

The build writes:

- `browser-data/`: regenerable pipeline staging output;
- `frontend/public/data/khong2017/`: the identical, deployable copy.

`browser-data/` is intentionally ignored by Git because the served copy is
the version needed for GitHub Pages.

## Inputs required to rebuild or adapt the dataset

| Input | Required information |
|---|---|
| Reference GTF | Full annotation matching the alignments; current build uses Ensembl 115 / GRCh38 |
| Sample sheet | Sample ID, condition, replicate, and display label |
| BAM + BAI | Coordinate-sorted alignment and matching index for every sample |
| STAR `SJ.out.tab` | Final observed splice-junction table for every sample |
| STAR `Log.final.out` | Uniquely mapped-read count used for coverage normalization |
| Salmon `quant.sf` | Transcript TPM and estimated counts for every sample |
| Swish/Fishpond CSV | Optional transcript effect sizes, q-values, and selection table |

Users specify gene symbols such as `CAPNS1`; the preparation script derives
gene IDs, chromosome, locus, strand, transcript IDs, and file paths. See the
pipeline guide for naming patterns and column mappings.

## Repository map

```text
Browser/
├── .github/workflows/       GitHub Pages deployment
├── config/                  Human-edited and generated pipeline configuration
├── docs/                    User guide and public-safe scientific handoff
├── frontend/                React/Vite application and deployable JSON
├── generated/gtf/           Small gene-specific annotation subsets
├── pipeline/                Preparation and source-to-JSON conversion scripts
│   └── schemas/             JSON schemas for browser outputs
├── tests/                   Independent source/output validation
└── requirements.txt         Pinned Python dependency
```

## Scientific conventions that must not change silently

- Browser coordinates are 1-based inclusive.
- Coverage starts as exact 1-bp BAM depth; 5-bp and 25-bp choices are display
  averages of that source, not different read counting.
- The alignment filter keeps `NH=1` and excludes SAM flags `0xF04`.
- Junction arcs use final STAR counts with
  `donor = STAR column 2 - 1` and `acceptor = STAR column 3 + 1`.
- BAM coverage and junctions are locus-level evidence; they do not assign a
  read to one transcript model.
- Salmon TPM is descriptive abundance. Swish log2 fold change and q-value are
  inferential statistics and are displayed separately.
- Isoform filter/sort/focus state is linked across abundance, differential
  evidence, and transcript architecture; it never rewrites source values.
- The one-isoform model under each sashimi track is an annotation overlay, not
  a BAM/STAR read-to-transcript assignment.
- The current significance rule is `abs(Swish log2FC) > 1` and `q < 0.01`.

## Source study

Khong *et al.* (2017), “The Stress Granule Transcriptome Reveals Principles
of mRNA Accumulation in Stress Granules”:
[paper](https://www.cell.com/molecular-cell/fulltext/S1097-2765(17)30790-6) ·
[GEO GSE99304](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99304)

Deployment references:
[Vite static deployment](https://vite.dev/guide/static-deploy.html) ·
[GitHub Pages custom workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages) ·
[GitHub Pages publishing source](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site)
