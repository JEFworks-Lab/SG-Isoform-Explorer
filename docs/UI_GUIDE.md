# Browser interface guide

This guide explains how the linked transcript views work and, just as
importantly, what the controls do **not** change.

## Page order

The browser presents evidence in this order:

1. **Abundance profile** — Salmon transcript TPM;
2. **Differential evidence** — Fishpond/Swish transcript statistics;
3. **Isoform models** — selected GTF transcript structures;
4. **Sashimi view** — locus-level BAM depth and final STAR junctions.

Isoform models are above the sashimi plots so the annotation vocabulary is
visible before the alignment evidence.

## Linked isoform view

The filter and sort controls in abundance, differential evidence, and isoform
models are synchronized. Changing a control in any one of those panels updates
all three panels.

Filters use the pipeline-authored transcript call:

- significant = enriched or depleted;
- SG enriched = significant positive Swish log2FC;
- SG depleted = significant negative Swish log2FC;
- not significant = metrics are present but the configured thresholds are not
  both met;
- not tested = the required inferential metrics are unavailable.

Sorting can use annotation order, SG or Total mean TPM, Swish log2FC, q-value,
or transcript name. Clicking or keyboard-selecting a transcript row focuses the
same transcript in abundance, differential evidence, and architecture. Focus,
filtering, and sorting are presentation state only: they do not recalculate TPM,
q-values, coverage, junction counts, or the descriptive all-isoform ratio.

## Abundance profile

Each transcript card has two vertical bars on one shared y-scale: Stress
Granule RNA and Total RNA. The bar is the active `Mean TPM`, `Replicate 1`,
`Replicate 2`, or `Replicate 3` value. Three labeled points remain overlaid at
the actual replicate TPM positions. When a replicate view is active, its point
is emphasized; the other points remain visible for context.

The SG/Total summary is descriptive. It uses all subset transcripts rather
than only the currently filtered cards and is not a gene-level significance
test.

## Differential evidence

Bars encode Swish log2FC, with negative values on the depleted side and
positive values on the enriched side. The q-value is written separately. The
call thresholds come from `enrichment.json` and are currently strict
`abs(Swish log2FC) > 1` and `q-value < 0.01`.

## Isoform models and sashimi tracks

The full isoform-model panel follows the linked filter and sort. Under every
visible sample sashimi plot is a separate selector that shows exactly one
transcript model at a time on the same genomic x-scale as that sample's depth
and junction arcs. Each sample selector is independent, which makes it possible
to compare different annotations under SG and Total tracks.

The one-isoform strip is an **annotation overlay**, not a read-assignment
filter. Selecting an isoform does not change BAM coverage, STAR arcs, arc
counts, or the junction table, and it does not claim that the displayed reads
came from that transcript. The linked call filter only limits which transcript
models are offered in the selector.

## Accessibility and keyboard use

- Filter buttons and sample toggles expose their selected state.
- Transcript cards, statistical rows, isoform models, arcs, and junction-table
  rows can be focused with the keyboard.
- Press Enter or Space on a transcript row to focus it across linked panels.
- Every sashimi isoform overlay uses a labeled native select control.
- Full junction coordinates and counts remain available in the table below
  each plot even when dense arc labels cannot all be drawn.

## Data contract

These interactions require no extra pipeline outputs. They join existing rows
by exact versioned `transcript_id` across:

- `transcripts.json` for exon geometry;
- `expression.json` for sample and condition TPM;
- `enrichment.json` for effect, q-value, and call;
- `coverage/*.json` and `junctions/*.json` for locus evidence.

`tests/validate_dataset.py` checks those transcript joins before release.
