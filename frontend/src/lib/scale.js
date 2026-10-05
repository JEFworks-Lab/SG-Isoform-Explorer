export const PLOT_WIDTH = 1100;
export const PLOT_LEFT = 148;
export const PLOT_RIGHT = 28;

export function locusX(value, start, end) {
  return PLOT_LEFT + ((value - start) / (end - start || 1)) * (PLOT_WIDTH - PLOT_LEFT - PLOT_RIGHT);
}

export function locusTicks(start, end) {
  return Array.from({ length: 5 }, (_, index) => start + (end - start) * index / 4);
}

// Shared piecewise coordinate transform for both sashimi and transcript models.
export function createLocusScale(gene, transcripts, start = gene.start, end = gene.end, compact = false) {
  const clipped = transcripts.flatMap((tx) => tx.exons)
    .map((exon) => [Math.max(start, exon.start), Math.min(end, exon.end)])
    .filter(([a, b]) => a <= b).sort((a, b) => a[0] - b[0]);
  const exons = [];
  for (const [a, b] of clipped) {
    if (exons.length && a <= exons[exons.length - 1][1] + 1) exons[exons.length - 1][1] = Math.max(b, exons[exons.length - 1][1]);
    else exons.push([a, b]);
  }
  const stops = [...new Set([start, end + 1, ...exons.flatMap(([a, b]) => [a, b + 1])])].sort((a, b) => a - b);
  const segments = [];
  let cumulative = 0;
  for (let i = 0; i < stops.length - 1; i++) {
    const a = stops[i], b = stops[i + 1];
    const isIntron = !exons.some(([lo, hi]) => lo <= a && hi >= b - 1);
    const units = compact && isIntron ? Math.min(b - a, 85) : b - a;
    segments.push({ a, b, compressed: units < b - a, units, offset: cumulative });
    cumulative += units;
  }
  const total = cumulative || 1;
  function x(coordinate) {
    const value = Math.max(start, Math.min(end + 1, coordinate));
    const part = segments.find((s) => value < s.b) || segments[segments.length - 1];
    return PLOT_LEFT + (part ? part.offset + (value - part.a) * part.units / (part.b - part.a) : 0) / total * (PLOT_WIDTH - PLOT_LEFT - PLOT_RIGHT);
  }
  function invert(pixel) {
    const distance = Math.max(0, Math.min(total, (pixel - PLOT_LEFT) / (PLOT_WIDTH - PLOT_LEFT - PLOT_RIGHT) * total));
    const part = segments.find((s) => distance < s.offset + s.units) || segments[segments.length - 1];
    return part ? Math.min(end, Math.max(start, Math.round(part.a + (distance - part.offset) * (part.b - part.a) / part.units))) : start;
  }
  const candidates = compact ? [start, ...exons.flatMap(([a, b]) => [a, b]), end] : locusTicks(start, end).map(Math.round);
  const ticks = [];
  for (const t of candidates) if (!ticks.length || x(t) - x(ticks[ticks.length - 1]) >= 125) ticks.push(t);
  if (x(end) - x(ticks[ticks.length - 1]) >= 90) ticks.push(end);
  return { x, invert, ticks, start, end, compressed: segments.filter((s) => s.compressed).map((s) => (x(s.a) + x(s.b)) / 2) };
}
