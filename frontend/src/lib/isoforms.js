export const ISOFORM_SORT_OPTIONS = [
  { value: 'annotation', label: 'Annotation order' },
  { value: 'sg_tpm', label: 'SG mean TPM' },
  { value: 'total_tpm', label: 'Total mean TPM' },
  { value: 'swish_log2fc', label: 'Swish log2FC' },
  { value: 'qvalue', label: 'q-value' },
  { value: 'transcript', label: 'Transcript name' },
];

export const DEFAULT_SORT_DIRECTIONS = {
  annotation: 'asc',
  sg_tpm: 'desc',
  total_tpm: 'desc',
  swish_log2fc: 'desc',
  qvalue: 'asc',
  transcript: 'asc',
};

export function callMatches(row, filter) {
  const call = row.call ?? row.statistics?.call ?? 'not_tested';
  if (filter === 'all') return true;
  if (filter === 'significant') return call === 'enriched' || call === 'depleted';
  return call === filter;
}

export function callLabel(call) {
  if (call === 'enriched') return 'SG enriched';
  if (call === 'depleted') return 'SG depleted';
  if (call === 'not_significant') return 'Not significant';
  return 'Not tested';
}

export function callCounts(rows) {
  const counts = { all: rows.length, significant: 0, enriched: 0, depleted: 0, not_significant: 0, not_tested: 0 };
  rows.forEach((row) => {
    const call = row.call ?? row.statistics?.call ?? 'not_tested';
    if (call === 'enriched' || call === 'depleted') counts.significant += 1;
    if (Object.hasOwn(counts, call)) counts[call] += 1;
  });
  return counts;
}

function rowValue(row, key) {
  if (key === 'annotation') return row.sourceIndex;
  if (key === 'sg_tpm') return row.expression?.conditions?.SG?.mean_tpm;
  if (key === 'total_tpm') return row.expression?.conditions?.Total?.mean_tpm;
  if (key === 'swish_log2fc') return row.enrichment?.swish_log2fc;
  if (key === 'qvalue') return row.enrichment?.qvalue;
  return row.transcript_name || row.transcript_id;
}

export function filterAndSortIsoforms(rows, filter, sortKey, direction) {
  return rows.filter((row) => callMatches(row, filter)).sort((a, b) => {
    const av = rowValue(a, sortKey);
    const bv = rowValue(b, sortKey);
    if (av === null || av === undefined || Number.isNaN(av)) return 1;
    if (bv === null || bv === undefined || Number.isNaN(bv)) return -1;
    const comparison = typeof av === 'string' ? av.localeCompare(bv) : av - bv;
    return comparison * (direction === 'asc' ? 1 : -1);
  });
}
