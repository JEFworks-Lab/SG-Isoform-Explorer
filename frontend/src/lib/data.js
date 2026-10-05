const baseUrl = import.meta.env.BASE_URL.endsWith('/')
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;

export const DATA_ROOT = `${baseUrl}data/khong2017`;

async function getJson(path) {
  const response = await fetch(path);
  if (!response.ok) {
    throw new Error(`Failed to load ${path}: ${response.status}`);
  }
  return response.json();
}

export async function loadCoreDataset() {
  const manifest = await getJson(`${DATA_ROOT}/manifest.json`);
  const [genes, transcripts, expression, enrichment, samples] = await Promise.all([
    getJson(`${DATA_ROOT}/${manifest.resources.genes}`),
    getJson(`${DATA_ROOT}/${manifest.resources.transcripts}`),
    getJson(`${DATA_ROOT}/${manifest.resources.expression}`),
    getJson(`${DATA_ROOT}/${manifest.resources.enrichment}`),
    getJson(`${DATA_ROOT}/${manifest.samples}`),
  ]);
  return { manifest, genes, transcripts, expression, enrichment, samples };
}

function fillPattern(pattern, geneId, sampleId) {
  return pattern.replace('{gene_id}', geneId).replace('{sample_id}', sampleId);
}

export async function loadGeneTracks(manifest, samples, geneId) {
  const results = await Promise.all(
    samples.map(async (sample) => {
      const coveragePath = fillPattern(manifest.resources.coverage_pattern, geneId, sample.sample_id);
      const junctionPath = fillPattern(manifest.resources.junction_pattern, geneId, sample.sample_id);
      const [coverage, junctions] = await Promise.all([
        getJson(`${DATA_ROOT}/${coveragePath}`),
        getJson(`${DATA_ROOT}/${junctionPath}`),
      ]);
      return { sample, coverage, junctions };
    }),
  );
  return results;
}

export function formatNumber(value, digits = 1) {
  if (value === null || value === undefined || Number.isNaN(value)) return '—';
  if (Math.abs(value) >= 1000) return Intl.NumberFormat('en', { notation: 'compact', maximumFractionDigits: 1 }).format(value);
  return Number(value).toFixed(digits).replace(/\.0$/, '');
}
