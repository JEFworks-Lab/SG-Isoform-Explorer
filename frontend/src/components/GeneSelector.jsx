import { useMemo, useState } from 'react';

export default function GeneSelector({ genes, selectedGeneId, onChange }) {
  const [query, setQuery] = useState('');
  const selected = genes.find((gene) => gene.gene_id === selectedGeneId);
  const visible = useMemo(() => {
    const q = query.trim().toLowerCase();
    if (!q) return genes;
    return genes.filter((gene) =>
      [gene.gene_id, gene.gene_name, gene.chrom].some((value) => String(value).toLowerCase().includes(q)),
    );
  }, [genes, query]);

  return (
    <div className="gene-picker panel-subtle">
      <div className="gene-picker__current">
        <span className="eyebrow">Selected locus</span>
        <div className="gene-picker__identity">
          <div>
            <strong>{selected?.gene_name || '—'}</strong>
            <span>{selected?.gene_id}</span>
          </div>
          <span className="pill">{selected?.strand} strand</span>
        </div>
      </div>
      <label className="gene-picker__search">
        <span>⌕</span>
        <input
          aria-label="Search genes"
          placeholder="Search gene name or ID…"
          value={query}
          onChange={(event) => setQuery(event.target.value)}
        />
      </label>
      <div className="gene-picker__list" role="listbox" aria-label="Genes">
        {visible.map((gene) => (
          <button
            key={gene.gene_id}
            className={`gene-option ${gene.gene_id === selectedGeneId ? 'is-active' : ''}`}
            onClick={() => {
              onChange(gene.gene_id);
              setQuery('');
            }}
          >
            <span>
              <b>{gene.gene_name}</b>
              <small>{gene.gene_id}</small>
            </span>
            <span className="gene-option__meta">{gene.transcript_ids.length} {gene.transcript_ids.length === 1 ? 'isoform' : 'isoforms'}</span>
          </button>
        ))}
      </div>
    </div>
  );
}
