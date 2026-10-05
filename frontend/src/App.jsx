import { useEffect, useMemo, useState } from 'react';
import GeneSelector from './components/GeneSelector.jsx';
import TranscriptTrack from './components/TranscriptTrack.jsx';
import ExpressionPanel from './components/ExpressionPanel.jsx';
import EnrichmentPanel from './components/EnrichmentPanel.jsx';
import SashimiPlot from './components/SashimiPlot.jsx';
import LoadingSkeleton from './components/LoadingSkeleton.jsx';
import { callCounts, DEFAULT_SORT_DIRECTIONS, filterAndSortIsoforms, ISOFORM_SORT_OPTIONS } from './lib/isoforms.js';
import { loadCoreDataset, loadGeneTracks } from './lib/data.js';
import { createLocusScale } from './lib/scale.js';

export default function App() {
  const [dataset, setDataset] = useState(null);
  const [selectedGeneId, setSelectedGeneId] = useState('ENSG00000126247');
  const [trackRows, setTrackRows] = useState([]);
  const [loading, setLoading] = useState(true);
  const [trackLoading, setTrackLoading] = useState(false);
  const [error, setError] = useState('');
  const [compact, setCompact] = useState(false);
  const [zoom, setZoom] = useState(null);
  const [startInput, setStartInput] = useState('');
  const [endInput, setEndInput] = useState('');
  const [zoomError, setZoomError] = useState('');
  const [isoformFilter, setIsoformFilter] = useState('all');
  const [isoformSort, setIsoformSort] = useState('annotation');
  const [isoformDirection, setIsoformDirection] = useState('asc');
  const [focusedTranscriptId, setFocusedTranscriptId] = useState(null);

  useEffect(() => {
    loadCoreDataset()
      .then((data) => {
        setDataset(data);
        if (!data.genes.some((gene) => gene.gene_id === selectedGeneId)) setSelectedGeneId(data.genes[0]?.gene_id);
      })
      .catch((err) => setError(err.message))
      .finally(() => setLoading(false));
  }, []);

  useEffect(() => {
    if (!dataset || !selectedGeneId) return;
    let ignore = false;
    setTrackLoading(true);
    setTrackRows([]);
    loadGeneTracks(dataset.manifest, dataset.samples, selectedGeneId)
      .then((rows) => { if (!ignore) setTrackRows(rows); })
      .catch((err) => { if (!ignore) setError(err.message); })
      .finally(() => { if (!ignore) setTrackLoading(false); });
    return () => { ignore = true; };
  }, [dataset, selectedGeneId]);

  const view = useMemo(() => {
    if (!dataset) return null;
    const gene = dataset.genes.find((item) => item.gene_id === selectedGeneId);
    if (!gene) return null;
    return {
      gene,
      transcripts: dataset.transcripts.filter((tx) => tx.gene_id === selectedGeneId),
      expressionRows: dataset.expression.transcripts.filter((row) => row.gene_id === selectedGeneId),
      geneEnrichment: dataset.enrichment.genes.find((row) => row.gene_id === selectedGeneId),
      transcriptEnrichment: dataset.enrichment.transcripts.filter((row) => row.gene_id === selectedGeneId),
    };
  }, [dataset, selectedGeneId]);
  const scale = useMemo(() => view && createLocusScale(view.gene, view.transcripts,
    zoom?.geneId === selectedGeneId ? zoom.start : view.gene.start,
    zoom?.geneId === selectedGeneId ? zoom.end : view.gene.end, compact),
    [view, zoom, selectedGeneId, compact]);
  const isoformRows = useMemo(() => {
    if (!view) return [];
    const expressionById = new Map(view.expressionRows.map((row) => [row.transcript_id, row]));
    const enrichmentById = new Map(view.transcriptEnrichment.map((row) => [row.transcript_id, row]));
    return view.transcripts.map((transcript, sourceIndex) => {
      const expression = expressionById.get(transcript.transcript_id);
      const enrichment = enrichmentById.get(transcript.transcript_id);
      return {
        ...transcript,
        expression,
        enrichment,
        call: enrichment?.call ?? expression?.statistics?.call ?? 'not_tested',
        sourceIndex,
      };
    });
  }, [view]);
  const visibleIsoforms = useMemo(() => filterAndSortIsoforms(
    isoformRows, isoformFilter, isoformSort, isoformDirection,
  ), [isoformRows, isoformFilter, isoformSort, isoformDirection]);
  const isoformCounts = useMemo(() => callCounts(isoformRows), [isoformRows]);

  useEffect(() => {
    if (!visibleIsoforms.length) {
      setFocusedTranscriptId(null);
      return;
    }
    if (!visibleIsoforms.some((row) => row.transcript_id === focusedTranscriptId)) {
      setFocusedTranscriptId(visibleIsoforms[0].transcript_id);
    }
  }, [visibleIsoforms, focusedTranscriptId]);

  if (loading) return <main className="app-shell"><LoadingSkeleton /></main>;
  if (error) return <main className="app-shell"><div className="error-card"><b>Dataset could not be loaded.</b><span>{error}</span></div></main>;
  if (!view) return null;

  function chooseGene(id) {
    setSelectedGeneId(id);
    setZoom(null);
    setZoomError('');
    setStartInput('');
    setEndInput('');
    setFocusedTranscriptId(null);
  }
  function changeIsoformSort(key) {
    setIsoformSort(key);
    setIsoformDirection(DEFAULT_SORT_DIRECTIONS[key]);
  }
  function applyZoom(event) {
    event.preventDefault();
    const start = Number(startInput.replaceAll(',', ''));
    const end = Number(endInput.replaceAll(',', ''));
    if (!Number.isInteger(start) || !Number.isInteger(end) || start < view.gene.start || end > view.gene.end || start >= end) {
      setZoomError(`Enter whole-base coordinates within ${view.gene.start.toLocaleString()}–${view.gene.end.toLocaleString()}, with start before end.`);
      return;
    }
    setZoomError('');
    setZoom({ geneId: selectedGeneId, start, end });
  }
  return <main className="app-shell">
    <header className="app-header">
      <div className="brand-lockup"><div className="brand-mark"><i /><i /><i /></div>
        <div><span>SG // ISOFORM EXPLORER</span><small>{dataset.manifest.dataset.name}</small></div></div>
      <div className="status-cluster"><span className="status-live"><i /> source validated</span><span>{dataset.manifest.dataset.genome}</span><span>{dataset.manifest.dataset.annotation}</span></div>
    </header>
    <div className="hero-grid">
      <section className="hero-copy">
        <span className="kicker">Khong 2017 · stress-granule fraction versus Total RNA</span>
        <h1>{view.gene.gene_name}<span> / {view.gene.gene_id}</span></h1>
        <p>Compare Ensembl 115 subset models with Salmon abundance, Fishpond/Swish statistics, exact uniquely mapped BAM depth, and final STAR splice-junctions.</p>
        <div className="locus-stat-row">
          <div><span>Locus</span><b>{view.gene.chrom}:{view.gene.start.toLocaleString()}–{view.gene.end.toLocaleString()}</b></div>
          <div><span>Subset models</span><b>{view.transcripts.length}</b></div>
          <div><span>Strand</span><b>{view.gene.strand}</b></div>
        </div>
      </section>
      <GeneSelector genes={dataset.genes} selectedGeneId={selectedGeneId} onChange={chooseGene} />
    </div>

    <div className="provenance-strip" aria-label="Data provenance">
      <span><i className="source-dot source-annotation" /> Ensembl 115 subset GTF</span>
      <span><i className="source-dot source-expression" /> Salmon 1.11.4 + Swish</span>
      <span><i className="source-dot source-alignment" /> STAR 2.7.11b final alignments</span>
      <span>Coordinates: 1-based inclusive · no <code>chr</code> prefix</span>
    </div>

    <div className="section-divider"><span>01 / EXPRESSION & DIFFERENTIAL EVIDENCE</span><p>One linked isoform view keeps abundance, annotation, and statistics in the same order.</p></div>
    <ExpressionPanel gene={view.gene} isoformRows={isoformRows} visibleRows={visibleIsoforms}
      geneEnrichment={view.geneEnrichment} samples={dataset.samples} controls={{
        filter: isoformFilter, onFilter: setIsoformFilter, counts: isoformCounts,
        sortKey: isoformSort, onSortKey: changeIsoformSort, direction: isoformDirection,
        onDirection: () => setIsoformDirection((current) => current === 'asc' ? 'desc' : 'asc'),
        sortOptions: ISOFORM_SORT_OPTIONS, visibleCount: visibleIsoforms.length,
      }} focusedTranscriptId={focusedTranscriptId} onFocusTranscript={setFocusedTranscriptId} />
    <EnrichmentPanel isoformRows={isoformRows} visibleRows={visibleIsoforms}
      callDefinition={dataset.enrichment.call_definition} controls={{
        filter: isoformFilter, onFilter: setIsoformFilter, counts: isoformCounts,
        sortKey: isoformSort, onSortKey: changeIsoformSort, direction: isoformDirection,
        onDirection: () => setIsoformDirection((current) => current === 'asc' ? 'desc' : 'asc'),
        sortOptions: ISOFORM_SORT_OPTIONS, visibleCount: visibleIsoforms.length,
      }} focusedTranscriptId={focusedTranscriptId} onFocusTranscript={setFocusedTranscriptId} />

    <div className="section-divider"><span>02 / GENOMIC EVIDENCE</span><p>Sample-specific read depth and splice support aligned with transcript architecture.</p></div>
    <div className="genomic-controls panel-subtle">
      <div className="genomic-controls__heading"><b>Genomic coordinates</b><span>{view.gene.chrom}:{scale.start.toLocaleString()}–{scale.end.toLocaleString()} · 1-based inclusive</span></div>
      <label className="compact-control"><input type="checkbox" checked={compact} onChange={(e) => setCompact(e.target.checked)} /> Compress introns</label>
      <form className="zoom-form" onSubmit={applyZoom}>
        <label>Start <input type="text" inputMode="numeric" placeholder={String(view.gene.start)} value={startInput} onChange={(e) => setStartInput(e.target.value)} /></label>
        <label>End <input type="text" inputMode="numeric" placeholder={String(view.gene.end)} value={endInput} onChange={(e) => setEndInput(e.target.value)} /></label>
        <button type="submit">Zoom to region</button>
        <button type="button" onClick={() => { setZoom(null); setZoomError(''); setStartInput(''); setEndInput(''); }}>Reset</button>
      </form>
      {zoomError && <span className="zoom-error" role="alert">{zoomError}</span>}
      {compact && <span className="genomic-caution">Introns are compressed on both plots; genomic distances are not to scale. Break marks show shortened introns.</span>}
    </div>
    <TranscriptTrack gene={view.gene} isoformRows={isoformRows} visibleRows={visibleIsoforms} locusScale={scale}
      controls={{
        filter: isoformFilter, onFilter: setIsoformFilter, counts: isoformCounts,
        sortKey: isoformSort, onSortKey: changeIsoformSort, direction: isoformDirection,
        onDirection: () => setIsoformDirection((current) => current === 'asc' ? 'desc' : 'asc'),
        sortOptions: ISOFORM_SORT_OPTIONS, visibleCount: visibleIsoforms.length,
      }} focusedTranscriptId={focusedTranscriptId} onFocusTranscript={setFocusedTranscriptId} />
    {trackLoading ? <LoadingSkeleton /> : <SashimiPlot trackRows={trackRows.filter((row) => row.coverage.gene_id === selectedGeneId)} samples={dataset.samples}
      scale={scale} transcripts={visibleIsoforms} focusedTranscriptId={focusedTranscriptId} />}
    <footer className="app-footer"><span>Khong <i>et al.</i> 2017 · U-2 OS stress-granule RNA versus Total RNA</span><span>GRCh38.p14 · Ensembl 115 · schema v{dataset.manifest.schema_version}</span></footer>
  </main>;
}
