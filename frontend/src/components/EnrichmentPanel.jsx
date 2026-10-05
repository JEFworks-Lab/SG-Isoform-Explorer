import { formatNumber } from '../lib/data.js';
import IsoformControls from './IsoformControls.jsx';
import { callLabel } from '../lib/isoforms.js';

function formatQValue(value) {
  const qValue = Number(value);
  if (!Number.isFinite(qValue)) return '—';
  if (qValue === 0) return '0';
  if (Math.abs(qValue) < 0.001) return qValue.toExponential(5).replace('e-', 'e−');
  return qValue.toPrecision(6);
}

function keyboardFocus(event, transcriptId, onFocusTranscript) {
  if (event.key !== 'Enter' && event.key !== ' ') return;
  event.preventDefault();
  onFocusTranscript(transcriptId);
}

export default function EnrichmentPanel({ isoformRows, visibleRows, callDefinition, controls,
  focusedTranscriptId, onFocusTranscript }) {
  const testedCount = isoformRows.filter((row) => row.enrichment?.metrics_available !== false).length;
  const maxAbs = Math.max(2.5, ...isoformRows.map((row) => Math.abs(row.enrichment?.swish_log2fc ?? 0)));
  const effectThreshold = callDefinition?.absolute_log2fc_gt ?? 1;
  const qThreshold = callDefinition?.qvalue_lt ?? 0.01;

  return <section className="panel enrichment-panel">
    <div className="section-heading">
      <div><span className="eyebrow">Inferential transcript statistics</span><h2>Differential evidence</h2></div>
      <span className="unit-chip">Fishpond / Swish</span>
    </div>
    <p className="section-intro">Swish log2 fold changes and multiple-testing-adjusted q-values. This panel follows the same isoform filter, order, and focused row as transcript architecture.</p>
    <div className="stats-summary" aria-label="Swish call summary">
      <div><span>Tested isoforms</span><strong>{testedCount}</strong><small>with Swish metrics</small></div>
      <div><span>Significant</span><strong>{controls.counts.significant}</strong><small>|log₂FC| &gt; {effectThreshold} · q &lt; {qThreshold}</small></div>
      <div className="is-positive"><span>SG enriched</span><strong>{controls.counts.enriched}</strong><small>significant positive effect</small></div>
      <div className="is-negative"><span>SG depleted</span><strong>{controls.counts.depleted}</strong><small>significant negative effect</small></div>
    </div>
    <IsoformControls {...controls} />
    <div className="enrichment-axis"><span aria-hidden="true" /><div className="enrichment-axis__effect"><span>SG depleted</span><i /><span>SG enriched</span></div><b>q-value</b></div>
    <div className="enrichment-list">
      {visibleRows.map((isoform) => {
        const row = isoform.enrichment ?? {};
        const value = row.swish_log2fc ?? 0;
        const width = (Math.abs(value) / maxAbs) * 40;
        const valueStyle = value >= 0 ? { left: `calc(50% + ${width}% + 8px)` } : { right: `calc(50% + ${width}% + 8px)` };
        return <div className={`enrichment-row call-row-${isoform.call} ${focusedTranscriptId === isoform.transcript_id ? 'is-focused' : ''}`}
          key={isoform.transcript_id} role="button" tabIndex="0" aria-pressed={focusedTranscriptId === isoform.transcript_id}
          onClick={() => onFocusTranscript(isoform.transcript_id)}
          onKeyDown={(event) => keyboardFocus(event, isoform.transcript_id, onFocusTranscript)}>
          <div className="enrichment-id"><b>{isoform.transcript_name || isoform.transcript_id}</b><small>{isoform.transcript_id}</small>
            <span className={`call-chip call-${isoform.call}`}>{callLabel(isoform.call)}</span></div>
          <div className="diverging-bar" title={`${isoform.transcript_name || isoform.transcript_id}: Swish log2FC ${value}; q-value ${row.qvalue}`}>
            <span className="center-line" />
            <i className={value >= 0 ? 'bar-positive' : 'bar-negative'} style={{ width: `${width}%`, [value >= 0 ? 'left' : 'right']: '50%' }} />
            <strong className={value >= 0 ? 'effect-value value-positive' : 'effect-value value-negative'} style={valueStyle}>
              {value >= 0 ? '+' : ''}{formatNumber(value, 3)}
            </strong>
          </div>
          <span className="q-value" aria-label={`q-value ${formatQValue(row.qvalue)}`} title={String(row.qvalue ?? 'Unavailable')}>q {formatQValue(row.qvalue)}</span>
        </div>;
      })}
      {!visibleRows.length && <div className="empty-filter-state"><b>No isoforms match this Swish call.</b><span>Choose another linked filter to restore transcript rows.</span></div>}
    </div>
    <p className="panel-note">Calls are transcript-level: <code>|Swish log2FC| &gt; {effectThreshold}</code> and <code>q-value &lt; {qThreshold}</code>. Selecting a row focuses the same isoform model and abundance card; it does not filter locus-level BAM or STAR evidence.</p>
  </section>;
}
