import { useMemo, useState } from 'react';
import { formatNumber } from '../lib/data.js';
import IsoformControls from './IsoformControls.jsx';
import { callLabel } from '../lib/isoforms.js';

function keyboardFocus(event, transcriptId, onFocusTranscript) {
  if (event.key !== 'Enter' && event.key !== ' ') return;
  event.preventDefault();
  onFocusTranscript(transcriptId);
}

export default function ExpressionPanel({ gene, isoformRows, visibleRows, geneEnrichment, samples, controls,
  focusedTranscriptId, onFocusTranscript }) {
  const [view, setView] = useState('mean');
  const byCondition = useMemo(() => Object.fromEntries(['SG', 'Total'].map((condition) => [condition,
    samples.filter((sample) => sample.condition === condition)
      .sort((a, b) => Number(a.replicate) - Number(b.replicate))])), [samples]);
  const getValue = (row, condition) => view === 'mean'
    ? row.expression?.conditions?.[condition]?.mean_tpm ?? 0
    : row.expression?.samples?.[byCondition[condition].find((sample) => String(sample.replicate) === view)?.sample_id]?.tpm ?? 0;
  const totals = Object.fromEntries(['SG', 'Total'].map((condition) => [condition,
    isoformRows.reduce((sum, row) => sum + getValue(row, condition), 0)]));
  const maxGene = Math.max(1, totals.SG, totals.Total);
  const allVisibleValues = visibleRows.flatMap((row) => ['SG', 'Total'].flatMap((condition) => [
    getValue(row, condition),
    ...byCondition[condition].map((sample) => row.expression?.samples?.[sample.sample_id]?.tpm ?? 0),
  ]));
  const maxTx = Math.max(1, ...allVisibleValues);
  const ratio = view === 'mean' && geneEnrichment?.ratio
    ? geneEnrichment.ratio
    : (totals.SG + 0.01) / (totals.Total + 0.01);
  const log2Ratio = view === 'mean' && geneEnrichment?.log2_ratio !== undefined
    ? geneEnrichment.log2_ratio
    : Math.log2(ratio);

  return <section className="panel expression-panel">
    <div className="section-heading">
      <div><span className="eyebrow">Expression · transcript TPM</span><h2>Abundance profile</h2></div>
      <div className="segmented-control" role="group" aria-label="Expression sample view">
        {['mean', '1', '2', '3'].map((option) => <button key={option} type="button" className={view === option ? 'is-active' : ''}
          aria-pressed={view === option} onClick={() => setView(option)}>{option === 'mean' ? 'Mean TPM' : 'Replicate ' + option}</button>)}
      </div>
    </div>
    <p className="section-intro">Bars show the selected mean or replicate on one shared TPM scale. Each bar retains all three replicate measurements as labeled points.</p>
    <div className="condition-summary condition-summary--with-ratio">
      {['SG', 'Total'].map((condition) => <div className={'condition-card ' + (condition === 'SG' ? 'is-sg' : '')} key={condition}>
        <span>{condition === 'SG' ? 'Stress Granule RNA' : 'Total RNA'} · summed isoform TPM</span>
        <strong>{formatNumber(totals[condition], 1)} <small>TPM</small></strong>
        <div className="meter"><i style={{ width: totals[condition] / maxGene * 100 + '%' }} /></div>
      </div>)}
      <div className="condition-card ratio-card">
        <span>Descriptive SG / Total · {view === 'mean' ? 'mean TPM' : `replicate ${view}`}</span>
        <strong>{formatNumber(ratio, 2)}× <small>log₂ {log2Ratio >= 0 ? '+' : ''}{formatNumber(log2Ratio, 2)}</small></strong>
        <p>(summed SG TPM + 0.01) / (summed Total TPM + 0.01)</p>
      </div>
    </div>
    <p className="ratio-explainer">This ratio uses all {isoformRows.length} subset isoforms and is not a gene-level significance test. Linked filtering changes the cards below, not this summary.</p>
    <IsoformControls {...controls} />
    <div className="abundance-scale-note"><span>Shared vertical scale</span><b>0–{formatNumber(maxTx, 1)} TPM</b><small>bar = {view === 'mean' ? 'mean of 3' : `replicate ${view}`} · points = individual replicates</small></div>
    <div className="abundance-list abundance-list--vertical">
      {visibleRows.map((row) => <article className={`abundance-card abundance-card--vertical ${focusedTranscriptId === row.transcript_id ? 'is-focused' : ''}`}
        key={row.transcript_id} role="button" tabIndex="0" aria-pressed={focusedTranscriptId === row.transcript_id}
        onClick={() => onFocusTranscript(row.transcript_id)}
        onKeyDown={(event) => keyboardFocus(event, row.transcript_id, onFocusTranscript)}>
        <div className="abundance-card__name"><div className="transcript-title-line"><b>{row.transcript_name || row.transcript_id}</b>
          <span className={`call-chip call-${row.call}`}>{callLabel(row.call)}</span></div>
          <small>{row.transcript_id}</small></div>
        <div className="vertical-condition-chart">
          {['SG', 'Total'].map((condition) => {
            const value = getValue(row, condition);
            const replicates = byCondition[condition].map((sample) => ({
              replicate: sample.replicate,
              sampleId: sample.sample_id,
              value: row.expression?.samples?.[sample.sample_id]?.tpm ?? 0,
            }));
            return <div className={`vertical-condition ${condition === 'SG' ? 'is-sg' : 'is-total'}`} key={condition}>
              <div className="vertical-bar-stage" aria-label={`${condition} TPM; ${view === 'mean' ? 'mean' : `replicate ${view}`} ${value.toFixed(1)}`}>
                <span className="vertical-bar-value" style={{ bottom: `calc(${value / maxTx * 100}% + 8px)` }}>{formatNumber(value, 1)}</span>
                <div className="vertical-bar-shell"><i style={{ height: `${value / maxTx * 100}%` }} /></div>
                {replicates.map((replicate, index) => <span key={replicate.sampleId}
                  className={`replicate-dot ${view === String(replicate.replicate) ? 'is-selected' : ''}`}
                  style={{ left: `${22 + index * 28}%`, bottom: `calc(${replicate.value / maxTx * 100}% - 4px)` }}
                  title={`${replicate.sampleId}: ${replicate.value.toFixed(1)} TPM`}>
                  <i /><b>{formatNumber(replicate.value, 1)}</b>
                </span>)}
              </div>
              <strong className="vertical-condition__name"><i className={`condition-dot ${condition === 'SG' ? 'sg' : 'total'}`} />{condition === 'SG' ? 'Stress Granule' : 'Total RNA'}</strong>
              <div className="replicate-ledger" aria-label={`${condition} replicate values`}>
                {replicates.map((replicate) => <span key={replicate.sampleId} className={view === String(replicate.replicate) ? 'is-selected' : ''}>
                  R{replicate.replicate} <b>{formatNumber(replicate.value, 1)}</b>
                </span>)}
              </div>
            </div>;
          })}
        </div>
      </article>)}
      {!visibleRows.length && <div className="empty-filter-state"><b>No {gene.gene_name} isoforms match this Swish call.</b><span>Choose another linked filter to restore transcript rows.</span></div>}
    </div>
    <p className="panel-note">TPM comes from Salmon transcript quantification. Bars and points share the same scale across visible isoforms. Selecting a card focuses the same transcript in architecture and differential evidence; read-depth scaling is separate.</p>
  </section>;
}
