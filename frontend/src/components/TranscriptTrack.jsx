import { PLOT_WIDTH as WIDTH } from '../lib/scale.js';
import IsoformControls from './IsoformControls.jsx';
import { callLabel } from '../lib/isoforms.js';

export default function TranscriptTrack({ gene, visibleRows, locusScale, controls, focusedTranscriptId, onFocusTranscript }) {
  const regionStart = locusScale.start;
  const regionEnd = locusScale.end;
  const rowHeight = 55;
  const height = 38 + Math.max(1, visibleRows.length) * rowHeight + 34;
  const ticks = locusScale.ticks;
  const scale = locusScale.x;

  function chooseWithKeyboard(event, transcriptId) {
    if (event.key !== 'Enter' && event.key !== ' ') return;
    event.preventDefault();
    onFocusTranscript(transcriptId);
  }

  return <section className="panel track-panel">
    <div className="section-heading">
      <div><span className="eyebrow">Transcript architecture</span><h2>Isoform models</h2></div>
      <span className="coordinate-chip">{gene.chrom}:{regionStart.toLocaleString()}–{regionEnd.toLocaleString()} · {gene.strand}</span>
    </div>
    <p className="section-intro">Filter and sort here or in differential evidence; the linked view updates both sections and abundance together. Select a row to focus the same transcript across those panels.</p>
    <IsoformControls {...controls} />
    {visibleRows.length ? <div className="track-scroll">
      <svg className="transcript-svg" viewBox={`0 0 ${WIDTH} ${height}`} role="img" aria-label="Filtered and sorted transcript exon structures">
        {ticks.map((tick) => {
          const x = scale(tick, regionStart, regionEnd);
          return <g key={tick}><line x1={x} x2={x} y1={20} y2={height - 30} className="grid-line" />
            <text x={x} y={height - 9} textAnchor="middle" className="axis-label">{Math.round(tick).toLocaleString()}</text></g>;
        })}
        {locusScale.compressed.map((x, index) => <g key={'break-' + index} className="axis-break" transform={`translate(${x},${height - 28})`}><path d="M -6 4 l 5 -8 M 1 4 l 5 -8" /></g>)}
        {visibleRows.map((tx, index) => {
          const y = 45 + index * rowHeight;
          const x1 = scale(Math.max(tx.start, regionStart));
          const x2 = scale(Math.min(tx.end + 1, regionEnd + 1));
          const focused = tx.transcript_id === focusedTranscriptId;
          return <g key={tx.transcript_id} className={`transcript-model ${focused ? 'is-focused' : ''}`}
            tabIndex="0" role="button" aria-pressed={focused}
            aria-label={`${tx.transcript_name || tx.transcript_id}, ${callLabel(tx.call)}, ${tx.exons.length} exons`}
            onClick={() => onFocusTranscript(tx.transcript_id)}
            onKeyDown={(event) => chooseWithKeyboard(event, tx.transcript_id)}>
            {focused && <rect x="3" y={y - 21} width={WIDTH - 6} height="42" rx="7" className="transcript-focus-bg" />}
            <text x="12" y={y + 5} className="tx-label">{tx.transcript_name || tx.transcript_id}<title>{tx.transcript_id}</title></text>
            {x1 < x2 && <line x1={x1} x2={x2} y1={y} y2={y} className="intron-line" />}
            {Array.from({ length: Math.max(0, Math.floor((x2 - x1) / 42)) }, (_, arrowIndex) => {
              const frac = (arrowIndex + 1) / (Math.max(1, Math.floor((x2 - x1) / 42)) + 1);
              const x = gene.strand === '-' ? x2 - frac * (x2 - x1) : x1 + frac * (x2 - x1);
              return <text key={arrowIndex} x={x} y={y - 5} textAnchor="middle" className="strand-arrow">{gene.strand === '-' ? '‹' : '›'}</text>;
            })}
            {tx.exons.filter((exon) => exon.start <= regionEnd && exon.end >= regionStart).map((exon, exonIndex) => {
              const exonX = scale(Math.max(exon.start, regionStart));
              const exonW = Math.max(2, scale(Math.min(exon.end + 1, regionEnd + 1)) - exonX);
              return <g key={`${tx.transcript_id}-${exon.start}`}><rect x={exonX} y={y - 11} width={exonW} height="22" rx="3" className="exon-block" />
                <title>{`${tx.transcript_id} exon ${exon.exon_number ?? exonIndex + 1}: ${exon.start}-${exon.end}`}</title></g>;
            })}
          </g>;
        })}
      </svg>
    </div> : <div className="empty-filter-state"><b>No transcript models match this Swish call.</b><span>Choose another linked filter to restore models.</span></div>}
    <p className="panel-note">■ Exon · connecting line Intron · › Transcription direction. Models come from the supplied Ensembl 115 subset GTFs; filtering and focusing change only the view, not the annotation or measurements.</p>
  </section>;
}
