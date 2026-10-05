import { useEffect, useMemo, useState } from 'react';
import { formatNumber } from '../lib/data.js';
import { displayedCoverage } from '../lib/coverage.js';
import { PLOT_WIDTH, PLOT_LEFT, PLOT_RIGHT } from '../lib/scale.js';

const HEIGHT = 370;
const BASELINE = 137;
const AMP = 103;
const PLOT_BOTTOM = 330;
const AXIS_Y = 359;
const BREAK_Y = 342;
const LABEL_HEIGHT = 22;

function junctionKey(junction) {
  return `${junction.donor}:${junction.acceptor}:${junction.strand_code}`;
}

function areaPath(points, scale, maxDepth) {
  if (!points.length) return '';
  return `M ${scale.x(points[0].position)} ${BASELINE} ` + points.map((point) =>
    `L ${scale.x(point.position)} ${BASELINE - point.value / maxDepth * AMP} L ${scale.x(Math.min(scale.end + 1, point.end + 1))} ${BASELINE - point.value / maxDepth * AMP}`).join(' ') +
    ` L ${scale.x(Math.min(scale.end + 1, points[points.length - 1].end + 1))} ${BASELINE} Z`;
}

function junctionGeometry(junction, scale) {
  const x1 = scale.x(junction.donor);
  const x2 = scale.x(junction.acceptor);
  const depth = 24 + Math.min(72, (x2 - x1) * .11);
  const labelText = junction.count === 0 ? `0u · ${junction.multi_mapping_count}m` : String(junction.count);
  return {
    junction,
    x1,
    x2,
    depth,
    anchorX: (x1 + x2) / 2,
    anchorY: BASELINE + depth,
    path: `M ${x1} ${BASELINE} Q ${(x1 + x2) / 2} ${BASELINE + depth * 2} ${x2} ${BASELINE}`,
    labelText,
    labelWidth: Math.max(44, labelText.length * 8 + 18),
  };
}

// Keep count badges on their own arc apex. When the plot is too dense to show
// every badge without overlap, the remaining counts stay available by
// hover/focus/click and in the complete junction table below the plot.
function layoutJunctionLabels(junctions, scale) {
  const geometry = junctions.map((junction) => junctionGeometry(junction, scale));
  const occupied = [];
  const visibleKeys = new Set();
  const ranked = [...geometry].sort((a, b) =>
    b.junction.count - a.junction.count || b.depth - a.depth || a.anchorX - b.anchorX);
  for (const item of ranked) {
    const collides = occupied.some((other) =>
      Math.abs(item.anchorX - other.anchorX) < (item.labelWidth + other.labelWidth) / 2 + 6 &&
      Math.abs(item.anchorY - other.anchorY) < LABEL_HEIGHT + 5,
    );
    if (!collides) {
      occupied.push(item);
      visibleKeys.add(junctionKey(item.junction));
    }
  }
  return geometry.map((item) => ({ ...item, showLabel: visibleKeys.has(junctionKey(item.junction)) }));
}

function IsoformOverlay({ transcript, transcripts, onChange, scale }) {
  return <div className="sample-isoform-overlay">
    <div className="sample-isoform-overlay__head">
      <div><span>One-isoform annotation overlay</span><small>Aligned to the coverage and arc coordinates above</small></div>
      <label>Isoform
        <select value={transcript?.transcript_id ?? ''} onChange={(event) => onChange(event.target.value)} disabled={!transcripts.length}>
          {!transcripts.length && <option value="">No isoforms in linked filter</option>}
          {transcripts.map((tx) => <option value={tx.transcript_id} key={tx.transcript_id}>{tx.transcript_name || tx.transcript_id}</option>)}
        </select>
      </label>
    </div>
    {transcript && <div className="sample-isoform-overlay__scroll"><svg viewBox={`0 0 ${PLOT_WIDTH} 70`} role="img"
      aria-label={`${transcript.transcript_name || transcript.transcript_id} exon annotation overlay`}>
      <line x1={PLOT_LEFT} x2={PLOT_WIDTH - PLOT_RIGHT} y1="35" y2="35" className="isoform-overlay-axis" />
      {transcript.exons.filter((exon) => exon.start <= scale.end && exon.end >= scale.start).map((exon, index) => {
        const x = scale.x(Math.max(exon.start, scale.start));
        const width = Math.max(2, scale.x(Math.min(exon.end + 1, scale.end + 1)) - x);
        return <g key={`${transcript.transcript_id}-${exon.start}`}><rect x={x} y="23" width={width} height="24" rx="3" className="isoform-overlay-exon" />
          <title>{`Exon ${exon.exon_number ?? index + 1}: ${exon.start.toLocaleString()}–${exon.end.toLocaleString()}`}</title></g>;
      })}
      {scale.compressed.map((x, index) => <g key={index} className="axis-break" transform={`translate(${x},35)`}><path d="M -6 4 l 5 -8 M 1 4 l 5 -8" /></g>)}
      <text x={PLOT_LEFT} y="64" className="isoform-overlay-label">{transcript.transcript_id} · {transcript.exons.length} exons</text>
    </svg></div>}
  </div>;
}

function SampleTrack({ row, scale, points, maxDepth, maxSupport, normalized, junctionMode, minimumUnique,
  transcripts, focusedTranscriptId }) {
  const [coverageHover, setCoverageHover] = useState(null);
  const [hoveredJunction, setHoveredJunction] = useState(null);
  const [focusedJunction, setFocusedJunction] = useState(null);
  const [pinnedJunction, setPinnedJunction] = useState(null);
  const [selectedTranscriptId, setSelectedTranscriptId] = useState(focusedTranscriptId ?? transcripts[0]?.transcript_id ?? '');
  const { sample, coverage, junctions } = row;
  const isSG = sample.condition === 'SG';
  const modeJunctions = useMemo(() => junctions.junctions
    .filter((junction) => junction.donor >= scale.start && junction.acceptor <= scale.end &&
      (junctionMode === 'all' || junction.matches_subset_annotation))
    .sort((a, b) => (b.acceptor - b.donor) - (a.acceptor - a.donor)),
  [junctions, scale, junctionMode]);
  const visibleJunctions = useMemo(() => modeJunctions.filter((junction) => junction.count >= minimumUnique),
    [modeJunctions, minimumUnique]);
  const labels = useMemo(() => layoutJunctionLabels(visibleJunctions, scale), [visibleJunctions, scale]);
  const activeKey = hoveredJunction || focusedJunction || pinnedJunction;
  const activeJunction = visibleJunctions.find((junction) => junctionKey(junction) === activeKey);
  const hiddenMultiOnly = modeJunctions.filter((junction) => junction.count === 0).length;
  const selectedTranscript = transcripts.find((transcript) => transcript.transcript_id === selectedTranscriptId);

  useEffect(() => {
    if (!transcripts.some((transcript) => transcript.transcript_id === selectedTranscriptId)) {
      const replacement = transcripts.find((transcript) => transcript.transcript_id === focusedTranscriptId) ?? transcripts[0];
      setSelectedTranscriptId(replacement?.transcript_id ?? '');
    }
  }, [transcripts, focusedTranscriptId, selectedTranscriptId]);

  function supportFor(junction) {
    return normalized ? junction.count * 1e6 / coverage.normalization_denominator : junction.count;
  }
  function togglePinnedJunction(key) {
    setPinnedJunction((current) => current === key ? null : key);
  }
  function activateWithKeyboard(event, key) {
    if (event.key !== 'Enter' && event.key !== ' ') return;
    event.preventDefault();
    togglePinnedJunction(key);
  }
  function accessibleLabel(junction) {
    const support = junction.count * 1e6 / coverage.normalization_denominator;
    const evidence = junction.count === 0
      ? `zero uniquely mapping reads and ${junction.multi_mapping_count} multimapping reads; multimapping-only evidence`
      : `${junction.count} uniquely mapping reads and ${junction.multi_mapping_count} multimapping reads`;
    return `${junction.donor.toLocaleString()} to ${junction.acceptor.toLocaleString()}: ${evidence}, ${formatNumber(support, 2)} unique reads per million STAR uniquely mapped reads, ${junction.matches_subset_annotation ? 'matches a displayed exon edge' : 'not represented in the displayed subset models'}`;
  }
  function move(event) {
    const rect = event.currentTarget.getBoundingClientRect();
    const svgX = (event.clientX - rect.left) / rect.width * PLOT_WIDTH;
    if (svgX < PLOT_LEFT || svgX > PLOT_WIDTH - PLOT_RIGHT) { setCoverageHover(null); return; }
    const coordinate = scale.invert(svgX);
    const index = coordinate - coverage.start;
    setCoverageHover({ coordinate, x: scale.x(coordinate), raw: coverage.values[index] ?? 0,
      adjusted: coverage.normalized_values[index] ?? 0 });
  }

  return <div className="sashimi-condition">
    <div className="sashimi-condition__head">
      <span className={'condition-dot ' + (isSG ? 'sg' : 'total')} />
      <div><b>{isSG ? 'Stress Granule RNA' : 'Total RNA'} · replicate {sample.replicate}</b>
        <small>{sample.sample_id} · {coverage.normalization_denominator.toLocaleString()} STAR uniquely mapped reads</small></div>
    </div>
    <div className={`junction-inspector ${activeJunction ? 'is-active' : ''}`} aria-live="polite">
      {activeJunction ? <>
        <span className="junction-inspector__eyebrow">Selected junction</span>
        <b>{activeJunction.donor.toLocaleString()} → {activeJunction.acceptor.toLocaleString()}</b>
        <strong>{activeJunction.count} unique</strong>
        <span>{activeJunction.multi_mapping_count} multi-map</span>
        <span>{activeJunction.matches_subset_annotation ? 'Subset edge' : 'Outside displayed subset'}</span>
        {pinnedJunction === activeKey && <button type="button" onClick={() => setPinnedJunction(null)}>Clear pin</button>}
      </> : <>
        <span>Hover, focus, or click an arc or table row to trace its endpoints and counts.</span>
        <small>{visibleJunctions.length} arcs shown{minimumUnique > 0 && hiddenMultiOnly ? ` · ${hiddenMultiOnly} multimapping-only row${hiddenMultiOnly === 1 ? '' : 's'} hidden` : ''}</small>
      </>}
    </div>
    <div className="chart-scroll"><svg viewBox={`0 0 ${PLOT_WIDTH} ${HEIGHT}`} className="sashimi-svg" role="group"
      aria-label={`${sample.sample_id} per-base coverage and splice junctions`} onMouseMove={move} onMouseLeave={() => setCoverageHover(null)}>
      {scale.ticks.map((tick) => <g key={tick}>
        <line x1={scale.x(tick)} x2={scale.x(tick)} y1="15" y2={PLOT_BOTTOM} className="grid-line" />
        <text x={scale.x(tick)} y={AXIS_Y} textAnchor="middle" className="axis-label">{tick.toLocaleString()}</text>
      </g>)}
      {[0, .5, 1].map((fraction) => <g key={fraction}>
        <text x={PLOT_LEFT - 12} y={BASELINE - fraction * AMP + 4} textAnchor="end" className="axis-label">{formatNumber(maxDepth * fraction, normalized ? 1 : 0)}</text>
      </g>)}
      <path d={areaPath(points, scale, maxDepth)} className={'coverage-area ' + (isSG ? 'coverage-area--sg' : 'coverage-area--total')} />
      <line x1={PLOT_LEFT} x2={PLOT_WIDTH - PLOT_RIGHT} y1={BASELINE} y2={BASELINE} className="locus-line" />

      {labels.map(({ junction, path, x1, x2, anchorX, anchorY, labelWidth, labelText, showLabel }) => {
        const key = junctionKey(junction);
        const isActive = activeKey === key;
        const support = supportFor(junction);
        return <g key={key}
          className={`junction-item ${isActive ? 'is-active' : ''} ${activeKey && !isActive ? 'is-muted' : ''} ${junction.count === 0 ? 'is-multi-only' : ''}`}
          tabIndex="0" role="button" aria-pressed={pinnedJunction === key} aria-label={accessibleLabel(junction)}
          onMouseEnter={() => setHoveredJunction(key)} onMouseLeave={() => setHoveredJunction(null)}
          onFocus={() => setFocusedJunction(key)} onBlur={() => setFocusedJunction(null)}
          onKeyDown={(event) => activateWithKeyboard(event, key)}
          onClick={(event) => { event.stopPropagation(); togglePinnedJunction(key); }}>
          <title>{accessibleLabel(junction)}</title>
          <path d={path} fill="none" className="junction-hit" />
          <path d={path} fill="none" className={'junction-arc ' + (isSG ? 'junction-arc--sg' : 'junction-arc--total')}
            style={{ strokeWidth: 1.4 + support / maxSupport * 5 }} />
          {isActive && <>
            <line x1={x1} x2={x1} y1={BASELINE - 4} y2={PLOT_BOTTOM} className="junction-endpoint-guide" />
            <line x1={x2} x2={x2} y1={BASELINE - 4} y2={PLOT_BOTTOM} className="junction-endpoint-guide" />
            <circle cx={x1} cy={BASELINE} r="4" className="junction-endpoint" />
            <circle cx={x2} cy={BASELINE} r="4" className="junction-endpoint" />
          </>}
          {(showLabel || isActive) && <g className="junction-label-group">
            <rect x={anchorX - labelWidth / 2} y={anchorY - LABEL_HEIGHT / 2} width={labelWidth} height={LABEL_HEIGHT} rx="10" className="junction-label-bg" />
            <text x={anchorX} y={anchorY + 4} textAnchor="middle" className="junction-label">{labelText}</text>
          </g>}
        </g>;
      })}
      {scale.compressed.map((x, index) => <g key={index} className="axis-break" transform={`translate(${x},${BREAK_Y})`}>
        <path d="M -6 4 l 5 -8 M 1 4 l 5 -8" /></g>)}
      {coverageHover && !activeKey && <g pointerEvents="none">
        <line x1={coverageHover.x} x2={coverageHover.x} y1="17" y2={BASELINE} className="hover-line" />
        <rect x={Math.min(coverageHover.x + 7, PLOT_WIDTH - 294)} y="16" width="283" height="45" rx="6" className="hover-box" />
        <text x={Math.min(coverageHover.x + 16, PLOT_WIDTH - 285)} y="34" className="hover-label">{`${sample.sample_id} · ${coverageHover.coordinate.toLocaleString()} bp`}</text>
        <text x={Math.min(coverageHover.x + 16, PLOT_WIDTH - 285)} y="51" className="hover-label">{`Raw ${coverageHover.raw} · per million ${formatNumber(coverageHover.adjusted, 2)}`}</text>
      </g>}
    </svg></div>

    <IsoformOverlay transcript={selectedTranscript} transcripts={transcripts} onChange={setSelectedTranscriptId} scale={scale} />

    <details className="junction-access">
      <summary>Junction table <span>{visibleJunctions.length} displayed arcs · coordinates and full counts</span></summary>
      <div className="junction-access__scroller">
        <div className="junction-access__table" role="table" aria-label={`${sample.sample_id} displayed splice junctions`}>
          <div className="junction-access__row junction-access__head" role="row">
            <span role="columnheader">Donor</span><span role="columnheader">Acceptor</span><span role="columnheader">Unique</span><span role="columnheader">Multi-map</span><span role="columnheader">Unique / million</span><span role="columnheader">Evidence</span><span role="columnheader">Subset edge</span>
          </div>
          {[...visibleJunctions].sort((a, b) => b.count - a.count).map((junction) => {
            const key = junctionKey(junction);
            return <div className={`junction-access__row ${activeKey === key ? 'is-active' : ''}`} role="row" key={key} tabIndex="0"
              aria-label={`${accessibleLabel(junction)}. Press Enter to pin this junction.`}
              onMouseEnter={() => setHoveredJunction(key)} onMouseLeave={() => setHoveredJunction(null)}
              onFocus={() => setFocusedJunction(key)} onBlur={() => setFocusedJunction(null)}
              onKeyDown={(event) => activateWithKeyboard(event, key)}
              onClick={() => togglePinnedJunction(key)}>
            <span role="cell">{junction.donor.toLocaleString()}</span><span role="cell">{junction.acceptor.toLocaleString()}</span>
            <strong role="cell">{junction.count.toLocaleString()}</strong><span role="cell">{junction.multi_mapping_count.toLocaleString()}</span><span role="cell">{formatNumber(junction.count * 1e6 / coverage.normalization_denominator, 3)}</span>
            <span role="cell" className={junction.count === 0 ? 'evidence-multi-only' : ''}>{junction.count === 0 ? 'Multi-map only' : 'Unique support'}</span>
            <span role="cell">{junction.matches_subset_annotation ? 'Yes' : 'No'}</span>
          </div>;})}
        </div>
      </div>
    </details>
  </div>;
}

export default function SashimiPlot({ trackRows, samples, scale, transcripts, focusedTranscriptId }) {
  const [selected, setSelected] = useState(() => ['SG', 'Total'].map((condition) =>
    samples.filter((sample) => sample.condition === condition).sort((a, b) => Number(a.replicate) - Number(b.replicate))[0]?.sample_id,
  ).filter(Boolean));
  const [normalized, setNormalized] = useState(true);
  const [binSize, setBinSize] = useState(1);
  const [junctionMode, setJunctionMode] = useState('subset');
  const [minimumUnique, setMinimumUnique] = useState(1);
  const chosen = trackRows.filter((row) => selected.includes(row.sample.sample_id))
    .sort((a, b) => (a.sample.condition === b.sample.condition ? Number(a.sample.replicate) - Number(b.sample.replicate) : a.sample.condition === 'SG' ? -1 : 1));
  const plotted = useMemo(() => chosen.map((row) => ({ row, points: displayedCoverage(row.coverage, scale, binSize, normalized) })),
    [trackRows, selected, scale, binSize, normalized]);
  const maxDepth = Math.max(1, ...plotted.flatMap(({ points }) => points.map((point) => point.value)));
  const maxSupport = Math.max(1, ...chosen.flatMap(({ coverage, junctions }) => junctions.junctions
    .filter((junction) => junction.donor >= scale.start && junction.acceptor <= scale.end && (junctionMode === 'all' || junction.matches_subset_annotation))
    .filter((junction) => junction.count >= minimumUnique)
    .map((junction) => normalized ? junction.count * 1e6 / coverage.normalization_denominator : junction.count)));
  function toggleSample(id) {
    setSelected((previous) => previous.includes(id) ? (previous.length > 1 ? previous.filter((item) => item !== id) : previous) : [...previous, id]);
  }
  return <section className="panel wide-panel">
    <div className="section-heading">
      <div><span className="eyebrow">Read depth + observed splice junctions</span><h2>Sashimi view</h2></div>
      <span className="unit-chip">{normalized ? 'Depth per million STAR unique reads' : 'Raw read depth'} · {binSize} bp</span>
    </div>
    <div className="sashimi-controls">
      <fieldset className="sample-fieldset"><legend>Show individual sample tracks</legend>
        {['SG', 'Total'].map((condition) => <div className="sample-group" key={condition}>
          <strong>{condition === 'SG' ? 'Stress Granule' : 'Total RNA'}</strong>
          {samples.filter((sample) => sample.condition === condition).sort((a, b) => Number(a.replicate) - Number(b.replicate)).map((sample) =>
            <label key={sample.sample_id}><input type="checkbox" checked={selected.includes(sample.sample_id)} onChange={() => toggleSample(sample.sample_id)} /> Replicate {sample.replicate}</label>)}
        </div>)}
      </fieldset>
      <div className="plot-settings"><label>Depth scale <select value={normalized ? 'normalized' : 'raw'} onChange={(event) => setNormalized(event.target.value === 'normalized')}>
        <option value="normalized">Per million STAR unique reads</option><option value="raw">Raw read depth</option>
      </select></label>
      <label>Display bin <select value={binSize} onChange={(event) => setBinSize(Number(event.target.value))}>
        <option value="1">1 bp · exact</option><option value="5">5 bp · mean depth</option><option value="25">25 bp · mean depth</option>
      </select></label>
      <label>Junctions <select value={junctionMode} onChange={(event) => setJunctionMode(event.target.value)}>
        <option value="subset">Subset-matched arcs</option><option value="all">All final STAR arcs</option>
      </select></label>
      <label>Unique-read filter <select value={minimumUnique} onChange={(event) => setMinimumUnique(Number(event.target.value))}>
        <option value="1">≥ 1 unique · default</option><option value="2">≥ 2 unique</option><option value="5">≥ 5 unique</option><option value="0">Include 0 unique</option>
      </select></label></div>
    </div>
    <div className="sashimi-grid">{plotted.map(({ row, points }) => <SampleTrack key={row.sample.sample_id} row={row} points={points}
      scale={scale} maxDepth={maxDepth} maxSupport={maxSupport} normalized={normalized} junctionMode={junctionMode} minimumUnique={minimumUnique}
      transcripts={transcripts} focusedTranscriptId={focusedTranscriptId} />)}</div>
    <p className="panel-note">Shared depth scale across selected samples. Each one-isoform strip is an annotation overlay aligned to the same coordinate scale; choosing it does not assign coverage or junction reads to that transcript. Non-overlapping count badges stay on their arc; hover, focus, or click an arc or table row to trace its endpoints. The default hides STAR rows with zero uniquely mapping reads. Arc labels are raw unique counts from final STAR <code>SJ.out.tab</code>, placed on adjacent exonic bases.</p>
  </section>;
}
