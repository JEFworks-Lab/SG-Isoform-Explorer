import { useState } from 'react';
import { CONDITION_ORDER, formatNumber } from '../lib/data.js';
import { PLOT_WIDTH, PLOT_LEFT, PLOT_RIGHT, locusX, locusTicks } from '../lib/scale.js';
import ReadPileup from './ReadPileup.jsx';

export function areaPath(row, start, end, baseline, amplitude, maxY) {
  if (!row?.positions?.length) return '';
  const points = row.positions.map((position, i) =>
    'L ' + locusX(position, start, end) + ' ' + (baseline - (row.values[i] / (maxY || 1)) * amplitude)).join(' ');
  return 'M ' + locusX(row.positions[0], start, end) + ' ' + baseline + ' ' + points +
    ' L ' + locusX(row.positions[row.positions.length - 1], start, end) + ' ' + baseline + ' Z';
}

export default function CoveragePlot({ gene, coverage, readRows, replicate, onChooseReplicate }) {
  const [mode, setMode] = useState('depth');
  if (!coverage.length) return null;
  const maxY = Math.max(1, ...coverage.flatMap((row) => row.values));
  const height = 235;
  const baseline = 187;
  const binSize = coverage[0]?.bin_size;

  return (
    <section className="panel wide-panel">
      <div className="section-heading">
        <div><span className="eyebrow">BAM evidence</span><h2>Unique-read coverage</h2></div>
        <div className="segmented-control" role="group" aria-label="Coverage display">
          <button className={mode === 'depth' ? 'is-active' : ''} aria-pressed={mode === 'depth'} onClick={() => setMode('depth')}>Normalized depth</button>
          <button className={mode === 'reads' ? 'is-active' : ''} aria-pressed={mode === 'reads'} onClick={() => { setMode('reads'); if (replicate === 'mean') onChooseReplicate('1'); }}>Individual reads</button>
        </div>
      </div>
      {mode === 'reads'
        ? <ReadPileup gene={gene} readRows={readRows} />
        : <div className="condition-tracks">
          {CONDITION_ORDER.map((condition) => {
            const row = coverage.find((item) => item.condition === condition);
            if (!row) return null;
            return <div className="plot-condition" key={condition}>
              <h3><span className={'condition-dot ' + (condition === 'SG' ? 'sg' : 'total')} />
                {condition === 'SG' ? 'Stress Granule RNA' : 'Total RNA'} · {row.sample_id || 'mean of ' + row.replicates + ' replicates'}</h3>
              <div className="chart-scroll"><svg className="coverage-svg" viewBox={'0 0 ' + PLOT_WIDTH + ' ' + height}
                role="img" aria-label={condition + ' coverage'}>
                {locusTicks(gene.start, gene.end).map((tick) => <g key={tick}>
                  <line x1={locusX(tick, gene.start, gene.end)} x2={locusX(tick, gene.start, gene.end)} y1="20" y2={baseline} className="grid-line" />
                  <text x={locusX(tick, gene.start, gene.end)} y={height - 12} textAnchor="middle" className="axis-label">{Math.round(tick).toLocaleString()}</text>
                </g>)}
                {[0, .5, 1].map((fraction) => <g key={fraction}>
                  <line x1={PLOT_LEFT} x2={PLOT_WIDTH - PLOT_RIGHT} y1={baseline - fraction * 145} y2={baseline - fraction * 145} className="grid-line" />
                  <text x={PLOT_LEFT - 12} y={baseline - fraction * 145 + 4} textAnchor="end" className="axis-label">{formatNumber(maxY * fraction)}</text>
                </g>)}
                <path d={areaPath(row, gene.start, gene.end, baseline, 145, maxY)}
                  className={'coverage-area ' + (condition === 'SG' ? 'coverage-area--sg' : 'coverage-area--total')} />
              </svg></div>
            </div>;
          })}
        </div>}
      {mode === 'depth' && <p className="panel-note">Mean uniquely mapped read depth per {binSize}-bp genomic bin, normalized per million STAR uniquely mapped reads; {replicate === 'mean' ? 'mean of three biological replicates' : 'replicate ' + replicate}. Both conditions share one y-axis and the same GRCh38 coordinate scale as every locus panel.</p>}
    </section>
  );
}
