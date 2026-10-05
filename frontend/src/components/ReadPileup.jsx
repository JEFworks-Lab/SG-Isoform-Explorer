import { CONDITION_ORDER } from '../lib/data.js';
import { PLOT_WIDTH, PLOT_LEFT, PLOT_RIGHT, locusX, locusTicks } from '../lib/scale.js';

export default function ReadPileup({ gene, readRows, junctionOnly = false }) {
  if (!readRows.length || !readRows.some((row) => row.reads)) {
    return <p className="panel-note">Individual alignments are unavailable for this dataset. Build read JSON from the matching BAM/SAM files to enable this view.</p>;
  }
  return <div className="read-condition-list">
    {CONDITION_ORDER.map((condition) => {
      const row = readRows.find((item) => item.sample.condition === condition);
      if (!row?.reads) return null;
      const allReads = junctionOnly ? (row.reads.spliced_reads || row.reads.reads.filter((read) => read.junctions.length)) : row.reads.reads;
      const shown = allReads.slice(0, 80);
      const height = 53 + shown.length * 13 + 32;
      return (
        <div className="read-condition" key={condition}>
          <h3><span className={'condition-dot ' + (condition === 'SG' ? 'sg' : 'total')} />
            {condition === 'SG' ? 'Stress Granule RNA' : 'Total RNA'} · {row.sample.sample_id}</h3>
          <p className="panel-note">{shown.length} displayed {junctionOnly ? 'spliced ' : ''}alignments; {junctionOnly ? row.reads.total_spliced_reads ?? shown.length : row.reads.total_reads} {junctionOnly ? 'spliced alignments' : 'alignments overlap this gene'} in this sample.
            {(junctionOnly ? row.reads.displayed_spliced_reads < row.reads.total_spliced_reads : row.reads.displayed_reads < row.reads.total_reads) && ' The JSON stores a reproducible random sample of at most 80 alignments.'}</p>
          {shown.length ? <div className="chart-scroll read-scroll">
            <svg viewBox={'0 0 ' + PLOT_WIDTH + ' ' + height} className="read-svg" role="img" aria-label={row.sample.sample_id + ' individual aligned reads'}>
              {locusTicks(gene.start, gene.end).map((tick) => (
                <g key={tick}>
                  <line x1={locusX(tick, gene.start, gene.end)} x2={locusX(tick, gene.start, gene.end)} y1="17" y2={height - 32} className="grid-line" />
                  <text x={locusX(tick, gene.start, gene.end)} y={height - 9} textAnchor="middle" className="axis-label">{Math.round(tick).toLocaleString()}</text>
                </g>
              ))}
              {shown.map((read, index) => {
                const y = 34 + index * 13;
                return <g className={'read-row ' + (condition === 'SG' ? 'read-row--sg' : 'read-row--total')} key={read.read_id + '-' + index}>
                  <title>{read.read_id + ' · ' + read.strand + ' strand · MAPQ ' + read.mapq + ' · ' + read.junctions.length + ' junction(s)'}</title>
                  {read.junctions.map((junction, j) =>
                    <line key={j} x1={locusX(junction.donor, gene.start, gene.end)} x2={locusX(junction.acceptor, gene.start, gene.end)} y1={y} y2={y} className="read-intron" />)}
                  {read.segments.map((segment, j) =>
                    <rect key={j} x={locusX(segment.start, gene.start, gene.end)} y={y - 3.5}
                      width={Math.max(2, locusX(segment.end, gene.start, gene.end) - locusX(segment.start, gene.start, gene.end) + 1)}
                      height="7" rx="1" className="read-segment" />)}
                </g>;
              })}
            </svg>
          </div> : <p className="panel-note">No spliced reads are in the displayed alignment subset.</p>}
        </div>
      );
    })}
    <p className="panel-note">Each bar is one aligned read segment. Thin connectors bridge spliced introns (N operations in CIGAR). Hover a read for its ID, strand, mapping quality and junction count. Reads shown are a capped sample, not a coverage estimate.</p>
  </div>;
}
