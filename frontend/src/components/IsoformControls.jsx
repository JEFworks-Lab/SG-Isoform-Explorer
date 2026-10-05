const CALL_FILTERS = [
  ['all', 'All'],
  ['significant', 'Significant'],
  ['enriched', 'SG enriched'],
  ['depleted', 'SG depleted'],
  ['not_significant', 'Not significant'],
  ['not_tested', 'Not tested'],
];

export default function IsoformControls({ filter, onFilter, counts, sortKey, onSortKey, direction, onDirection, sortOptions, visibleCount }) {
  return <div className="isoform-toolbar">
    <div className="linked-view-label"><span>Linked isoform view</span><small>Architecture · abundance · statistics</small></div>
    <div className="isoform-filter" role="group" aria-label="Filter isoforms by Swish call">
      {CALL_FILTERS.map(([value, label]) => <button type="button" key={value}
        className={filter === value ? 'is-active' : ''} aria-pressed={filter === value}
        onClick={() => onFilter(value)}>
        <span>{label}</span><b>{counts[value] ?? 0}</b>
      </button>)}
    </div>
    <div className="isoform-sort">
      <label>Sort isoforms
        <select value={sortKey} onChange={(event) => onSortKey(event.target.value)}>
          {sortOptions.map((option) => <option value={option.value} key={option.value}>{option.label}</option>)}
        </select>
      </label>
      <button type="button" className="sort-direction" onClick={onDirection}
        aria-label={`Sort ${direction === 'desc' ? 'ascending' : 'descending'}`}
        title={`Currently ${direction === 'desc' ? 'descending' : 'ascending'}`}>
        <span aria-hidden="true">{direction === 'desc' ? '↓' : '↑'}</span>
        {direction === 'desc' ? 'Descending' : 'Ascending'}
      </button>
      <span className="visible-count">Showing <b>{visibleCount}</b> of {counts.all}</span>
    </div>
  </div>;
}
