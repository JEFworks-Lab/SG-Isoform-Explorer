export default function LoadingSkeleton() {
  return (
    <div className="loading-grid">
      {[1, 2, 3].map((item) => <div key={item} className="panel skeleton-card"><i /><i /><i /></div>)}
    </div>
  );
}
