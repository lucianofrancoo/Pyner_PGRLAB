import type { MineroResponse, ProAnalysisResponse, RelevanceLabel } from '../types';
import { PaperNetwork } from './PaperNetwork';

interface AnalysisViewProps {
  response: MineroResponse | null;
  proResponse: ProAnalysisResponse | null;
}

const EVIDENCE_DISPLAY: Record<string, string> = {
  directa: 'Direct',
  indirecta: 'Indirect',
  'débil': 'Weak',
};

function percent(value: number, total: number): number {
  if (!total) return 0;
  return Math.round((value / total) * 100);
}

export function AnalysisView({ response, proResponse }: AnalysisViewProps) {
  if (!response) {
    return (
      <section className="panel">
        <header className="panel-header">
          <h2>Analytics</h2>
          <p>Run a search to enable relevance and evidence metrics.</p>
        </header>
      </section>
    );
  }

  const total = response.results.length;
  const byLabel: Record<RelevanceLabel, number> = { alta: 0, media: 0, baja: 0 };
  const byEvidence: Record<string, number> = { directa: 0, indirecta: 0, 'débil': 0 };

  response.results.forEach((result) => {
    byLabel[result.classification.relevance_label] += 1;
    byEvidence[result.classification.evidence_level] =
      (byEvidence[result.classification.evidence_level] ?? 0) + 1;
  });

  const authorByPmid: Record<string, string> = {};
  response.results.forEach((result) => {
    if ('pmid' in result && Array.isArray(result.authors) && result.authors.length > 0) {
      authorByPmid[result.pmid] = result.authors[0];
    }
  });

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Analytics</h2>
        <p>Compact quality overview for fast scientific decision-making.</p>
      </header>

      <div className="kpi-grid">
        <article>
          <h3>Total results</h3>
          <p>{total}</p>
        </article>
        <article>
          <h3>Classification mode</h3>
          <p>{response.metadata.model_default}</p>
        </article>
        <article>
          <h3>LLM available</h3>
          <p>{response.metadata.llm_runtime_available ? 'Yes' : 'No'}</p>
        </article>
      </div>

      <div className="chart-grid">
        <article className="bar-panel">
          <h3>Relevance</h3>
          {(['alta', 'media', 'baja'] as RelevanceLabel[]).map((label) => {
            const count = byLabel[label];
            return (
              <div key={label} className="bar-row">
                <span>{label}</span>
                <div className="bar-track">
                  <div className={`bar-fill ${label}`} style={{ width: `${percent(count, total)}%` }} />
                </div>
                <strong>{count}</strong>
              </div>
            );
          })}
        </article>

        <article className="bar-panel">
          <h3>Evidence level</h3>
          {Object.entries(byEvidence).map(([level, count]) => (
            <div key={level} className="bar-row">
              <span>{EVIDENCE_DISPLAY[level] ?? level}</span>
              <div className="bar-track">
                <div className="bar-fill evidencia" style={{ width: `${percent(count, total)}%` }} />
              </div>
              <strong>{count}</strong>
            </div>
          ))}
        </article>
      </div>

      {proResponse?.results?.length ? (
        <article className="analytics-network">
          <header className="analytics-network-header">
            <h3>Paper Network (Pro)</h3>
            <p>
              Integrated graph with native filters for relevance, year and journal.
            </p>
          </header>
          <div className="analytics-network-frame">
            <PaperNetwork
              papers={proResponse.results}
              authorByPmid={authorByPmid}
              title="PYNER Paper Network"
              subtitle="Integrated in Minero web app (no embedded page)."
              compact
            />
          </div>
        </article>
      ) : null}
    </section>
  );
}
