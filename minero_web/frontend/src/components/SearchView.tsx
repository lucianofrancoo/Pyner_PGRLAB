import { useEffect, useMemo, useRef, useState } from 'react';
import { FlaskConical, Search, Pickaxe } from 'lucide-react';
import type { QueryGeneration, SearchPayload, SourceMode } from '../types';

interface SearchViewProps {
  loading: boolean;
  isGenerating: boolean;
  isRunning: boolean;
  llmAvailable: boolean;
  pendingQuery: QueryGeneration | null;
  onGenerate: (payload: SearchPayload) => Promise<void>;
  onConfirm: () => Promise<void>;
  onDiscard: () => void;
}

const DEFAULT_QUERY =
  'Arabidopsis thaliana RNA-Seq under drought and water stress in root';

const MAX_HARD_LIMIT = 10000;

const SOURCE_DESCRIPTIONS: Record<string, string> = {
  pmc: 'PMC: full-text search — more results, slower.',
  bioproject: 'BioProject: slower, focused on projects and SRA hierarchy.',
};

const LINE_WINDOW_POINTS = 28;

type QuerySegment = {
  text: string;
  type: 'plain' | 'group';
  groupIndex?: number;
};

function splitByTopLevelParentheses(query: string): QuerySegment[] {
  if (!query.trim()) return [{ text: query, type: 'plain' }];

  const segments: QuerySegment[] = [];
  let cursor = 0;
  let groupIndex = 0;

  while (cursor < query.length) {
    if (query[cursor] !== '(') {
      const nextParen = query.indexOf('(', cursor);
      const end = nextParen === -1 ? query.length : nextParen;
      segments.push({ text: query.slice(cursor, end), type: 'plain' });
      cursor = end;
      continue;
    }

    let depth = 0;
    let end = cursor;
    while (end < query.length) {
      const char = query[end];
      if (char === '(') depth += 1;
      if (char === ')') depth -= 1;
      end += 1;
      if (depth === 0) break;
    }

    segments.push({
      text: query.slice(cursor, end),
      type: 'group',
      groupIndex,
    });
    groupIndex += 1;
    cursor = end;
  }

  return segments;
}

export function SearchView({
  loading,
  isGenerating,
  isRunning,
  llmAvailable,
  pendingQuery,
  onGenerate,
  onConfirm,
  onDiscard,
}: SearchViewProps) {
  const [naturalQuery, setNaturalQuery] = useState(DEFAULT_QUERY);
  const [source, setSource] = useState<SourceMode>('pmc');
  const [useLlm, setUseLlm] = useState(true);
  // '' = unlimted input display; internally sent as MAX_HARD_LIMIT
  const [maxResultsInput, setMaxResultsInput] = useState<string>('');

  const resolvedMaxResults = useMemo(() => {
    const parsed = parseInt(maxResultsInput, 10);
    if (!maxResultsInput.trim() || isNaN(parsed) || parsed <= 0) return MAX_HARD_LIMIT;
    return Math.min(parsed, MAX_HARD_LIMIT);
  }, [maxResultsInput]);

  const sourceDescription = SOURCE_DESCRIPTIONS[source];
  const isSearchRunning = isRunning && Boolean(pendingQuery);
  const querySegments = useMemo(
    () => splitByTopLevelParentheses(pendingQuery?.ncbi_query ?? ''),
    [pendingQuery?.ncbi_query]
  );
  const compactQuery = useMemo(
    () => (pendingQuery?.ncbi_query ?? '').replace(/\s+/g, ' ').trim(),
    [pendingQuery?.ncbi_query]
  );

  const [processedCount, setProcessedCount] = useState(0);
  const [lineSeries, setLineSeries] = useState<number[]>([0]);
  const processedRef = useRef(0);
  const estimatedTarget = useMemo(() => {
    if (resolvedMaxResults >= MAX_HARD_LIMIT) {
      return source === 'bioproject' ? 3500 : 2200;
    }
    return Math.max(120, Math.min(resolvedMaxResults, 8000));
  }, [resolvedMaxResults, source]);

  const linePoints = useMemo(() => {
    const safeSeries = lineSeries.length ? lineSeries : [0];
    const maxValue = Math.max(estimatedTarget, ...safeSeries, 1);
    return safeSeries
      .map((value, index) => {
        const x = safeSeries.length === 1 ? 0 : (index / (safeSeries.length - 1)) * 100;
        const y = 38 - (value / maxValue) * 32;
        return `${x.toFixed(2)},${y.toFixed(2)}`;
      })
      .join(' ');
  }, [lineSeries, estimatedTarget]);

  const areaPoints = useMemo(() => {
    if (!linePoints) return '';
    return `0,38 ${linePoints} 100,38`;
  }, [linePoints]);

  useEffect(() => {
    if (!isSearchRunning) {
      processedRef.current = 0;
      setProcessedCount(0);
      setLineSeries([0]);
      return;
    }

    const timer = window.setInterval(() => {
      const current = processedRef.current;
      const baseJump = Math.max(2, Math.floor(estimatedTarget * 0.018));
      const next =
        current >= estimatedTarget
          ? estimatedTarget
          : Math.min(estimatedTarget, current + baseJump + Math.floor(Math.random() * 12));

      processedRef.current = next;
      setProcessedCount(next);
      setLineSeries((prev) => {
        const appended = [...prev, next];
        if (appended.length > LINE_WINDOW_POINTS) appended.shift();
        return appended;
      });
    }, 650);

    return () => window.clearInterval(timer);
  }, [isSearchRunning, estimatedTarget]);

  async function handleGenerate(): Promise<void> {
    await onGenerate({
      natural_query: naturalQuery.trim(),
      source,
      max_results: resolvedMaxResults,
      use_llm: useLlm,
    });
  }

  async function handleSubmit(event: React.FormEvent): Promise<void> {
    event.preventDefault();
    await handleGenerate();
  }

  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Search</h2>
        <p>Write your biological question in natural language and run Minero.</p>
      </header>

      <form className={`search-form${isSearchRunning ? ' search-form--minimized' : ''}`} onSubmit={handleSubmit}>
        <label className="field">
          <span>Biological query</span>
          <textarea
            value={naturalQuery}
            onChange={(event) => setNaturalQuery(event.target.value)}
            placeholder="Ex: tomato RNA-Seq under phosphate deficiency"
            rows={6}
            required
          />
        </label>

        <div className="field-row">
          <label className="field">
            <span>Data source</span>
            <select value={source} onChange={(event) => setSource(event.target.value as SourceMode)}>
              <option value="pmc">PMC (full-text)</option>
              <option value="bioproject">BioProject</option>
            </select>
            <small>{sourceDescription}</small>
          </label>

          <label className="field">
            <span>Max results</span>
            <input
              type="number"
              min={1}
              max={MAX_HARD_LIMIT}
              value={maxResultsInput}
              onChange={(e) => setMaxResultsInput(e.target.value)}
              placeholder={`Unlimited (up to ${MAX_HARD_LIMIT.toLocaleString()})`}
            />
            <small>
              {maxResultsInput.trim() && resolvedMaxResults < MAX_HARD_LIMIT
                ? `Will fetch up to ${resolvedMaxResults} records.`
                : 'Leave blank for unlimited (max 10,000).'}
            </small>
          </label>
        </div>

        <label className="switch-field">
          <input
            type="checkbox"
            checked={useLlm}
            onChange={(event) => setUseLlm(event.target.checked)}
          />
          <span>
            LLM Classification (Ollama)
            <em>
              {llmAvailable
                ? 'Available. If it fails, Minero uses heuristic fallback.'
                : 'Unavailable right now. Heuristic fallback will be used.'}
            </em>
          </span>
        </label>

        <div className="actions">
          <button type="submit" className="primary" disabled={loading || naturalQuery.trim().length < 3}>
            {isGenerating ? <Pickaxe size={16} className="chop" /> : <Search size={16} />}
            {isGenerating ? 'Generating query...' : 'Generate query'}
          </button>
        </div>
      </form>

      {pendingQuery ? (
        <section className="panel compact" style={{ marginTop: 12 }}>
          <header className="panel-header">
            <h2>Confirmation</h2>
            <p>
              {isSearchRunning
                ? 'Search is running. Query preview collapsed while records are being processed.'
                : 'Review the generated NCBI query and confirm to run the real search.'}
            </p>
          </header>
          {isSearchRunning ? (
            <>
              <div className="query-mini">
                <small>Active query</small>
                <p title={compactQuery}>{compactQuery}</p>
              </div>

              <div className="loading-workbench">
                <div className="loading-chart">
                  <div className="loading-chart-header">
                    <strong>Processing Stream</strong>
                    <small>Records processed over time</small>
                  </div>
                  <div className="loading-line-chart">
                    <svg viewBox="0 0 100 40" preserveAspectRatio="none" aria-hidden="true">
                      <defs>
                        <linearGradient id="lineAreaGrad" x1="0" y1="0" x2="0" y2="1">
                          <stop offset="0%" stopColor="rgba(74, 222, 128, 0.45)" />
                          <stop offset="100%" stopColor="rgba(74, 222, 128, 0.02)" />
                        </linearGradient>
                      </defs>
                      <polygon points={areaPoints} fill="url(#lineAreaGrad)" />
                      <polyline points={linePoints} fill="none" stroke="rgba(74, 222, 128, 0.95)" strokeWidth="1.6" />
                    </svg>
                  </div>
                  <div className="loading-chart-metrics">
                    <article>
                      <small>Processed</small>
                      <strong>{processedCount.toLocaleString()}</strong>
                    </article>
                    <article>
                      <small>Target window</small>
                      <strong>{estimatedTarget.toLocaleString()}</strong>
                    </article>
                  </div>
                </div>

                <div className="loading-skeleton" aria-hidden="true">
                  <div className="skeleton-line w-70" />
                  <div className="skeleton-line w-55" />
                  <div className="skeleton-line w-85" />
                  <div className="skeleton-row">
                    <div className="skeleton-cell w-25" />
                    <div className="skeleton-cell w-45" />
                    <div className="skeleton-cell w-30" />
                  </div>
                  <div className="skeleton-row">
                    <div className="skeleton-cell w-20" />
                    <div className="skeleton-cell w-40" />
                    <div className="skeleton-cell w-35" />
                  </div>
                </div>
              </div>
            </>
          ) : (
            <>
              <div className="query-block query-block-colored">
                {querySegments.map((segment, index) => {
                  if (segment.type === 'plain') {
                    return (
                      <span key={`plain-${index}`} className="query-segment-plain">
                        {segment.text}
                      </span>
                    );
                  }
                  return (
                    <span
                      key={`group-${index}`}
                      className={`query-segment query-segment--${(segment.groupIndex ?? 0) % 6}`}
                    >
                      {segment.text}
                    </span>
                  );
                })}
              </div>
              <div className="actions">
                <button type="button" className="ghost" onClick={handleGenerate} disabled={loading}>
                  <FlaskConical size={16} /> Regenerate
                </button>
                <button type="button" className="ghost" onClick={onDiscard} disabled={loading}>
                  Discard
                </button>
                <button type="button" className="primary" onClick={onConfirm} disabled={loading}>
                  {isRunning ? <Pickaxe size={16} className="chop" /> : <Search size={16} />}
                  {isRunning ? 'Searching...' : 'Confirm and search'}
                </button>
              </div>
            </>
          )}
        </section>
      ) : null}
    </section>
  );
}
