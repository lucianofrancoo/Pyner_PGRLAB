import { useMemo, useState } from 'react';
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

const SOURCE_DESCRIPTIONS: Record<SourceMode, string> = {
  pubmed: 'PubMed: fast, focused on papers and abstracts.',
  pmc: 'PMC: full-text search — more results, slower.',
  bioproject: 'BioProject: slower, focused on projects and SRA hierarchy.',
};

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
  const [source, setSource] = useState<SourceMode>('pubmed');
  const [useLlm, setUseLlm] = useState(true);
  // '' = unlimted input display; internally sent as MAX_HARD_LIMIT
  const [maxResultsInput, setMaxResultsInput] = useState<string>('');

  const resolvedMaxResults = useMemo(() => {
    const parsed = parseInt(maxResultsInput, 10);
    if (!maxResultsInput.trim() || isNaN(parsed) || parsed <= 0) return MAX_HARD_LIMIT;
    return Math.min(parsed, MAX_HARD_LIMIT);
  }, [maxResultsInput]);

  const sourceDescription = SOURCE_DESCRIPTIONS[source];

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

      <form className="search-form" onSubmit={handleSubmit}>
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
              <option value="pubmed">PubMed</option>
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
          <button type="button" className="ghost" onClick={() => setNaturalQuery(DEFAULT_QUERY)}>
            <FlaskConical size={16} /> Load sample
          </button>

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
            <p>Review the generated NCBI query and confirm to run the real search.</p>
          </header>
          <pre className="query-block">{pendingQuery.ncbi_query}</pre>
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
        </section>
      ) : null}
    </section>
  );
}
