import { useEffect, useMemo, useState } from 'react';
import { BookOpenText, ChartSpline, Database, Search } from 'lucide-react';
import { SearchView } from './components/SearchView';
import { ResultsView } from './components/ResultsView';
import { ProAnalysisView } from './components/ProAnalysisView';
import { AnalysisView } from './components/AnalysisView';
import { HelpView } from './components/HelpView';
import { StatusBanner } from './components/StatusBanner';
import { analyzePapers, checkHealth, generateQuery, runSearchWithQuery } from './lib/api';
import type {
  AppStatus,
  AppView,
  MineroResponse,
  ProAnalysisResponse,
  QueryGeneration,
  SearchHistoryEntry,
  SearchPayload,
} from './types';

const API_BASE = import.meta.env.VITE_MINERO_API_URL ?? 'http://127.0.0.1:8010';

// ─── Nombre del proyecto ────────────────────────────────────────────────────
// Cambiar este valor cuando se defina el nombre final (Pyner / MAIner / Minero)
const APP_NAME = 'Minero';
const HISTORY_STORAGE_KEY = 'minero.search_history.v1';
const MAX_HISTORY_ENTRIES = 20;

const NAV_ITEMS: Array<{ id: AppView; label: string; icon: typeof Search }> = [
  { id: 'buscar', label: 'Search', icon: Search },
  { id: 'resultados', label: 'Classified Results', icon: Database },
  { id: 'analisis', label: 'Analytics', icon: ChartSpline },
  { id: 'ayuda', label: 'Help', icon: BookOpenText },
];

const LOADING_STEPS = {
  generating: ['Reading your question...', 'Building search query...'],
  running: ['Searching records...', 'Sorting results...', 'Scoring relevance...', 'Preparing summary...'],
  pro: ['Reading abstracts...', 'Extracting key entities...', 'Creating graphs...', 'Scoring evidence...'],
} as const;

function statusMessage(status: AppStatus): string {
  switch (status) {
    case 'loading':
      return 'Processing query and classifying records...';
    case 'success':
      return 'Search completed with classified results.';
    case 'partial-success':
      return 'Search completed with warnings. Review flagged records.';
    case 'empty':
      return 'No records found for the current filters.';
    case 'error':
      return 'The operation could not be completed.';
    default:
      return 'System ready. Enter a biological query to begin.';
  }
}

interface PendingRun {
  source: SearchPayload['source'];
  max_results: number;
  use_llm: boolean;
  query_generation: QueryGeneration;
}

function PickaxeLogo() {
  return (
    <img
      src="/logos/gemini.png"
      alt="Minero project logo"
      className="pickaxe-logo"
    />
  );
}

/** Logos institucionales: PGRLab y Phytolearning en el footer del sidebar. */
function InstitutionalLogos() {
  return (
    <div className="inst-logos">
      <a
        className="inst-logo-wrap"
        href="https://pgrlab.cl/"
        target="_blank"
        rel="noreferrer"
        aria-label="Visit Plant Genome Regulation Lab website"
      >
        <img
          src="/logos/pgrlab.png"
          alt="PGRLab"
          title="Plant Genome Regulation Lab (PGRLab)"
          className="inst-logo"
          onError={(e) => { (e.target as HTMLImageElement).parentElement!.style.display = 'none'; }}
        />
      </a>
      <a
        className="inst-logo-wrap inst-logo-wrap--phyto"
        href="https://phytolearning.cl/"
        target="_blank"
        rel="noreferrer"
        aria-label="Visit Phytolearning website"
      >
        <img
          src="/logos/phytolearning.png"
          alt="Phytolearning"
          title="Phytolearning — Núcleo Milenio en Ciencia de Datos y Resiliencia Vegetal"
          className="inst-logo"
          onError={(e) => { (e.target as HTMLImageElement).parentElement!.style.display = 'none'; }}
        />
      </a>
    </div>
  );
}

export default function App() {
  const [view, setView] = useState<AppView>('buscar');
  const [status, setStatus] = useState<AppStatus>('idle');
  const [response, setResponse] = useState<MineroResponse | null>(null);
  const [searchHistory, setSearchHistory] = useState<SearchHistoryEntry[]>([]);
  const [error, setError] = useState<string>('');
  const [llmAvailable, setLlmAvailable] = useState(false);
  const [pendingRun, setPendingRun] = useState<PendingRun | null>(null);
  const [step, setStep] = useState<'idle' | 'generating' | 'running'>('idle');
  const [loadingStepIndex, setLoadingStepIndex] = useState(0);

  // ─── Estado modo Pro ──────────────────────────────────────────
  const [proResponse, setProResponse] = useState<ProAnalysisResponse | null>(null);
  const [proLoading, setProLoading] = useState(false);
  const [showPro, setShowPro] = useState(false);

  useEffect(() => {
    checkHealth()
      .then((health) => setLlmAvailable(Boolean(health.llm_runtime_available)))
      .catch(() => setLlmAvailable(false));
  }, []);

  useEffect(() => {
    try {
      const raw = window.localStorage.getItem(HISTORY_STORAGE_KEY);
      if (!raw) return;
      const parsed = JSON.parse(raw);
      if (Array.isArray(parsed)) {
        setSearchHistory(parsed.slice(0, MAX_HISTORY_ENTRIES));
      }
    } catch {
      setSearchHistory([]);
    }
  }, []);

  useEffect(() => {
    window.localStorage.setItem(HISTORY_STORAGE_KEY, JSON.stringify(searchHistory));
  }, [searchHistory]);

  async function handleGenerate(payload: SearchPayload): Promise<void> {
    setStep('generating');
    setStatus('loading');
    setError('');

    try {
      const query_generation = await generateQuery({
        natural_query: payload.natural_query,
        use_llm: payload.use_llm,
      });
      setPendingRun({
        source: payload.source,
        max_results: payload.max_results,
        use_llm: payload.use_llm,
        query_generation,
      });
      setStatus('idle');
      setView('buscar');
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Unexpected error';
      setError(message);
      setStatus('error');
    } finally {
      setStep('idle');
    }
  }

  async function handleConfirm(): Promise<void> {
    if (!pendingRun) return;

    setStep('running');
    setStatus('loading');
    setError('');
    // Reset pro state on new search
    setProResponse(null);
    setShowPro(false);

    try {
      const result = await runSearchWithQuery({
        source: pendingRun.source,
        max_results: pendingRun.max_results,
        use_llm: pendingRun.use_llm,
        ncbi_query: pendingRun.query_generation.ncbi_query,
        query_generation: pendingRun.query_generation,
      });
      setResponse(result);
      const historyEntry: SearchHistoryEntry = {
        id: `${Date.now()}-${Math.random().toString(36).slice(2, 8)}`,
        created_at: new Date().toISOString(),
        source: pendingRun.source,
        natural_query: pendingRun.query_generation.user_input ?? '',
        ncbi_query: pendingRun.query_generation.ncbi_query,
        status: result.metadata.status,
        total_results: result.metadata.total_results,
        response: result,
      };
      setSearchHistory((prev) => [historyEntry, ...prev].slice(0, MAX_HISTORY_ENTRIES));
      setPendingRun(null);
      setStatus(result.metadata.status);
      setView('resultados');
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Unexpected error';
      setError(message);
      setStatus('error');
    } finally {
      setStep('idle');
    }
  }

  const message = useMemo(() => {
    if (status === 'error' && error) {
      return error;
    }
    return statusMessage(status);
  }, [status, error]);

  const activeLoadingStage = useMemo<'generating' | 'running' | 'pro' | null>(() => {
    if (proLoading) return 'pro';
    if (step === 'generating') return 'generating';
    if (step === 'running' || status === 'loading') return 'running';
    return null;
  }, [proLoading, step, status]);

  useEffect(() => {
    setLoadingStepIndex(0);

    if (!activeLoadingStage) return;

    const stages = LOADING_STEPS[activeLoadingStage];
    const interval = window.setInterval(() => {
      setLoadingStepIndex((prev) => (prev + 1) % stages.length);
    }, 1400);

    return () => window.clearInterval(interval);
  }, [activeLoadingStage]);

  const loadingDetail = activeLoadingStage
    ? LOADING_STEPS[activeLoadingStage][loadingStepIndex]
    : null;

  const authorByPmid = useMemo<Record<string, string>>(() => {
    if (!response) return {};
    const map: Record<string, string> = {};
    response.results.forEach((item) => {
      if ('pmid' in item && Array.isArray(item.authors) && item.authors.length > 0) {
        map[item.pmid] = item.authors[0];
      }
    });
    return map;
  }, [response]);

  async function handleRunPro(
    publications: Record<string, unknown>[],
    query: string
  ): Promise<void> {
    setProLoading(true);
    try {
      const result = await analyzePapers(publications, query);
      setProResponse(result);
      setShowPro(true);
    } catch (err) {
      const message = err instanceof Error ? err.message : 'Pro analysis failed';
      setError(message);
      setStatus('error');
    } finally {
      setProLoading(false);
    }
  }

  function handleOpenHistory(entryId: string): void {
    const match = searchHistory.find((entry) => entry.id === entryId);
    if (!match) return;
    setResponse(match.response);
    setPendingRun(null);
    setProResponse(null);
    setShowPro(false);
    setError('');
    setStatus(match.status);
    setView('resultados');
  }

  function handleDeleteHistory(entryId: string): void {
    setSearchHistory((prev) => prev.filter((entry) => entry.id !== entryId));
  }

  function handleClearHistory(): void {
    setSearchHistory([]);
  }

  return (
    <div className="app-shell">
      <div className="bg-layer bg-grid" />
      <div className="bg-layer bg-mesh" />
      <div className="scanline" />

      <aside className="sidebar">
        <header className="sidebar-header">
          <div className="brand-cluster">
            <PickaxeLogo />
            <div>
              <p className="brand-eyebrow">PGRLAB · PHYTOLEARNING</p>
              <h1>{APP_NAME}</h1>
              <span>Guided scientific literature mining</span>
            </div>
          </div>
        </header>

        <nav>
          {NAV_ITEMS.map(({ id, label, icon: Icon }) => (
            <button
              key={id}
              className={view === id ? 'nav-item active' : 'nav-item'}
              onClick={() => setView(id)}
            >
              <Icon size={18} />
              <span>{label}</span>
            </button>
          ))}
        </nav>

        <footer>
          <InstitutionalLogos />
          <small>Build stable 1.0.0</small>
          <small>LLM: {llmAvailable ? 'Ollama available' : 'Heuristic fallback'}</small>
          <small>API: {API_BASE}</small>
        </footer>
      </aside>

      <main className="main-content">
        <div className="main-frame">
          <StatusBanner
            status={status}
            message={message}
            detail={loadingDetail}
            metadata={response?.metadata ?? null}
          />

          <section className="workspace">
            {view === 'buscar' ? (
              <SearchView
                loading={status === 'loading'}
                isGenerating={step === 'generating'}
                isRunning={step === 'running'}
                llmAvailable={llmAvailable}
                pendingQuery={pendingRun?.query_generation ?? null}
                searchHistory={searchHistory}
                onGenerate={handleGenerate}
                onConfirm={handleConfirm}
                onDiscard={() => setPendingRun(null)}
                onOpenHistory={handleOpenHistory}
                onDeleteHistory={handleDeleteHistory}
                onClearHistory={handleClearHistory}
              />
            ) : null}
            {view === 'resultados' ? (
              showPro && proResponse ? (
                <ProAnalysisView
                  proResponse={proResponse}
                  authorByPmid={authorByPmid}
                  onClose={() => setShowPro(false)}
                />
              ) : (
                <ResultsView
                  response={response}
                  onRunPro={handleRunPro}
                  proLoading={proLoading}
                />
              )
            ) : null}
            {view === 'analisis' ? (
              <AnalysisView response={response} proResponse={proResponse} />
            ) : null}
            {view === 'ayuda' ? <HelpView /> : null}
          </section>
        </div>
      </main>
    </div>
  );
}
