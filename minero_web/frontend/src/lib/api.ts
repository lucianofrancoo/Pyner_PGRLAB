import type { MineroResponse, ProAnalysisResponse, QueryGeneration, SearchPayload } from '../types';

const API_BASE = import.meta.env.VITE_MINERO_API_URL ?? '';
const DEFAULT_TIMEOUT_MS = 180000;
const PUBMED_SEARCH_TIMEOUT_MS = 180000;
const BIOPROJECT_SEARCH_TIMEOUT_MS = 1200000;

async function fetchWithTimeout(
  input: string,
  init: RequestInit = {},
  timeoutMs: number = DEFAULT_TIMEOUT_MS
): Promise<Response> {
  const controller = new AbortController();
  const timeoutId = window.setTimeout(() => controller.abort(), timeoutMs);

  try {
    return await fetch(input, { ...init, signal: controller.signal });
  } catch (error) {
    if (error instanceof DOMException && error.name === 'AbortError') {
      throw new Error(
        `Timeout while waiting for backend response (${API_BASE || 'proxy /api'}).`
      );
    }
    throw error;
  } finally {
    window.clearTimeout(timeoutId);
  }
}

export async function generateQuery(payload: Pick<SearchPayload, 'natural_query' | 'use_llm'>): Promise<QueryGeneration> {
  const response = await fetchWithTimeout(`${API_BASE}/api/minero/generate-query`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload),
  });

  if (!response.ok) {
    let detail = 'Could not generate query';
    try {
      const data = await response.json();
      detail = data.detail ?? detail;
    } catch {
      // ignore
    }
    throw new Error(detail);
  }

  const data = (await response.json()) as { query_generation: QueryGeneration };
  return data.query_generation;
}

export async function runSearchWithQuery(
  payload: Omit<SearchPayload, 'natural_query'> & { query_generation: QueryGeneration; ncbi_query: string }
): Promise<MineroResponse> {
  const timeoutMs =
    payload.source === 'bioproject' ? BIOPROJECT_SEARCH_TIMEOUT_MS : PUBMED_SEARCH_TIMEOUT_MS;

  const response = await fetchWithTimeout(`${API_BASE}/api/minero/run-search`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload),
  }, timeoutMs);

  if (!response.ok) {
    let detail = 'Could not run search';
    try {
      const data = await response.json();
      detail = data.detail ?? detail;
    } catch {
      // ignore
    }
    throw new Error(detail);
  }

  return response.json();
}

export async function runMineroSearch(payload: SearchPayload): Promise<MineroResponse> {
  const timeoutMs =
    payload.source === 'bioproject' ? BIOPROJECT_SEARCH_TIMEOUT_MS : PUBMED_SEARCH_TIMEOUT_MS;

  const response = await fetchWithTimeout(`${API_BASE}/api/minero/search`, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(payload),
  }, timeoutMs);

  if (!response.ok) {
    let detail = 'Could not complete search';
    try {
      const data = await response.json();
      detail = data.detail ?? detail;
    } catch {
      // ignore JSON parsing errors
    }
    throw new Error(detail);
  }

  return response.json();
}

export async function checkHealth(): Promise<{ llm_runtime_available: boolean }> {
  const response = await fetchWithTimeout(`${API_BASE}/api/minero/health`);
  if (!response.ok) {
    throw new Error('Could not verify backend health');
  }
  return response.json();
}

// Modo Pro — timeout generoso: LLM procesa cada paper individualmente
const PRO_ANALYSIS_TIMEOUT_MS = 30 * 60 * 1000; // 30 min

export async function analyzePapers(
  publications: Record<string, unknown>[],
  query: string
): Promise<ProAnalysisResponse> {
  const response = await fetchWithTimeout(
    `${API_BASE}/api/minero/analyze-papers`,
    {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ publications, query }),
    },
    PRO_ANALYSIS_TIMEOUT_MS
  );

  if (!response.ok) {
    let detail = 'Pro analysis failed';
    try {
      const data = await response.json();
      detail = data.detail ?? detail;
    } catch {
      // ignore
    }
    throw new Error(detail);
  }

  return response.json();
}
