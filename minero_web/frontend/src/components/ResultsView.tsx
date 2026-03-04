import { useEffect, useMemo, useState } from 'react';
import { Download, FlaskConical, Search as SearchIcon } from 'lucide-react';
import type { BioprojectResult, MineroResponse, MineroResult, PubmedResult } from '../types';
import { exportCsv, exportJson } from '../lib/exporters';

interface ResultsViewProps {
  response: MineroResponse | null;
  onRunPro?: (publications: Record<string, unknown>[], query: string) => void;
  proLoading?: boolean;
}

function isPubmed(source: string, item: MineroResult): item is PubmedResult {
  return source === 'pubmed' && 'pmid' in item;
}

function isBioproject(source: string, item: MineroResult): item is BioprojectResult {
  return source === 'bioproject' && 'bioproject' in item;
}

function shortText(value: string, max: number = 54): string {
  if (value.length <= max) return value;
  return `${value.slice(0, max - 1)}…`;
}

function extractTagByPrefix(tags: string[], prefix: string): string {
  const match = tags.find((tag) => tag.startsWith(prefix));
  if (!match) return '-';
  return match.replace(prefix, '').trim() || '-';
}

function formatTag(tag: string): string {
  return tag
    .replace(/^organismo:/, 'organism: ')
    .replace(/^condicion:/, 'condition: ')
    .replace(/^estrategia:/, 'strategy: ')
    .replace(/^tejido:/, 'tissue: ')
    .replace(/^evidencia:/, 'evidence: ')
    .replace(/^alineacion:/, 'alignment: ');
}

function pubmedTypeLabel(publicationType: string): string {
  if (publicationType.toLowerCase().includes('review')) return 'Review';
  return 'Article';
}

function shortAuthors(authors: string[] | undefined, max: number = 3): string {
  if (!authors || authors.length === 0) return '-';
  if (authors.length <= max) return authors.join(', ');
  return `${authors.slice(0, max).join(', ')} +${authors.length - max}`;
}

function parseCount(value: unknown): number {
  const n = Number(value);
  return Number.isFinite(n) ? n : 0;
}

type TissueLabel = 'root' | 'leaf' | 'seed' | 'fruit' | 'whole-plant' | 'unknown';

const TISSUE_COLORS: Record<TissueLabel, string> = {
  root: '#19d3a2',
  leaf: '#3b82f6',
  seed: '#9b5de5',
  fruit: '#f59e0b',
  'whole-plant': '#10b981',
  unknown: '#334155',
};

function inferTissue(text: string): TissueLabel {
  const t = text.toLowerCase();
  if (/\broot|roots|radic/i.test(t)) return 'root';
  if (/\bleaf|leaves|foliar/i.test(t)) return 'leaf';
  if (/\bseed|seedling/i.test(t)) return 'seed';
  if (/\bfruit|fruits|tomato fruit/i.test(t)) return 'fruit';
  if (/\bplant|whole plant|seedling\b/i.test(t)) return 'whole-plant';
  return 'unknown';
}

export function ResultsView({ response, onRunPro, proLoading }: ResultsViewProps) {
  const [searchText, setSearchText] = useState('');
  const [selectedIndex, setSelectedIndex] = useState(0);

  useEffect(() => {
    setSelectedIndex(0);
  }, [response]);

  const filteredResults = useMemo(() => {
    if (!response) return [];
    return response.results.filter((item) => {
      const haystack = JSON.stringify(item).toLowerCase();
      return searchText.trim() === '' || haystack.includes(searchText.toLowerCase());
    });
  }, [response, searchText]);

  const selected = filteredResults[selectedIndex] ?? null;
  const isPubmedSource =
    response != null &&
    (response.metadata.source === 'pubmed' || response.metadata.source === 'pmc');


  const streamStats = useMemo(() => {
    if (!response) {
      return {
        total: 0,
        withDoi: 0,
        withPmcid: 0,
        reviews: 0,
        journals: 0,
        years: [] as string[],
        projects: 0,
        experiments: 0,
        biosamples: 0,
        runs: 0,
        withLinkedPapers: 0,
        tissueCounts: {} as Record<TissueLabel, number>,
      };
    }

    if (response.metadata.source !== 'pubmed' && response.metadata.source !== 'pmc') {
      const rows = response.results.filter((item): item is BioprojectResult => 'bioproject' in item);
      const tissueCounts: Record<TissueLabel, number> = {
        root: 0,
        leaf: 0,
        seed: 0,
        fruit: 0,
        'whole-plant': 0,
        unknown: 0,
      };
      let experiments = 0;
      let biosamples = 0;
      let runs = 0;
      let withLinkedPapers = 0;

      rows.forEach((row) => {
        const haystack = `${row.title ?? ''} ${row.description ?? ''}`;
        tissueCounts[inferTissue(haystack)] += 1;
        experiments += parseCount(row.sra_experiments_count);
        biosamples += parseCount(row.biosamples_count);
        runs += parseCount(row.sra_runs_count);
        if (parseCount(row.publications_found) > 0) withLinkedPapers += 1;
      });

      return {
        total: rows.length,
        withDoi: 0,
        withPmcid: 0,
        reviews: 0,
        journals: 0,
        years: [],
        projects: rows.length,
        experiments,
        biosamples,
        runs,
        withLinkedPapers,
        tissueCounts,
      };
    }

    const rows = response.results.filter((item): item is PubmedResult => 'pmid' in item);
    const journals = new Set<string>();
    const years = new Set<string>();
    const tissueCounts: Record<TissueLabel, number> = {
      root: 0,
      leaf: 0,
      seed: 0,
      fruit: 0,
      'whole-plant': 0,
      unknown: 0,
    };
    let withDoi = 0;
    let withPmcid = 0;
    let reviews = 0;

    rows.forEach((row) => {
      const doi = String(row.doi ?? '').trim().toUpperCase();
      const pmcid = String(row.pmcid ?? '').trim().toUpperCase();
      const publicationType = String(row.publication_type ?? '').toLowerCase();
      if (doi && doi !== 'NA') withDoi += 1;
      if (pmcid && pmcid !== 'NA') withPmcid += 1;
      if (publicationType.includes('review')) reviews += 1;
      if (row.journal) journals.add(row.journal);
      if (row.year) years.add(row.year);
      tissueCounts[inferTissue(`${row.title ?? ''} ${row.abstract ?? ''}`)] += 1;
    });

    return {
      total: rows.length,
      withDoi,
      withPmcid,
      reviews,
      journals: journals.size,
      years: Array.from(years).sort((a, b) => Number(a) - Number(b)),
      projects: 0,
      experiments: 0,
      biosamples: 0,
      runs: 0,
      withLinkedPapers: 0,
      tissueCounts,
    };
  }, [response]);

  const tissueEntries = useMemo(() => {
    const ordered = Object.entries(streamStats.tissueCounts || {})
      .filter(([, count]) => count > 0)
      .sort((a, b) => b[1] - a[1]) as Array<[TissueLabel, number]>;
    return ordered.length ? ordered : ([['unknown', streamStats.total]] as Array<[TissueLabel, number]>);
  }, [streamStats]);

  const tissueDonutStyle = useMemo(() => {
    const total = Math.max(1, streamStats.total);
    let start = 0;
    const slices = tissueEntries.map(([label, count]) => {
      const pct = (count / total) * 100;
      const end = start + pct;
      const slice = `${TISSUE_COLORS[label]} ${start}% ${end}%`;
      start = end;
      return slice;
    });
    if (start < 100) {
      slices.push(`#334155 ${start}% 100%`);
    }
    return `conic-gradient(${slices.join(', ')})`;
  }, [tissueEntries, streamStats.total]);

  if (!response) {
    return (
      <section className="panel">
        <header className="panel-header">
          <h2>Classified Results</h2>
          <p>Run a search to inspect relevance-ranked scientific records.</p>
        </header>
      </section>
    );
  }

  return (
    <section className="panel repo-panel">
      <header className="repo-header">
        <div className="repo-title-block">
          <h2>Mined Repository</h2>
          <p>SYNCHRONIZED SRA METADATA</p>
        </div>
        <div className="repo-tools">
          <label className="repo-filter">
            <SearchIcon size={15} />
            <input
              value={searchText}
              onChange={(event) => setSearchText(event.target.value)}
              placeholder="Quick Filter..."
            />
          </label>
          <button className="primary repo-export" type="button" onClick={() => exportCsv(filteredResults, response.metadata.source)}>
            <Download size={14} /> Export CSV
          </button>
          <button className="ghost repo-export-secondary" type="button" onClick={() => exportJson(response)}>
            <Download size={14} /> JSON
          </button>
          {/* Botón Pro — solo para PubMed / PMC */}
          {onRunPro && isPubmedSource && (
            <button
              className="ghost repo-export-secondary"
              type="button"
              disabled={proLoading}
              onClick={() => {
                const pubs = (response?.results ?? []) as unknown as Record<string, unknown>[];
                const query = response?.query_generation?.user_input ?? '';
                onRunPro(pubs, query);
              }}
            >
              <FlaskConical size={14} />
              {proLoading ? 'Analyzing…' : 'Run Pro Analysis'}
            </button>
          )}
        </div>
      </header>

      <section className="repo-metadata-stream">
        <div className="stream-column">
          <p>{isPubmedSource ? 'PUBMED / PMC SNAPSHOT' : 'BIOPROJECT SNAPSHOT'}</p>
          <div className="stream-kpis">
            <article>
              <strong>{streamStats.total}</strong>
              <span>{isPubmedSource ? 'papers' : 'projects'}</span>
            </article>
            <article>
              <strong>{isPubmedSource ? streamStats.withDoi : streamStats.experiments}</strong>
              <span>{isPubmedSource ? 'with DOI' : 'experiments'}</span>
            </article>
            <article>
              <strong>{isPubmedSource ? streamStats.withPmcid : streamStats.runs}</strong>
              <span>{isPubmedSource ? 'with PMCID' : 'runs'}</span>
            </article>
            <article>
              <strong>{isPubmedSource ? streamStats.reviews : streamStats.biosamples}</strong>
              <span>{isPubmedSource ? 'reviews' : 'biosamples'}</span>
            </article>
          </div>
        </div>
        <div className="stream-column">
          <p>TISSUE EXPLORED</p>
          <div className="tissue-donut-wrap">
            <div className="tissue-donut" style={{ background: tissueDonutStyle }} />
            <div className="tissue-legend">
              {tissueEntries.slice(0, 5).map(([label, count]) => (
                <div key={label}>
                  <span style={{ background: TISSUE_COLORS[label] }} />
                  <small>{label}</small>
                  <strong>{count}</strong>
                </div>
              ))}
            </div>
          </div>
        </div>
        <div className="stream-column">
          <p>{isPubmedSource ? 'SOURCE QUALITY' : 'PROJECT CONTEXT'}</p>
          <div className="stream-bars">
            <article>
              <strong>{isPubmedSource ? streamStats.journals : streamStats.withLinkedPapers}</strong>
              <span>{isPubmedSource ? 'journals' : 'linked papers'}</span>
            </article>
            <article>
              <strong>{isPubmedSource ? streamStats.years.length : filteredResults.length}</strong>
              <span>{isPubmedSource ? 'year span' : 'filtered'}</span>
            </article>
          </div>
        </div>
      </section>

      <div className="repo-table-wrap">
        <table className="repo-table">
          <thead>
            {isPubmedSource ? (
              <tr>
                <th>PUBMED ID</th>
                <th>YEAR</th>
                <th>TITLE / CONTEXT</th>
                <th>AUTHORS</th>
                <th>JOURNAL</th>
                <th>DOI</th>
                <th>PMCID</th>
                <th>TYPE</th>
              </tr>
            ) : (
              <tr>
                <th>BIOPROJECT ID</th>
                <th>TITLE / CONTEXT</th>
                <th>ORGANISM</th>
                <th>TISSUE</th>
                <th>EXPS</th>
                <th>RUNS</th>
                <th>STRATEGY</th>
              </tr>
            )}
          </thead>
          <tbody>
            {filteredResults.map((item, index) => {
              const active = index === selectedIndex;
              const tags = item.classification.tags ?? [];
              const tissueTag = extractTagByPrefix(tags, 'tejido:');
              const strategyTag = extractTagByPrefix(tags, 'estrategia:');
              if (isPubmed(response.metadata.source, item)) {
                return (
                  <tr
                    key={`${index}-${item.pmid}`}
                    className={active ? 'active' : ''}
                    onClick={() => setSelectedIndex(index)}
                  >
                    <td className="repo-id">{item.pmid}</td>
                    <td>{item.year || '-'}</td>
                    <td className="repo-title">{shortText(item.title || 'Untitled record', 66)}</td>
                    <td className="repo-authors">{shortText(shortAuthors(item.authors), 44)}</td>
                    <td className="repo-organism">{shortText(item.journal || 'Unknown journal', 34)}</td>
                    <td className="repo-doi">{item.doi && item.doi !== 'NA' ? item.doi : '-'}</td>
                    <td>{item.pmcid && item.pmcid !== 'NA' ? item.pmcid : '-'}</td>
                    <td>
                      <span className="repo-chip">{pubmedTypeLabel(item.publication_type ?? 'Journal Article').toUpperCase()}</span>
                    </td>
                  </tr>
                );
              }

              return (
                <tr
                  key={`${index}-${isBioproject(response.metadata.source, item) ? item.bioproject : index}`}
                  className={active ? 'active' : ''}
                  onClick={() => setSelectedIndex(index)}
                >
                  <td className="repo-id">
                    {isBioproject(response.metadata.source, item) ? item.bioproject : '-'}
                  </td>
                  <td className="repo-title">{shortText(item.title || 'Untitled record')}</td>
                  <td className="repo-organism">
                    {isBioproject(response.metadata.source, item) ? item.organism || 'Unknown organism' : '-'}
                  </td>
                  <td>
                    <span className="repo-chip">{tissueTag.toUpperCase()}</span>
                  </td>
                  <td>{isBioproject(response.metadata.source, item) ? parseCount(item.sra_experiments_count) : '-'}</td>
                  <td>{isBioproject(response.metadata.source, item) ? parseCount(item.sra_runs_count) : '-'}</td>
                  <td>{strategyTag === '-' ? '-' : shortText(strategyTag, 20)}</td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      <aside className="details-panel">
        <h3>Why this record matters</h3>
        {selected ? (
          <>
            <p>{selected.classification.reason_short}</p>
            {isPubmed(response.metadata.source, selected) ? (
              <dl>
                <div>
                  <dt>DOI</dt>
                  <dd>{selected.doi && selected.doi !== 'NA' ? selected.doi : '-'}</dd>
                </div>
                <div>
                  <dt>PMCID</dt>
                  <dd>{selected.pmcid && selected.pmcid !== 'NA' ? selected.pmcid : '-'}</dd>
                </div>
                <div>
                  <dt>PubMed URL</dt>
                  <dd>{selected.url ? <a className="repo-link" href={selected.url} target="_blank" rel="noreferrer">Open</a> : '-'}</dd>
                </div>
              </dl>
            ) : isBioproject(response.metadata.source, selected) ? (
              <dl>
                <div>
                  <dt>Experiments</dt>
                  <dd>{parseCount(selected.sra_experiments_count)}</dd>
                </div>
                <div>
                  <dt>Biosamples</dt>
                  <dd>{parseCount(selected.biosamples_count)}</dd>
                </div>
                <div>
                  <dt>Runs</dt>
                  <dd>{parseCount(selected.sra_runs_count)}</dd>
                </div>
              </dl>
            ) : null}
            <ul>
              {selected.classification.tags.map((tag) => (
                <li key={tag}>{formatTag(tag)}</li>
              ))}
            </ul>
            <dl>
              <div>
                <dt>Model source</dt>
                <dd>{selected.classification.model_source}</dd>
              </div>
            </dl>
          </>
        ) : (
          <p>No records match the current filters.</p>
        )}
      </aside>
    </section>
  );
}
