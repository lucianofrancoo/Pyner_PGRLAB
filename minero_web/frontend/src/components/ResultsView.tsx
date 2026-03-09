import { useEffect, useMemo, useState } from 'react';
import { Download, FlaskConical, Search as SearchIcon } from 'lucide-react';
import type { BioprojectResult, MineroResponse, MineroResult, PubmedResult } from '../types';
import { exportCsv, exportJson } from '../lib/exporters';
import {
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer,
  PieChart, Pie, Cell, Legend
} from 'recharts';

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

function BioProjectDashboard({ results }: { results: BioprojectResult[] }) {
  let totalExps = 0;
  let totalRuns = 0;
  let totalBiosamples = 0;

  const tissueMap: Record<string, number> = {};
  const conditionMap: Record<string, number> = {};
  const strategyMap: Record<string, number> = {};

  results.forEach(r => {
    totalExps += parseCount(r.sra_experiments_count);
    totalRuns += parseCount(r.sra_runs_count);
    totalBiosamples += parseCount(r.biosamples_count);

    const tags = r.classification.tags ?? [];
    const tissue = extractTagByPrefix(tags, 'tejido:');
    if (tissue !== '-') tissueMap[tissue] = (tissueMap[tissue] || 0) + 1;

    const condition = extractTagByPrefix(tags, 'condicion:');
    if (condition !== '-') conditionMap[condition] = (conditionMap[condition] || 0) + 1;

    const strategy = extractTagByPrefix(tags, 'estrategia:');
    if (strategy !== '-') strategyMap[strategy] = (strategyMap[strategy] || 0) + 1;
  });

  const totalsData = [
    { name: 'Projects', count: results.length },
    { name: 'Biosamples', count: totalBiosamples },
    { name: 'Runs', count: totalRuns },
  ].reverse(); // reverse to show larger categories at the top/bottom depending on layout

  const toChartData = (map: Record<string, number>) =>
    Object.entries(map).map(([name, value]) => ({ name, value })).sort((a, b) => b.value - a.value).slice(0, 10);

  const tissueData = toChartData(tissueMap);
  const conditionData = toChartData(conditionMap);
  const strategyData = toChartData(strategyMap);

  const COLORS = ['#19d3a2', '#3b82f6', '#9b5de5', '#f59e0b', '#10b981', '#ef4444', '#d4a017', '#60a5fa', '#a78bfa', '#f87171'];

  return (
    <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(320px, 1fr))', gap: '16px', marginBottom: '24px' }}>
      <div className="panel flow-panel" style={{ height: 280, padding: '16px', display: 'flex', flexDirection: 'column' }}>
        <h4 style={{ margin: '0 0 10px 0', fontSize: '0.75rem', color: '#8ba5c4', textTransform: 'uppercase', letterSpacing: '0.1em' }}>Volume Overview</h4>
        <div style={{ flex: 1, minHeight: 0 }}>
          <ResponsiveContainer width="100%" height="100%">
            <BarChart data={totalsData} layout="vertical" margin={{ top: 5, right: 30, left: 10, bottom: 5 }}>
              <XAxis type="number" stroke="#4f6a88" fontSize={11} allowDecimals={false} />
              <YAxis dataKey="name" type="category" stroke="#dce8f7" fontSize={11} width={80} fontWeight={600} />
              <Tooltip contentStyle={{ backgroundColor: '#1c2d42', borderColor: '#4f6a88', borderRadius: '8px', color: '#dce8f7', fontSize: '12px' }} cursor={{ fill: 'rgba(255,255,255,0.05)' }} />
              <Bar dataKey="count" fill="#4ade80" radius={[0, 4, 4, 0]} barSize={20} />
            </BarChart>
          </ResponsiveContainer>
        </div>
      </div>

      {tissueData.length > 0 && (
        <div className="panel flow-panel" style={{ height: 280, padding: '16px', display: 'flex', flexDirection: 'column' }}>
          <h4 style={{ margin: '0 0 4px 0', fontSize: '0.75rem', color: '#8ba5c4', textTransform: 'uppercase', letterSpacing: '0.1em' }}>Tissues Detected</h4>
          <div style={{ flex: 1, minHeight: 0 }}>
            <ResponsiveContainer width="100%" height="100%">
              <PieChart margin={{ top: 0, right: 0, left: 0, bottom: 0 }}>
                <Pie data={tissueData} innerRadius="55%" outerRadius="80%" paddingAngle={2} dataKey="value">
                  {tissueData.map((_, index) => <Cell key={`cell-${index}`} fill={COLORS[index % COLORS.length]} />)}
                </Pie>
                <Tooltip contentStyle={{ backgroundColor: '#1c2d42', borderColor: '#4f6a88', borderRadius: '8px', color: '#dce8f7', fontSize: '12px' }} />
                <Legend verticalAlign="bottom" height={24} wrapperStyle={{ fontSize: '11px', color: '#dce8f7', paddingTop: '10px' }} />
              </PieChart>
            </ResponsiveContainer>
          </div>
        </div>
      )}

      {conditionData.length > 0 && (
        <div className="panel flow-panel" style={{ height: 280, padding: '16px', display: 'flex', flexDirection: 'column' }}>
          <h4 style={{ margin: '0 0 10px 0', fontSize: '0.75rem', color: '#8ba5c4', textTransform: 'uppercase', letterSpacing: '0.1em' }}>Top Conditions</h4>
          <div style={{ flex: 1, minHeight: 0 }}>
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={conditionData} margin={{ top: 5, right: 10, left: -20, bottom: 40 }}>
                <XAxis dataKey="name" stroke="#4f6a88" fontSize={10} angle={-35} textAnchor="end" interval={0} />
                <YAxis stroke="#4f6a88" fontSize={11} allowDecimals={false} />
                <Tooltip contentStyle={{ backgroundColor: '#1c2d42', borderColor: '#4f6a88', borderRadius: '8px', color: '#dce8f7', fontSize: '12px' }} cursor={{ fill: 'rgba(255,255,255,0.05)' }} />
                <Bar dataKey="value" fill="#d4a017" radius={[4, 4, 0, 0]} barSize={26} />
              </BarChart>
            </ResponsiveContainer>
          </div>
        </div>
      )}

      {strategyData.length > 0 && (
        <div className="panel flow-panel" style={{ height: 280, padding: '16px', display: 'flex', flexDirection: 'column' }}>
          <h4 style={{ margin: '0 0 4px 0', fontSize: '0.75rem', color: '#8ba5c4', textTransform: 'uppercase', letterSpacing: '0.1em' }}>Strategies</h4>
          <div style={{ flex: 1, minHeight: 0 }}>
            <ResponsiveContainer width="100%" height="100%">
              <PieChart margin={{ top: 0, right: 0, left: 0, bottom: 0 }}>
                <Pie data={strategyData} innerRadius="30%" outerRadius="80%" paddingAngle={1} dataKey="value">
                  {strategyData.map((_, index) => <Cell key={`cell-${index}`} fill={COLORS[(index + 4) % COLORS.length]} />)}
                </Pie>
                <Tooltip contentStyle={{ backgroundColor: '#1c2d42', borderColor: '#4f6a88', borderRadius: '8px', color: '#dce8f7', fontSize: '12px' }} />
                <Legend verticalAlign="bottom" height={24} wrapperStyle={{ fontSize: '11px', color: '#dce8f7', paddingTop: '10px' }} />
              </PieChart>
            </ResponsiveContainer>
          </div>
        </div>
      )}
    </div>
  );
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
      let biosamples = 0;
      let runs = 0;

      rows.forEach((row) => {
        const haystack = `${row.title ?? ''} ${row.description ?? ''}`;
        tissueCounts[inferTissue(haystack)] += 1;
        biosamples += parseCount(row.biosamples_count);
        runs += parseCount(row.sra_runs_count);
      });

      return {
        total: rows.length,
        withDoi: 0,
        withPmcid: 0,
        reviews: 0,
        journals: 0,
        years: [],
        projects: rows.length,
        biosamples,
        runs,
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
              <strong>{isPubmedSource ? streamStats.withDoi : streamStats.biosamples}</strong>
              <span>{isPubmedSource ? 'with DOI' : 'biosamples'}</span>
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
              <strong>{isPubmedSource ? streamStats.journals : filteredResults.length}</strong>
              <span>{isPubmedSource ? 'journals' : 'projects'}</span>
            </article>
            <article>
              <strong>{isPubmedSource ? streamStats.years.length : streamStats.runs}</strong>
              <span>{isPubmedSource ? 'year span' : 'total runs'}</span>
            </article>
          </div>
        </div>
      </section>

      {!isPubmedSource && filteredResults.length > 0 && (
        <BioProjectDashboard results={filteredResults as BioprojectResult[]} />
      )}

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
              <>
                <dl>
                  <div>
                    <dt>Biosamples</dt>
                    <dd>{parseCount(selected.biosamples_count)}</dd>
                  </div>
                  <div>
                    <dt>Runs</dt>
                    <dd>{parseCount(selected.sra_runs_count)}</dd>
                  </div>
                </dl>
              </>
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
