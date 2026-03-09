import { useMemo, useState } from 'react';
import { Download, Search as SearchIcon, Network, Table as TableIcon } from 'lucide-react';
import type { ProAnalysisResponse, ProAnalysisResult } from '../types';

interface ProAnalysisViewProps {
    proResponse: ProAnalysisResponse;
    onClose: () => void;
}

// Grupos de columnas para mostrar en el detalle
const FIELD_GROUPS: Array<{ label: string; fields: (keyof ProAnalysisResult)[] }> = [
    {
        label: 'Organism',
        fields: ['Organisms', 'Species', 'Strain_Variety', 'Genotype'],
    },
    {
        label: 'Sample',
        fields: ['Tissues_Organs', 'Source_Tissue_Origin', 'Cell_Type', 'Developmental_Stage', 'Organism_Age', 'Growth_Phase'],
    },
    {
        label: 'Conditions',
        fields: ['Conditions', 'Environmental_Stress', 'Temperature_Range', 'Light_Conditions', 'Growth_Medium', 'Sample_Collection_Conditions'],
    },
    {
        label: 'Molecules',
        fields: ['Molecules_Extracted', 'RNA_Type', 'DNA_Type', 'Protein_Type', 'Other_Molecules'],
    },
    {
        label: 'Techniques',
        fields: ['Strategies', 'Measurement_Tools', 'Detection_Method'],
    },
    {
        label: 'Time course',
        fields: ['Time_Course_Design', 'Time_Points', 'Time_Intervals', 'Time_Duration'],
    },
    {
        label: 'Replication',
        fields: ['Sample_Size', 'Biological_Replicates', 'Technical_Replicates', 'Replication_Design'],
    },
    {
        label: 'Treatment',
        fields: ['Treatment_Groups', 'Control_Type', 'Dose_Range'],
    },
    {
        label: 'Quality',
        fields: ['Quality_Metrics', 'Contamination_Check'],
    },
    {
        label: 'Biology',
        fields: ['Pathway_Focus', 'Biomarkers_Measured', 'Disease_Model'],
    },
    {
        label: 'Analysis',
        fields: ['Normalization_Method', 'Statistical_Method', 'Differential_Expression_Threshold', 'Raw_Data_Available'],
    },
];

function scoreColor(score: number): string {
    if (score >= 7) return '#19d3a2';
    if (score >= 4) return '#f59e0b';
    return '#ef4444';
}

function exportProCsv(results: ProAnalysisResult[]): void {
    if (!results.length) return;
    const cols = Object.keys(results[0]) as (keyof ProAnalysisResult)[];
    const header = cols.join('\t');
    const rows = results.map((r) =>
        cols.map((c) => {
            const val = r[c] ?? '';
            return String(val).replace(/\t/g, ' ');
        }).join('\t')
    );
    const tsv = [header, ...rows].join('\n');
    const blob = new Blob([tsv], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `pro_analysis_${Date.now()}.tsv`;
    a.click();
    URL.revokeObjectURL(url);
}

export function ProAnalysisView({ proResponse, onClose }: ProAnalysisViewProps) {
    const [searchText, setSearchText] = useState('');
    const [selectedIndex, setSelectedIndex] = useState(0);
    const [activeGroup, setActiveGroup] = useState(0);
    const [viewMode, setViewMode] = useState<'table' | 'network'>('table');

    const filtered = useMemo(() => {
        if (!searchText.trim()) return proResponse.results;
        const q = searchText.toLowerCase();
        return proResponse.results.filter((r) =>
            JSON.stringify(r).toLowerCase().includes(q)
        );
    }, [proResponse.results, searchText]);

    const selected: ProAnalysisResult | null = filtered[selectedIndex] ?? null;
    const meta = proResponse.metadata;

    return (
        <section className="panel repo-panel pro-panel">
            <header className="repo-header">
                <div className="repo-title-block">
                    <h2>Pro Analysis</h2>
                    <p>53 EXPERIMENTAL METADATA COLUMNS · OLLAMA {meta.model.toUpperCase()}</p>
                </div>
                <div className="repo-tools">
                    <div className="toggle-group" style={{ display: 'flex', gap: '0.25rem', background: '#0f172a', padding: '0.25rem', borderRadius: '4px', border: '1px solid #1e293b' }}>
                        <button
                            type="button"
                            className={`ghost ${viewMode === 'table' ? 'active' : ''}`}
                            onClick={() => setViewMode('table')}
                            style={{
                                background: viewMode === 'table' ? '#1e293b' : 'transparent',
                                padding: '0.25rem 0.75rem',
                                color: viewMode === 'table' ? '#f8fafc' : '#94a3b8'
                            }}
                        >
                            <TableIcon size={14} /> Data
                        </button>
                        {proResponse.network_html && (
                            <button
                                type="button"
                                className={`ghost ${viewMode === 'network' ? 'active' : ''}`}
                                onClick={() => setViewMode('network')}
                                style={{
                                    background: viewMode === 'network' ? '#1e293b' : 'transparent',
                                    padding: '0.25rem 0.75rem',
                                    color: viewMode === 'network' ? '#f8fafc' : '#94a3b8'
                                }}
                            >
                                <Network size={14} /> Graph
                            </button>
                        )}
                    </div>
                    {viewMode === 'table' && (
                        <label className="repo-filter">
                            <SearchIcon size={15} />
                            <input
                                value={searchText}
                                onChange={(e) => setSearchText(e.target.value)}
                                placeholder="Quick filter..."
                            />
                        </label>
                    )}
                    <button
                        className="primary repo-export"
                        type="button"
                        onClick={() => exportProCsv(filtered)}
                    >
                        <Download size={14} /> Export TSV
                    </button>
                    <button className="ghost" type="button" onClick={onClose}>
                        ← Back to Results
                    </button>
                </div>
            </header>

            {/* Stats strip */}
            <section className="repo-metadata-stream">
                <div className="stream-column">
                    <p>PRO SUMMARY</p>
                    <div className="stream-kpis">
                        <article>
                            <strong>{meta.total_analyzed}</strong>
                            <span>analyzed</span>
                        </article>
                        <article>
                            <strong style={{ color: '#19d3a2' }}>{meta.total_relevant}</strong>
                            <span>relevant</span>
                        </article>
                        <article>
                            <strong style={{ color: '#ef4444' }}>{meta.total_errors}</strong>
                            <span>errors</span>
                        </article>
                        <article>
                            <strong>{meta.pmc_full_text_used}</strong>
                            <span>full text</span>
                        </article>
                    </div>
                </div>
                <div className="stream-column">
                    <p>ENHANCEMENTS</p>
                    <div className="stream-kpis">
                        <article>
                            <strong>{meta.techniques_enhanced}</strong>
                            <span>techniques+</span>
                        </article>
                        <article>
                            <strong>{meta.abstract_only}</strong>
                            <span>abstract only</span>
                        </article>
                    </div>
                </div>
            </section>

            {viewMode === 'network' && proResponse.network_html ? (
                <div style={{ flex: 1, position: 'relative', width: '100%', minHeight: '600px', backgroundColor: '#fff', borderRadius: '4px', overflow: 'hidden' }}>
                    <iframe
                        title="Interactive Network Graph"
                        srcDoc={proResponse.network_html}
                        style={{ width: '100%', height: '100%', border: 'none' }}
                        sandbox="allow-scripts allow-downloads allow-same-origin"
                    />
                </div>
            ) : (
                <>
                    {/* Table */}
                    <div className="repo-table-wrap">
                        <table className="repo-table">
                            <thead>
                                <tr>
                                    <th>PMID</th>
                                    <th>SCORE</th>
                                    <th>TITLE</th>
                                    <th>ORGANISMS</th>
                                    <th>STRATEGIES</th>
                                    <th>CONDITIONS</th>
                                    <th>RELEVANT</th>
                                </tr>
                            </thead>
                            <tbody>
                                {filtered.map((item, idx) => (
                                    <tr
                                        key={`${idx}-${item.PMID}`}
                                        className={idx === selectedIndex ? 'active' : ''}
                                        onClick={() => { setSelectedIndex(idx); setActiveGroup(0); }}
                                    >
                                        <td className="repo-id">{item.PMID}</td>
                                        <td>
                                            <span style={{ color: scoreColor(item.Relevance_Score), fontWeight: 700 }}>
                                                {item.Relevance_Score}/10
                                            </span>
                                        </td>
                                        <td className="repo-title">
                                            {item.Title.length > 60 ? item.Title.slice(0, 59) + '…' : item.Title}
                                        </td>
                                        <td className="repo-organism">
                                            {item.Organisms.length > 28 ? item.Organisms.slice(0, 27) + '…' : item.Organisms}
                                        </td>
                                        <td>
                                            {item.Strategies.length > 22 ? item.Strategies.slice(0, 21) + '…' : item.Strategies}
                                        </td>
                                        <td>
                                            {item.Conditions.length > 22 ? item.Conditions.slice(0, 21) + '…' : item.Conditions}
                                        </td>
                                        <td>
                                            <span className={`repo-chip ${item.Is_Relevant === 'Yes' ? 'chip-yes' : item.Is_Relevant === 'Error' ? 'chip-error' : ''}`}>
                                                {item.Is_Relevant.toUpperCase()}
                                            </span>
                                        </td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>

                    {/* Detail panel */}
                    <aside className="details-panel pro-details">
                        {selected ? (
                            <>
                                <h3>{selected.Title}</h3>
                                <p className="pro-summary">{selected.Summary}</p>
                                <p className="pro-explanation"><em>{selected.Relevance_Explanation}</em></p>

                                {/* Group tabs */}
                                <div className="pro-group-tabs">
                                    {FIELD_GROUPS.map((g, gi) => (
                                        <button
                                            key={g.label}
                                            className={gi === activeGroup ? 'pro-tab active' : 'pro-tab'}
                                            onClick={() => setActiveGroup(gi)}
                                        >
                                            {g.label}
                                        </button>
                                    ))}
                                </div>

                                {/* Active group fields */}
                                <dl className="pro-fields">
                                    {FIELD_GROUPS[activeGroup].fields.map((field) => {
                                        const val = selected[field];
                                        const strVal = String(val ?? '-');
                                        const isEmpty = !strVal || strVal === 'not described' || strVal === '-';
                                        return (
                                            <div key={field} className={isEmpty ? 'pro-field empty' : 'pro-field'}>
                                                <dt>{String(field).replace(/_/g, ' ')}</dt>
                                                <dd>{isEmpty ? <span className="muted">not described</span> : strVal}</dd>
                                            </div>
                                        );
                                    })}
                                </dl>

                                <p className="pro-abstract-preview">{selected.Abstract_Preview}</p>
                            </>
                        ) : (
                            <p>No records match the current filter.</p>
                        )}
                    </aside>
                </>
            )}
        </section>
    );
}
