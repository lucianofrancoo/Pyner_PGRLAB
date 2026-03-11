import { useEffect, useMemo, useRef, useState } from 'react';
import type { PointerEvent as ReactPointerEvent } from 'react';
import type { ProAnalysisResult } from '../types';

type EdgeType = 'journal' | 'organism';

interface PaperNetworkProps {
  papers: ProAnalysisResult[];
  authorByPmid?: Record<string, string>;
  title?: string;
  subtitle?: string;
  compact?: boolean;
}

interface NetworkNode {
  id: string;
  pmid: string;
  title: string;
  label: string;
  year: number | null;
  yearRaw: string;
  journal: string;
  score: number;
  isRelevant: boolean;
  organisms: string;
}

interface NetworkEdge {
  source: string;
  target: string;
  type: EdgeType;
}

interface NodeState {
  x: number;
  y: number;
  vx: number;
  vy: number;
}

function shortTitle(value: string, words = 4): string {
  const parts = value.trim().split(/\s+/).filter(Boolean);
  if (parts.length <= words) return value;
  return `${parts.slice(0, words).join(' ')}…`;
}

function parseYear(raw?: string): number | null {
  if (!raw) return null;
  const m = raw.match(/\d{4}/);
  if (!m) return null;
  const n = Number(m[0]);
  if (!Number.isFinite(n)) return null;
  if (n < 1900 || n > 2100) return null;
  return n;
}

function scoreColor(score: number): string {
  if (score >= 7) return '#19d3a2';
  if (score >= 4) return '#d4a017';
  return '#ef4444';
}

function organismTokens(raw: string): string[] {
  const cleaned = raw
    .toLowerCase()
    .replace(/not described|n\/a|none/gi, '')
    .trim();
  if (!cleaned) return [];
  return cleaned
    .split(/[;,|]/)
    .map((v) => v.trim())
    .filter((v) => v.length >= 3);
}

function clamp(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}

export function PaperNetwork({
  papers,
  authorByPmid,
  title = 'Paper Network',
  subtitle = 'Integrated graph view of relevance, journal and publication year.',
  compact = false,
}: PaperNetworkProps) {
  const baseNodes = useMemo<NetworkNode[]>(() => {
    return papers.map((paper, idx) => {
      const pmid = paper.PMID || `paper-${idx + 1}`;
      const year = parseYear(paper.Year);
      const journal = paper.Journal?.trim() || 'Unknown';
      const score = Number.isFinite(Number(paper.Relevance_Score))
        ? Number(paper.Relevance_Score)
        : 0;

      const firstAuthor = (authorByPmid?.[pmid] || '').trim();
      const fallbackLabel = `PMID ${pmid}`;
      const titleLabel = shortTitle(paper.Title || fallbackLabel, 3);

      return {
        id: pmid,
        pmid,
        title: paper.Title || 'Untitled',
        label: firstAuthor || titleLabel,
        year,
        yearRaw: paper.Year || 'N/A',
        journal,
        score,
        isRelevant: String(paper.Is_Relevant).toLowerCase() === 'yes',
        organisms: paper.Organisms || '',
      };
    });
  }, [papers, authorByPmid]);

  const years = useMemo(() => {
    const unique = Array.from(new Set(baseNodes.map((n) => n.year).filter((v): v is number => v !== null)));
    unique.sort((a, b) => a - b);
    return unique;
  }, [baseNodes]);

  const journals = useMemo(() => {
    return Array.from(new Set(baseNodes.map((n) => n.journal))).sort((a, b) => a.localeCompare(b));
  }, [baseNodes]);

  const [minScore, setMinScore] = useState(0);
  const [selectedJournal, setSelectedJournal] = useState('all');
  const [yearMin, setYearMin] = useState<number | null>(null);
  const [yearMax, setYearMax] = useState<number | null>(null);
  const [selectedNodeId, setSelectedNodeId] = useState<string>('');

  const svgRef = useRef<SVGSVGElement | null>(null);
  const draggingIdRef = useRef<string | null>(null);

  const [positions, setPositions] = useState<Record<string, NodeState>>({});

  useEffect(() => {
    if (!years.length) {
      setYearMin(null);
      setYearMax(null);
      return;
    }
    setYearMin((prev) => (prev === null ? years[0] : prev));
    setYearMax((prev) => (prev === null ? years[years.length - 1] : prev));
  }, [years]);

  const filteredNodes = useMemo(() => {
    return baseNodes.filter((node) => {
      if (node.score < minScore) return false;
      if (selectedJournal !== 'all' && node.journal !== selectedJournal) return false;

      if (node.year !== null && yearMin !== null && node.year < yearMin) return false;
      if (node.year !== null && yearMax !== null && node.year > yearMax) return false;
      return true;
    });
  }, [baseNodes, minScore, selectedJournal, yearMin, yearMax]);

  const edges = useMemo<NetworkEdge[]>(() => {
    const ids = new Set(filteredNodes.map((n) => n.id));
    const allEdges: NetworkEdge[] = [];
    const seen = new Set<string>();

    const byJournal = new Map<string, string[]>();
    for (const node of filteredNodes) {
      const key = node.journal.toLowerCase();
      const current = byJournal.get(key) || [];
      current.push(node.id);
      byJournal.set(key, current);
    }

    for (const [, nodeIds] of byJournal) {
      if (nodeIds.length < 2 || nodeIds.length > 12) continue;
      for (let i = 0; i < nodeIds.length; i += 1) {
        for (let j = i + 1; j < nodeIds.length; j += 1) {
          const a = nodeIds[i];
          const b = nodeIds[j];
          const key = a < b ? `${a}::${b}` : `${b}::${a}`;
          if (seen.has(key)) continue;
          seen.add(key);
          allEdges.push({ source: a, target: b, type: 'journal' });
        }
      }
    }

    const byOrganism = new Map<string, string[]>();
    for (const node of filteredNodes) {
      for (const token of organismTokens(node.organisms)) {
        const current = byOrganism.get(token) || [];
        current.push(node.id);
        byOrganism.set(token, current);
      }
    }

    for (const [, nodeIds] of byOrganism) {
      if (nodeIds.length < 2 || nodeIds.length > 8) continue;
      for (let i = 0; i < nodeIds.length; i += 1) {
        for (let j = i + 1; j < nodeIds.length; j += 1) {
          const a = nodeIds[i];
          const b = nodeIds[j];
          const key = a < b ? `${a}::${b}` : `${b}::${a}`;
          if (seen.has(key)) continue;
          seen.add(key);
          allEdges.push({ source: a, target: b, type: 'organism' });
        }
      }
    }

    return allEdges.slice(0, 320).filter((edge) => ids.has(edge.source) && ids.has(edge.target));
  }, [filteredNodes]);

  const canvasSize = useMemo(() => {
    const n = Math.max(filteredNodes.length, 1);
    const side = Math.ceil(Math.sqrt(n));
    const width = Math.max(compact ? 860 : 980, side * (compact ? 170 : 190));
    const height = Math.max(compact ? 360 : 470, side * (compact ? 130 : 150));
    return { width, height };
  }, [filteredNodes.length, compact]);

  useEffect(() => {
    setPositions((prev) => {
      const next: Record<string, NodeState> = {};
      const pad = 40;
      for (const node of filteredNodes) {
        const existing = prev[node.id];
        if (existing) {
          next[node.id] = existing;
          continue;
        }
        next[node.id] = {
          x: pad + Math.random() * Math.max(20, canvasSize.width - pad * 2),
          y: pad + Math.random() * Math.max(20, canvasSize.height - pad * 2),
          vx: (Math.random() - 0.5) * 0.8,
          vy: (Math.random() - 0.5) * 0.8,
        };
      }
      return next;
    });
  }, [filteredNodes, canvasSize]);

  useEffect(() => {
    if (!filteredNodes.length) return;
    const idSet = new Set(filteredNodes.map((n) => n.id));

    let raf = 0;
    const edgeStrength = 0.0011;
    const repulsion = compact ? 1100 : 1300;
    const damping = 0.9;
    const centerStrength = 0.00055;

    const step = () => {
      setPositions((prev) => {
        const next: Record<string, NodeState> = { ...prev };
        const ids = Object.keys(next).filter((id) => idSet.has(id));

        for (let i = 0; i < ids.length; i += 1) {
          const aId = ids[i];
          const a = next[aId];
          if (!a || draggingIdRef.current === aId) continue;

          for (let j = i + 1; j < ids.length; j += 1) {
            const bId = ids[j];
            const b = next[bId];
            if (!b || draggingIdRef.current === bId) continue;

            const dx = a.x - b.x;
            const dy = a.y - b.y;
            const d2 = dx * dx + dy * dy + 0.01;
            const force = repulsion / d2;

            a.vx += (dx / Math.sqrt(d2)) * force * 0.015;
            a.vy += (dy / Math.sqrt(d2)) * force * 0.015;
            b.vx -= (dx / Math.sqrt(d2)) * force * 0.015;
            b.vy -= (dy / Math.sqrt(d2)) * force * 0.015;
          }
        }

        for (const edge of edges) {
          const a = next[edge.source];
          const b = next[edge.target];
          if (!a || !b) continue;
          const dx = b.x - a.x;
          const dy = b.y - a.y;
          const dist = Math.sqrt(dx * dx + dy * dy) || 1;
          const target = edge.type === 'organism' ? 160 : 130;
          const pull = (dist - target) * edgeStrength;
          const fx = (dx / dist) * pull;
          const fy = (dy / dist) * pull;

          if (draggingIdRef.current !== edge.source) {
            a.vx += fx;
            a.vy += fy;
          }
          if (draggingIdRef.current !== edge.target) {
            b.vx -= fx;
            b.vy -= fy;
          }
        }

        const cx = canvasSize.width / 2;
        const cy = canvasSize.height / 2;
        const pad = 22;

        for (const id of ids) {
          const n = next[id];
          if (!n) continue;
          if (draggingIdRef.current === id) continue;

          n.vx += (cx - n.x) * centerStrength;
          n.vy += (cy - n.y) * centerStrength;

          n.vx *= damping;
          n.vy *= damping;
          n.x += n.vx;
          n.y += n.vy;

          n.x = clamp(n.x, pad, canvasSize.width - pad);
          n.y = clamp(n.y, pad, canvasSize.height - pad);
        }

        return next;
      });

      raf = window.requestAnimationFrame(step);
    };

    raf = window.requestAnimationFrame(step);
    return () => window.cancelAnimationFrame(raf);
  }, [filteredNodes, edges, canvasSize, compact]);

  const visibleNodes = useMemo(() => {
    return filteredNodes.map((node) => {
      const p = positions[node.id];
      return {
        ...node,
        x: p?.x ?? canvasSize.width * 0.5,
        y: p?.y ?? canvasSize.height * 0.5,
      };
    });
  }, [filteredNodes, positions, canvasSize]);

  const nodeById = useMemo(() => {
    const map = new Map<string, (NetworkNode & { x: number; y: number })>();
    for (const n of visibleNodes) map.set(n.id, n);
    return map;
  }, [visibleNodes]);

  useEffect(() => {
    if (!visibleNodes.length) {
      setSelectedNodeId('');
      return;
    }
    if (!selectedNodeId || !visibleNodes.some((n) => n.id === selectedNodeId)) {
      setSelectedNodeId(visibleNodes[0].id);
    }
  }, [visibleNodes, selectedNodeId]);

  const selected = selectedNodeId ? nodeById.get(selectedNodeId) || null : null;

  function pointerToSvg(clientX: number, clientY: number): { x: number; y: number } | null {
    const svg = svgRef.current;
    if (!svg) return null;
    const rect = svg.getBoundingClientRect();
    if (!rect.width || !rect.height) return null;
    const x = ((clientX - rect.left) / rect.width) * canvasSize.width;
    const y = ((clientY - rect.top) / rect.height) * canvasSize.height;
    return { x, y };
  }

  function onPointerMove(event: ReactPointerEvent<SVGSVGElement>): void {
    const dragging = draggingIdRef.current;
    if (!dragging) return;
    const p = pointerToSvg(event.clientX, event.clientY);
    if (!p) return;
    setPositions((prev) => ({
      ...prev,
      [dragging]: {
        ...(prev[dragging] || { x: p.x, y: p.y, vx: 0, vy: 0 }),
        x: clamp(p.x, 18, canvasSize.width - 18),
        y: clamp(p.y, 18, canvasSize.height - 18),
        vx: 0,
        vy: 0,
      },
    }));
  }

  function stopDragging(): void {
    draggingIdRef.current = null;
  }

  if (!papers.length) {
    return (
      <article className="paper-network-inline">
        <header className="paper-network-top">
          <h3>{title}</h3>
          <p>{subtitle}</p>
        </header>
        <p className="paper-network-empty">Run Pro analysis to generate the integrated paper network.</p>
      </article>
    );
  }

  return (
    <article className="paper-network-inline">
      <header className="paper-network-top">
        <div>
          <h3>{title}</h3>
          <p>{subtitle}</p>
        </div>
        <div className="paper-network-statline">
          <span>Papers: {filteredNodes.length}/{papers.length}</span>
          <span>Edges: {edges.length}</span>
        </div>
      </header>

      <section className="paper-network-controls">
        <label>
          <span>Min score</span>
          <input
            type="range"
            min={0}
            max={10}
            step={1}
            value={minScore}
            onChange={(e) => setMinScore(Number(e.target.value))}
          />
          <strong>&ge; {minScore}</strong>
        </label>

        <label>
          <span>Journal</span>
          <select value={selectedJournal} onChange={(e) => setSelectedJournal(e.target.value)}>
            <option value="all">All journals</option>
            {journals.map((journal) => (
              <option key={journal} value={journal}>{journal}</option>
            ))}
          </select>
        </label>

        <label>
          <span>Year range</span>
          <div className="paper-network-year-range">
            <input
              type="number"
              value={yearMin ?? ''}
              placeholder={years[0] ? String(years[0]) : 'min'}
              onChange={(e) => setYearMin(e.target.value ? Number(e.target.value) : null)}
            />
            <input
              type="number"
              value={yearMax ?? ''}
              placeholder={years[years.length - 1] ? String(years[years.length - 1]) : 'max'}
              onChange={(e) => setYearMax(e.target.value ? Number(e.target.value) : null)}
            />
          </div>
        </label>
      </section>

      {visibleNodes.length ? (
        <div className="paper-network-layout">
          <div className="paper-network-canvas-wrap">
            <svg
              ref={svgRef}
              width={canvasSize.width}
              height={canvasSize.height}
              viewBox={`0 0 ${canvasSize.width} ${canvasSize.height}`}
              role="img"
              aria-label="Paper relevance network"
              onPointerMove={onPointerMove}
              onPointerUp={stopDragging}
              onPointerLeave={stopDragging}
            >
              {edges.map((edge, idx) => {
                const source = nodeById.get(edge.source);
                const target = nodeById.get(edge.target);
                if (!source || !target) return null;
                return (
                  <line
                    key={`${edge.source}-${edge.target}-${idx}`}
                    x1={source.x}
                    y1={source.y}
                    x2={target.x}
                    y2={target.y}
                    className={`paper-network-edge ${edge.type}`}
                  />
                );
              })}

              {visibleNodes.map((node) => {
                const radius = 12 + Math.max(0, node.score) * (compact ? 1.1 : 1.25);
                const active = selectedNodeId === node.id;
                return (
                  <g
                    key={node.id}
                    onClick={() => setSelectedNodeId(node.id)}
                    onPointerDown={(event) => {
                      draggingIdRef.current = node.id;
                      setSelectedNodeId(node.id);
                      (event.currentTarget as SVGGElement).setPointerCapture(event.pointerId);
                    }}
                    className={`paper-network-node ${active ? 'active' : ''}`}
                  >
                    <circle
                      cx={node.x}
                      cy={node.y}
                      r={radius}
                      fill={scoreColor(node.score)}
                      opacity={node.isRelevant ? 0.92 : 0.68}
                    />
                    <text x={node.x} y={node.y + radius + 14} textAnchor="middle" className="paper-network-node-label">
                      {node.label}
                    </text>
                    {node.year ? (
                      <text x={node.x} y={node.y + radius + 28} textAnchor="middle" className="paper-network-node-year">
                        {node.year}
                      </text>
                    ) : null}
                  </g>
                );
              })}
            </svg>
          </div>

          <aside className="paper-network-details">
            {selected ? (
              <>
                <h4>{selected.title}</h4>
                <dl>
                  <div>
                    <dt>Author label</dt>
                    <dd>{selected.label}</dd>
                  </div>
                  <div>
                    <dt>PMID</dt>
                    <dd>{selected.pmid}</dd>
                  </div>
                  <div>
                    <dt>Score</dt>
                    <dd>{selected.score}/10</dd>
                  </div>
                  <div>
                    <dt>Year</dt>
                    <dd>{selected.yearRaw}</dd>
                  </div>
                  <div>
                    <dt>Journal</dt>
                    <dd>{selected.journal}</dd>
                  </div>
                  <div>
                    <dt>Relevant</dt>
                    <dd>{selected.isRelevant ? 'Yes' : 'No'}</dd>
                  </div>
                </dl>
              </>
            ) : (
              <p>Select a node to inspect metadata.</p>
            )}
          </aside>
        </div>
      ) : (
        <p className="paper-network-empty">No papers match current filters.</p>
      )}

      <footer className="paper-network-legend">
        <span><i className="dot high" /> High (7-10)</span>
        <span><i className="dot mid" /> Mid (4-6)</span>
        <span><i className="dot low" /> Low (0-3)</span>
        <span><i className="line journal" /> Shared journal</span>
        <span><i className="line organism" /> Shared organism</span>
      </footer>
    </article>
  );
}
