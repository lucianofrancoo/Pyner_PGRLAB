"""
Paper Visualizer — Interactive HTML Network Graph v2
====================================================
Generates a self-contained HTML file with a D3.js force-directed graph.
Features:
  - Hierarchical layout by relevance score (if classified data available)
  - Ordered by year (recent on top)
  - Side panel with filters (year, journal, relevance)
  - Auto-detection of classified TSV
  - Dark/Light mode toggle
  - Node labels with first 5 words of title
"""

import json
import csv
import os
import re
import html
import logging
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict

logger = logging.getLogger(__name__)


class PaperVisualizer:
    """Generates interactive HTML visualizations from paper data"""
    
    def __init__(self):
        self.papers = []
        self.classified_data = {}
        self.query = ""
    
    def load_data(self, json_path: str) -> int:
        """Load paper data from fetcher JSON output."""
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        self.query = data.get('metadata', {}).get('query', '')
        self.papers = data.get('publications', [])
        logger.info(f"Loaded {len(self.papers)} papers from {json_path}")
        return len(self.papers)
    
    def auto_detect_classified(self, json_path: str) -> Optional[str]:
        """
        Auto-detect associated classified TSV from JSON path.
        pubmed_fetch_TIMESTAMP.json → classified_papers_TIMESTAMP.tsv
        """
        basename = os.path.basename(json_path)
        match = re.search(r'pubmed_fetch_(\d+_\d+)', basename)
        if match:
            timestamp = match.group(1)
            classified_path = os.path.join(
                os.path.dirname(json_path),
                f'classified_papers_{timestamp}.tsv'
            )
            if os.path.exists(classified_path):
                logger.info(f"Auto-detected classified TSV: {classified_path}")
                return classified_path
        return None
    
    def load_classified_data(self, tsv_path: str) -> int:
        """Load classified data from analyzer TSV output."""
        with open(tsv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                pmid = row.get('PMID', '')
                if pmid:
                    self.classified_data[pmid] = row
        
        logger.info(f"Loaded {len(self.classified_data)} classified papers from {tsv_path}")
        return len(self.classified_data)
    
    def _short_title(self, title: str, n_words: int = 5) -> str:
        """Get first N words of title + ellipsis"""
        words = (title or 'No title').split()
        if len(words) <= n_words:
            return title
        return ' '.join(words[:n_words]) + '…'
    
    def _build_nodes(self) -> List[Dict]:
        """Build node data for each paper"""
        nodes = []
        for i, paper in enumerate(self.papers):
            pmid = paper.get('pmid', f'unknown_{i}')
            classified = self.classified_data.get(pmid, {})
            
            try:
                score = int(classified.get('Relevance_Score', 5))
            except (ValueError, TypeError):
                score = 5
            
            is_relevant = classified.get('Is_Relevant', 'Unknown')
            summary = classified.get('Summary', '')
            explanation = classified.get('Relevance_Explanation', '')
            organisms = classified.get('Organisms', '')
            tissues = classified.get('Tissues_Organs', '')
            conditions = classified.get('Conditions', '')
            strategies = classified.get('Strategies', '')
            
            authors = paper.get('authors', [])
            if len(authors) > 3:
                author_str = ', '.join(authors[:3]) + f' et al. (+{len(authors)-3})'
            else:
                author_str = ', '.join(authors) if authors else 'Unknown'
            
            title = paper.get('title', 'No title') or 'No title'
            year = paper.get('year', '2020')
            journal = paper.get('journal', 'Unknown')
            
            node = {
                'id': pmid,
                'label': self._short_title(title),
                'full_title': title,
                'authors': author_str,
                'year': year,
                'journal': journal,
                'doi': paper.get('doi', ''),
                'pmcid': paper.get('pmcid', ''),
                'url': paper.get('url', f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'),
                'abstract': (paper.get('abstract', '') or '')[:400],
                'score': score,
                'is_relevant': is_relevant,
                'summary': summary,
                'explanation': explanation,
                'organisms': organisms,
                'tissues': tissues,
                'conditions': conditions,
                'strategies': strategies,
                'group': self._get_group(score),
                'radius': max(14, 10 + score * 3),
            }
            nodes.append(node)
        
        return nodes
    
    def _get_group(self, score: int) -> int:
        if score >= 7:
            return 0
        elif score >= 4:
            return 1
        else:
            return 2
    
    def _build_edges(self, nodes: List[Dict]) -> List[Dict]:
        """Build edges between papers sharing organisms or journals"""
        edges = []
        seen = set()
        
        organism_papers = defaultdict(list)
        for node in nodes:
            orgs = node.get('organisms', '')
            if orgs and orgs != 'not described':
                for org in orgs.split(' ; '):
                    org = org.strip().lower()
                    if org and org != 'not described':
                        organism_papers[org].append(node['id'])
        
        for org, pmids in organism_papers.items():
            if len(pmids) > 1:
                for i in range(len(pmids)):
                    for j in range(i+1, len(pmids)):
                        key = tuple(sorted([pmids[i], pmids[j]]))
                        if key not in seen:
                            seen.add(key)
                            edges.append({'source': pmids[i], 'target': pmids[j], 'type': 'organism', 'label': org})
        
        journal_papers = defaultdict(list)
        for node in nodes:
            j = node.get('journal', '').strip().lower()
            if j:
                journal_papers[j].append(node['id'])
        
        for journal, pmids in journal_papers.items():
            if 2 <= len(pmids) <= 8:
                for i in range(len(pmids)):
                    for j in range(i+1, len(pmids)):
                        key = tuple(sorted([pmids[i], pmids[j]]))
                        if key not in seen:
                            seen.add(key)
                            edges.append({'source': pmids[i], 'target': pmids[j], 'type': 'journal', 'label': journal})
        
        return edges
    
    def generate_html(self, output_path: str) -> str:
        """Generate self-contained interactive HTML visualization."""
        nodes = self._build_nodes()
        edges = self._build_edges(nodes)
        
        has_scores = bool(self.classified_data)
        total = len(nodes)
        relevant = sum(1 for n in nodes if n['is_relevant'] == 'Yes')
        
        # Collect unique years & journals for filters
        years = sorted(set(n['year'] for n in nodes if n['year']), reverse=True)
        journals = sorted(set(n['journal'] for n in nodes if n['journal']))
        
        html_content = self._render_html(nodes, edges, has_scores, total, relevant, years, journals)
        
        os.makedirs(os.path.dirname(output_path) if os.path.dirname(output_path) else '.', exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        logger.info(f"✓ Generated visualization: {output_path}")
        return output_path
    
    def _render_html(self, nodes, edges, has_scores, total, relevant, years, journals) -> str:
        nodes_json = json.dumps(nodes, ensure_ascii=False)
        edges_json = json.dumps(edges, ensure_ascii=False)
        years_json = json.dumps(years)
        journals_json = json.dumps(journals, ensure_ascii=False)
        query_escaped = html.escape(self.query[:200])
        
        return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>PYNER — Paper Network</title>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<style>
:root {{
  --bg: #f8fafc;
  --bg2: #ffffff;
  --bg3: #f1f5f9;
  --text: #1e293b;
  --text2: #64748b;
  --text3: #94a3b8;
  --border: rgba(0,0,0,0.08);
  --accent: #3b82f6;
  --green: #059669;
  --yellow: #d97706;
  --red: #dc2626;
  --panel-bg: rgba(255,255,255,0.97);
  --card-bg: rgba(241,245,249,0.9);
}}

[data-theme="dark"] {{
  --bg: #0a0e17;
  --bg2: #111827;
  --bg3: #1e293b;
  --text: #e0e6f0;
  --text2: #8892a8;
  --text3: #4a5568;
  --border: rgba(255,255,255,0.06);
  --panel-bg: rgba(15,20,35,0.96);
  --card-bg: rgba(30,40,60,0.8);
}}

* {{ margin: 0; padding: 0; box-sizing: border-box; }}

body {{
    font-family: 'Inter', system-ui, sans-serif;
    background: var(--bg);
    color: var(--text);
    overflow: hidden;
    height: 100vh;
    transition: background 0.3s, color 0.3s;
}}

/* Header */
.header {{
    position: fixed;
    top: 0;
    left: 0;
    right: 0;
    z-index: 100;
    padding: 10px 24px;
    background: var(--panel-bg);
    backdrop-filter: blur(12px);
    border-bottom: 1px solid var(--border);
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 16px;
}}

.header h1 {{
    font-size: 16px;
    font-weight: 700;
    background: linear-gradient(135deg, #6ee7b7, #3b82f6);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    white-space: nowrap;
}}

.header-stats {{
    display: flex;
    gap: 16px;
    font-size: 12px;
    color: var(--text2);
}}

.header-stats .val {{ color: var(--text); font-weight: 600; }}

.header-right {{
    display: flex;
    align-items: center;
    gap: 8px;
}}

.legend {{
    display: flex;
    gap: 12px;
    font-size: 11px;
    color: var(--text2);
}}

.legend-item {{ display: flex; align-items: center; gap: 4px; }}
.legend-dot {{ width: 8px; height: 8px; border-radius: 50%; }}

/* Theme toggle */
.theme-btn {{
    padding: 6px 12px;
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--text2);
    font-size: 12px;
    cursor: pointer;
    transition: all 0.2s;
}}
.theme-btn:hover {{ color: var(--text); }}

/* Side panel */
.side-panel {{
    position: fixed;
    top: 48px;
    left: 0;
    bottom: 0;
    width: 240px;
    background: var(--panel-bg);
    backdrop-filter: blur(12px);
    border-right: 1px solid var(--border);
    z-index: 90;
    padding: 16px;
    overflow-y: auto;
    transition: transform 0.3s;
}}

.side-panel.collapsed {{
    transform: translateX(-240px);
}}

.filter-section {{
    margin-bottom: 16px;
}}

.filter-title {{
    font-size: 10px;
    font-weight: 700;
    color: var(--green);
    text-transform: uppercase;
    letter-spacing: 1px;
    margin-bottom: 8px;
}}

.filter-option {{
    display: flex;
    align-items: center;
    gap: 6px;
    padding: 4px 0;
    font-size: 12px;
    color: var(--text2);
    cursor: pointer;
    transition: color 0.15s;
}}
.filter-option:hover {{ color: var(--text); }}
.filter-option input {{ accent-color: var(--accent); }}
.filter-option.active {{ color: var(--text); }}

.filter-count {{
    margin-left: auto;
    font-size: 10px;
    color: var(--text3);
    background: var(--bg3);
    padding: 1px 6px;
    border-radius: 8px;
}}

/* Relevance slider */
.relevance-slider {{
    width: 100%;
    margin: 8px 0;
    accent-color: var(--accent);
}}

.relevance-display {{
    font-size: 12px;
    color: var(--text2);
    text-align: center;
}}

/* SVG */
#graph {{
    width: 100%;
    height: 100vh;
    display: block;
}}

/* Node labels */
.node-label {{
    font-size: 9px;
    fill: var(--text2);
    pointer-events: none;
    text-anchor: middle;
    font-weight: 500;
}}

.node-year {{
    font-size: 8px;
    fill: var(--text3);
    pointer-events: none;
    text-anchor: middle;
}}

/* Tooltip */
.tooltip {{
    position: fixed;
    padding: 16px 20px;
    background: var(--panel-bg);
    backdrop-filter: blur(16px);
    border: 1px solid var(--border);
    border-radius: 12px;
    font-size: 13px;
    line-height: 1.5;
    max-width: 440px;
    pointer-events: none;
    opacity: 0;
    transition: opacity 0.15s;
    z-index: 200;
    box-shadow: 0 8px 32px rgba(0,0,0,0.3);
}}
.tooltip.visible {{ opacity: 1; }}
.tooltip .tt-title {{ font-weight: 700; font-size: 14px; color: var(--text); margin-bottom: 6px; line-height: 1.4; }}
.tooltip .tt-meta {{ color: var(--text2); font-size: 11px; margin-bottom: 6px; }}
.tooltip .tt-score {{ display: inline-block; padding: 2px 10px; border-radius: 12px; font-weight: 700; font-size: 12px; margin-bottom: 6px; }}
.tooltip .tt-score.high {{ background: rgba(16,185,129,0.2); color: #34d399; }}
.tooltip .tt-score.mid {{ background: rgba(245,158,11,0.2); color: #fbbf24; }}
.tooltip .tt-score.low {{ background: rgba(239,68,68,0.2); color: #f87171; }}
.tooltip .tt-section {{ margin-top: 6px; padding-top: 6px; border-top: 1px solid var(--border); }}
.tooltip .tt-label {{ font-size: 10px; font-weight: 600; color: var(--green); text-transform: uppercase; letter-spacing: 0.5px; margin-bottom: 2px; }}
.tooltip .tt-value {{ color: var(--text2); font-size: 12px; }}

/* Edges */
.edge {{ stroke-opacity: 0.12; stroke-width: 1; }}
.edge.organism {{ stroke: var(--accent); }}
.edge.journal {{ stroke: var(--text3); stroke-dasharray: 4 4; }}

/* Controls */
.controls {{
    position: fixed;
    bottom: 16px;
    right: 16px;
    display: flex;
    gap: 6px;
    z-index: 100;
}}
.controls button {{
    padding: 6px 14px;
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--text2);
    font-size: 11px;
    cursor: pointer;
    backdrop-filter: blur(8px);
    transition: all 0.2s;
}}
.controls button:hover {{ color: var(--text); }}

.info-panel {{
    position: fixed;
    bottom: 16px;
    left: 250px;
    font-size: 10px;
    color: var(--text3);
    z-index: 100;
    transition: left 0.3s;
}}
.info-panel.shifted {{ left: 16px; }}
</style>
</head>
<body>

<div class="header">
    <h1>🔬 PYNER — Paper Network</h1>
    <div class="header-stats">
        <span>Papers: <span class="val">{total}</span></span>
        {"<span>Relevant: <span class='val'>" + str(relevant) + " / " + str(total) + "</span></span>" if has_scores else ""}
    </div>
    <div class="header-right">
        <div class="legend">
            {"".join([
                '<div class="legend-item"><div class="legend-dot" style="background:#34d399;"></div>≥7</div>',
                '<div class="legend-item"><div class="legend-dot" style="background:#fbbf24;"></div>4-6</div>',
                '<div class="legend-item"><div class="legend-dot" style="background:#f87171;"></div>≤3</div>',
            ]) if has_scores else '<div class="legend-item"><div class="legend-dot" style="background:#3b82f6;"></div>Papers</div>'}
        </div>
        <button class="theme-btn" onclick="toggleTheme()" title="Toggle dark/light mode">🌓</button>
    </div>
</div>

<div class="side-panel" id="sidePanel">
    {"" if not has_scores else '''
    <div class="filter-section">
        <div class="filter-title">Relevance Score</div>
        <input type="range" class="relevance-slider" id="relevanceSlider" min="0" max="10" value="0" oninput="applyFilters()">
        <div class="relevance-display" id="relevanceDisplay">≥ 0</div>
    </div>
    '''}
    <div class="filter-section" id="yearFilters">
        <div class="filter-title"><label><input type="checkbox" checked id="yearSelectAll" onchange="toggleAllYears(this)"> Year</label></div>
    </div>
    <div class="filter-section" id="journalFilters">
        <div class="filter-title"><label><input type="checkbox" checked id="journalSelectAll" onchange="toggleAllJournals(this)"> Journal</label></div>
    </div>
</div>

<div class="tooltip" id="tooltip"></div>
<svg id="graph"></svg>

<div class="controls">
    <button onclick="togglePanel()">☰ Filters</button>
    <button onclick="resetZoom()">Reset View</button>
    <button onclick="toggleLabels()">Labels</button>
</div>

<div class="info-panel" id="infoPanel">
    Drag nodes · Scroll to zoom · Click to open paper
</div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
(function() {{
    const nodesData = {nodes_json};
    const edgesData = {edges_json};
    const hasScores = {str(has_scores).lower()};
    const allYears = {years_json};
    const allJournals = {journals_json};
    
    // State
    let activeYears = new Set(allYears);
    let activeJournals = new Set(allJournals);
    let minRelevance = 0;
    let showLabels = true;
    let panelOpen = true;
    
    // Colors
    const colorScale = d3.scaleOrdinal().domain([0, 1, 2]).range(['#34d399', '#fbbf24', '#f87171']);
    const defaultColor = '#3b82f6';
    function nodeColor(d) {{ return hasScores ? colorScale(d.group) : defaultColor; }}
    
    // Build filter UI
    const yearDiv = document.getElementById('yearFilters');
    allYears.forEach(y => {{
        const count = nodesData.filter(n => n.year === y).length;
        yearDiv.innerHTML += `<label class="filter-option active"><input type="checkbox" checked value="${{y}}" class="year-cb" onchange="applyFilters()"> ${{y}} <span class="filter-count">${{count}}</span></label>`;
    }});
    
    const journalDiv = document.getElementById('journalFilters');
    allJournals.forEach(j => {{
        const count = nodesData.filter(n => n.journal === j).length;
        const short = j.length > 28 ? j.substring(0, 26) + '…' : j;
        journalDiv.innerHTML += `<label class="filter-option active"><input type="checkbox" checked value="${{j}}" class="journal-cb" onchange="applyFilters()"> ${{short}} <span class="filter-count">${{count}}</span></label>`;
    }});
    
    // Select All toggles
    window.toggleAllYears = (master) => {{
        document.querySelectorAll('.year-cb').forEach(cb => {{ cb.checked = master.checked; }});
        applyFilters();
    }};
    window.toggleAllJournals = (master) => {{
        document.querySelectorAll('.journal-cb').forEach(cb => {{ cb.checked = master.checked; }});
        applyFilters();
    }};
    
    // Dimensions
    const panelW = 240;
    let width = window.innerWidth - panelW;
    let height = window.innerHeight;
    const offsetX = panelW;
    
    const svg = d3.select('#graph').attr('width', window.innerWidth).attr('height', height);
    const g = svg.append('g');
    
    const zoom = d3.zoom().scaleExtent([0.15, 5]).on('zoom', e => g.attr('transform', e.transform));
    svg.call(zoom);
    
    // Layout bounds — keep everything inside a reasonable area
    const bndLeft = offsetX + 60;
    const bndRight = offsetX + width - 60;
    const bndTop = 80;
    const bndBottom = height - 60;
    const cx = (bndLeft + bndRight) / 2;
    const cy = (bndTop + bndBottom) / 2;
    
    // Year scale for vertical positioning (recent → top)
    const yearExtent = d3.extent(nodesData, d => +d.year || 2020);
    const yScale = d3.scaleLinear()
        .domain([yearExtent[1], yearExtent[0]])  // Recent on top
        .range([bndTop + 40, bndBottom - 40]);
    
    // Scale forces based on number of nodes — gentle repulsion, strong centering
    const n = nodesData.length;
    const chargeStr = n > 40 ? -80 : n > 20 ? -120 : -180;
    const linkDist = n > 40 ? 80 : n > 20 ? 100 : 120;
    
    // Simulation
    const simulation = d3.forceSimulation(nodesData)
        .force('link', d3.forceLink(edgesData).id(d => d.id).distance(linkDist).strength(0.2))
        .force('charge', d3.forceManyBody().strength(chargeStr))
        .force('x', d3.forceX(cx).strength(0.06))
        .force('y', d3.forceY(d => {{
            if (hasScores) {{
                const scoreY = (10 - d.score) / 10;
                const yearY = yScale(+d.year || 2020);
                return scoreY * 150 + yearY * 0.7;
            }}
            return yScale(+d.year || 2020);
        }}).strength(0.07))
        .force('collision', d3.forceCollide().radius(d => d.radius + 6))
        .alphaDecay(0.025);
    
    // Draw edges
    const link = g.append('g').selectAll('line').data(edgesData).join('line')
        .attr('class', d => 'edge ' + d.type);
    
    // Draw nodes
    const node = g.append('g').selectAll('g').data(nodesData).join('g')
        .attr('cursor', 'pointer')
        .call(d3.drag().on('start', dragstarted).on('drag', dragged).on('end', dragended));
    
    node.append('circle')
        .attr('r', d => d.radius)
        .attr('fill', nodeColor)
        .attr('fill-opacity', 0.7)
        .attr('stroke', nodeColor)
        .attr('stroke-opacity', 0.3)
        .attr('stroke-width', 2);
    
    node.append('circle')
        .attr('r', d => d.radius + 4)
        .attr('fill', 'none')
        .attr('stroke', nodeColor)
        .attr('stroke-opacity', 0.08)
        .attr('stroke-width', 6);
    
    // Labels: first 5 words of title
    const labels = node.append('text')
        .attr('class', 'node-label')
        .attr('dy', d => d.radius + 13)
        .text(d => d.label);
    
    // Year sub-label
    node.append('text')
        .attr('class', 'node-year')
        .attr('dy', d => d.radius + 23)
        .text(d => d.year);
    
    // Tooltip
    const tooltip = document.getElementById('tooltip');
    
    node.on('mouseover', (event, d) => {{
        const connected = new Set([d.id]);
        edgesData.forEach(e => {{
            const s = typeof e.source === 'object' ? e.source.id : e.source;
            const t = typeof e.target === 'object' ? e.target.id : e.target;
            if (s === d.id) connected.add(t);
            if (t === d.id) connected.add(s);
        }});
        
        node.select('circle:first-child').attr('fill-opacity', n => connected.has(n.id) ? 0.95 : 0.1);
        link.attr('stroke-opacity', e => {{
            const s = typeof e.source === 'object' ? e.source.id : e.source;
            const t = typeof e.target === 'object' ? e.target.id : e.target;
            return (s === d.id || t === d.id) ? 0.5 : 0.02;
        }});
        
        let sc = d.score >= 7 ? 'high' : d.score >= 4 ? 'mid' : 'low';
        let h = `<div class="tt-title">${{d.full_title}}</div>`;
        h += `<div class="tt-meta">${{d.authors}} · ${{d.journal}} · ${{d.year}}</div>`;
        if (hasScores) h += `<span class="tt-score ${{sc}}">Relevance: ${{d.score}}/10</span>`;
        if (d.explanation) h += `<div class="tt-section"><div class="tt-label">Why relevant</div><div class="tt-value">${{d.explanation}}</div></div>`;
        if (d.summary) h += `<div class="tt-section"><div class="tt-label">Summary</div><div class="tt-value">${{d.summary.substring(0,200)}}${{d.summary.length>200?'…':''}}</div></div>`;
        if (d.organisms && d.organisms !== 'not described') h += `<div class="tt-section"><div class="tt-label">Organisms</div><div class="tt-value">${{d.organisms}}</div></div>`;
        if (d.conditions && d.conditions !== 'not described') h += `<div class="tt-section"><div class="tt-label">Conditions</div><div class="tt-value">${{d.conditions}}</div></div>`;
        if (d.strategies && d.strategies !== 'not described') h += `<div class="tt-section"><div class="tt-label">Strategies</div><div class="tt-value">${{d.strategies}}</div></div>`;
        h += `<div class="tt-section"><div class="tt-meta">PMID: ${{d.id}}${{d.doi ? ' · DOI: '+d.doi : ''}}</div></div>`;
        
        tooltip.innerHTML = h;
        tooltip.classList.add('visible');
    }})
    .on('mousemove', event => {{
        const x = event.clientX + 16;
        const y = event.clientY - 16;
        const r = tooltip.getBoundingClientRect();
        tooltip.style.left = (x + r.width > window.innerWidth ? event.clientX - r.width - 16 : x) + 'px';
        tooltip.style.top = (y + r.height > window.innerHeight ? window.innerHeight - r.height - 8 : y) + 'px';
    }})
    .on('mouseout', () => {{
        tooltip.classList.remove('visible');
        node.select('circle:first-child').attr('fill-opacity', 0.7);
        link.attr('stroke-opacity', 0.12);
    }})
    .on('click', (e, d) => {{ if (d.url) window.open(d.url, '_blank'); }});
    
    // Tick — clamp positions within bounds so no node escapes
    simulation.on('tick', () => {{
        nodesData.forEach(d => {{
            if (!d.fx) d.x = Math.max(bndLeft, Math.min(bndRight, d.x));
            if (!d.fy) d.y = Math.max(bndTop, Math.min(bndBottom, d.y));
        }});
        link.attr('x1', d => d.source.x).attr('y1', d => d.source.y)
            .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
        node.attr('transform', d => `translate(${{d.x}},${{d.y}})`);
    }});
    
    // Auto-fit after simulation stabilizes
    simulation.on('end', () => {{ fitToView(); }});
    
    // Drag — nodes stay pinned where you drop them; double-click to release
    function dragstarted(e, d) {{ if (!e.active) simulation.alphaTarget(0.1).restart(); d.fx = d.x; d.fy = d.y; }}
    function dragged(e, d) {{ d.fx = e.x; d.fy = e.y; }}
    function dragended(e, d) {{ if (!e.active) simulation.alphaTarget(0); /* keep pinned: fx/fy stay */ }}
    
    // Double-click to unpin a node
    node.on('dblclick', (e, d) => {{ d.fx = null; d.fy = null; e.stopPropagation(); }});
    
    // Filters
    window.applyFilters = () => {{
        activeYears.clear();
        document.querySelectorAll('#yearFilters input:checked').forEach(cb => activeYears.add(cb.value));
        activeJournals.clear();
        document.querySelectorAll('#journalFilters input:checked').forEach(cb => activeJournals.add(cb.value));
        
        const slider = document.getElementById('relevanceSlider');
        if (slider) {{
            minRelevance = +slider.value;
            document.getElementById('relevanceDisplay').textContent = '≥ ' + minRelevance;
        }}
        
        node.attr('visibility', d => {{
            const yearOk = activeYears.has(d.year);
            const journalOk = activeJournals.has(d.journal);
            const scoreOk = !hasScores || d.score >= minRelevance;
            return (yearOk && journalOk && scoreOk) ? 'visible' : 'hidden';
        }});
        
        const visibleIds = new Set();
        node.each(function(d) {{ if (d3.select(this).attr('visibility') !== 'hidden') visibleIds.add(d.id); }});
        
        link.attr('visibility', d => {{
            const s = typeof d.source === 'object' ? d.source.id : d.source;
            const t = typeof d.target === 'object' ? d.target.id : d.target;
            return visibleIds.has(s) && visibleIds.has(t) ? 'visible' : 'hidden';
        }});
    }};
    
    // Controls
    function fitToView() {{
        const pad = 60;
        let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
        nodesData.forEach(d => {{
            if (d.x - d.radius < minX) minX = d.x - d.radius;
            if (d.x + d.radius > maxX) maxX = d.x + d.radius;
            if (d.y - d.radius < minY) minY = d.y - d.radius;
            if (d.y + d.radius > maxY) maxY = d.y + d.radius;
        }});
        const bw = maxX - minX + pad * 2;
        const bh = maxY - minY + pad * 2;
        const vw = window.innerWidth;
        const vh = window.innerHeight;
        const scale = Math.min(vw / bw, vh / bh, 1.5);
        const tx = vw / 2 - (minX + maxX) / 2 * scale;
        const ty = vh / 2 - (minY + maxY) / 2 * scale;
        svg.transition().duration(600).call(zoom.transform,
            d3.zoomIdentity.translate(tx, ty).scale(scale)
        );
    }}
    window.resetZoom = fitToView;
    window.toggleLabels = () => {{
        showLabels = !showLabels;
        const v = showLabels ? 'visible' : 'hidden';
        node.selectAll('.node-label, .node-year').attr('visibility', v);
    }};
    window.togglePanel = () => {{
        panelOpen = !panelOpen;
        document.getElementById('sidePanel').classList.toggle('collapsed');
        document.getElementById('infoPanel').classList.toggle('shifted');
    }};
    window.toggleTheme = () => {{
        const body = document.body;
        const current = body.getAttribute('data-theme');
        body.setAttribute('data-theme', current === 'dark' ? '' : 'dark');
    }};
    
    // Responsive
    window.addEventListener('resize', () => {{
        width = window.innerWidth - (panelOpen ? panelW : 0);
        height = window.innerHeight;
        svg.attr('width', window.innerWidth).attr('height', height);
        simulation.force('x', d3.forceX(d => (panelOpen ? panelW : 0) + width / 2).strength(0.05));
        simulation.alpha(0.3).restart();
    }});
}})();
</script>
</body>
</html>'''

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Generate interactive HTML paper network')
    parser.add_argument('json_path', help='Path to pubmed_fetch_*.json')
    parser.add_argument('--classified', '-c', help='Path to classified_papers_*.tsv (auto-detected if omitted)', default=None)
    parser.add_argument('--output', '-o', help='Output HTML path', default=None)
    
    args = parser.parse_args()
    
    if not os.path.exists(args.json_path):
        print(f"Error: File not found: {args.json_path}")
        sys.exit(1)
    
    # Output path
    if args.output:
        output_path = args.output
    else:
        basename = os.path.basename(args.json_path)
        match = re.search(r'pubmed_fetch_(\d+_\d+)', basename)
        timestamp = match.group(1) if match else 'output'
        output_path = os.path.join(os.path.dirname(args.json_path) or 'output', f'paper_network_{timestamp}.html')
    
    # Generate
    print("\n" + "=" * 60)
    print("🔬 PYNER — Paper Network Visualization")
    print("=" * 60)
    
    viz = PaperVisualizer()
    n_papers = viz.load_data(args.json_path)
    print(f"  📄 Loaded {n_papers} papers from JSON")
    
    # Auto-detect or use explicit classified path
    classified_path = args.classified
    if not classified_path:
        classified_path = viz.auto_detect_classified(args.json_path)
    
    if classified_path and os.path.exists(classified_path):
        n_classified = viz.load_classified_data(classified_path)
        print(f"  📊 Auto-loaded {n_classified} classified results")
        print(f"  ✓ Mode: Enriched (hierarchical by relevance)")
    else:
        print(f"  ✓ Mode: Basic (no classified data found)")
    
    final_path = viz.generate_html(output_path)
    print(f"\n  ✓ Generated: {final_path}")
    print(f"  📦 Size: {os.path.getsize(final_path) / 1024:.1f} KB")
    print("=" * 60)
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main())

