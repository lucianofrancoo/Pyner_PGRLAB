#!/usr/bin/env python3
"""
Generate a boolean NCBI query from semantic search results.

Usage:
  python3 Query_generator/generate_boolean_query.py "Mouse RNAseq"

The script calls the local API `/search`, collects the top results and
maps them to synonym groups to produce a boolean query ready for NCBI/Entrez.
"""

import sys
import json
import re
import subprocess
from typing import List


COMMON_MAP = {
    'mouse': 'Mus musculus',
    'mice': 'Mus musculus',
    'rat': 'Rattus norvegicus',
    'rattus': 'Rattus norvegicus',
    'human': 'Homo sapiens',
    'homo sapiens': 'Homo sapiens',
    'arabidopsis': 'Arabidopsis thaliana',
}

STRATEGY_SYNS = [ '"RNA-Seq"', '"RNA sequencing"', 'transcriptome', 'rnaseq' ]
GENE_EXPR_SYNS = [ '"gene expression"', 'transcriptome', '"RNA-Seq"', 'expression' ]
CHIP_SYNS = [ '"ChIP-Seq"', '"ChIP sequencing"', '"chromatin immunoprecipitation"' ]
GUT_SYNS = [ 'gut', 'intestinal', 'fecal' ]
METAGENOME_SYNS = [ 'metagenome', 'metatranscriptome' ]

# Keywords that suggest a term is NOT an organism
NON_ORGANISM_KEYWORDS = {
    'drought', 'stress', 'disease', 'infection', 'treatment', 'condition',
    'root', 'leaf', 'flower', 'tissue', 'organ', 'development', 'growth',
    'sequencing', 'rnaseq', 'chip-seq', 'strategy', 'method', 'analysis',
}


def call_search_api(query: str, top_k: int = 5, expand: bool = False) -> dict:
    cmd = [
        'curl', '-s', '-X', 'POST',
        'http://localhost:8000/search',
        '-H', 'Content-Type: application/json',
        '-d', json.dumps({'query': query, 'top_k': top_k, 'expand': expand})
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
    if proc.returncode != 0:
        raise RuntimeError(f"API call failed: {proc.stderr}")
    return json.loads(proc.stdout)

def extract_species_from_text(text: str) -> List[str]:
    text_l = text.lower()
    found = []
    # look for latin name pattern
    m = re.search(r'([A-Z][a-z]+)\s+([a-z]+)', text)
    if m:
        latin = f"{m.group(1)} {m.group(2)}"
        found.append(latin)
    # look for common names (whole-word match)
    for key, latin in COMMON_MAP.items():
        if re.search(rf"\b{re.escape(key)}\b", text_l):
            if latin not in found:
                found.append(latin)
    return found


def load_kb_organisms():
    """Load organisms from the copied KB to validate organism detection."""
    try:
        from pathlib import Path
        kb_path = Path(__file__).resolve().parents[0] / 'phases' / 'phase1' / 'output' / 'stage3_knowledge_base.json'
        if not kb_path.exists():
            return set()
        with open(kb_path, 'r') as f:
            kb = json.load(f)
        orgs = set(k.lower() for k in kb.get('organisms', {}).keys())
        return orgs
    except Exception:
        return set()


def build_boolean(results: List[dict], seed_organisms: List[str] = None) -> str:
    organisms = []
    facets = []
    kb_orgs = load_kb_organisms()

    for r in results:
        qtext = r.get('query_text', '')
        qtype = r.get('query_type', '').lower()

        # Always add the raw query phrase as fallback (quoted)
        raw_phrase = f'"{qtext}"'

        if qtype == 'organism':
            # Use KB organisms to decide whether this text really names an organism
            species = extract_species_from_text(qtext)
            matched = []
            for s in species:
                if s.lower() in kb_orgs:
                    matched.append(s)
            # Also check for common name occurrences from KB
            for org in kb_orgs:
                if org in qtext.lower() and org not in matched:
                    matched.append(org)

            if matched:
                for s in matched:
                    label = s
                    # if KB provided latin name, keep it; else use as-is
                    organisms.append(f"{label}[Organism]")
            else:
                # Not a validated organism — treat as All Fields
                facets.append([raw_phrase])

        elif qtype == 'strategy':
            # Strategy: map known sequencing strategies
            if 'chip' in qtext.lower():
                facets.append(CHIP_SYNS)
            elif 'rna' in qtext.lower() or 'rnaseq' in qtext.lower() or 'rna-seq' in qtext.lower():
                facets.append(STRATEGY_SYNS)
            else:
                facets.append([raw_phrase])

        elif qtype == 'gene_expression':
            facets.append(GENE_EXPR_SYNS)

        elif 'metagenome' in qtext.lower() or 'gut' in qtext.lower():
            facets.append(METAGENOME_SYNS + GUT_SYNS)

        else:
            # Generic fallback: include quoted phrase
            facets.append([raw_phrase])

    # Build clauses
    clause_list = []
    # include seed organisms detected from the original query
    if seed_organisms:
        for s in seed_organisms:
            if s not in organisms:
                organisms.append(s + '[Organism]')

    if organisms:
        clause_list.append('(' + ' OR '.join(organisms) + ')')

    # merge facets into a single AND chain; within each facet use OR
    # deduplicate facet groups (by tuple of terms)
    uniq_facets = []
    seen = set()
    for facet in facets:
        key = tuple(facet)
        if key in seen:
            continue
        seen.add(key)
        uniq_facets.append(facet)

    for facet in uniq_facets:
        uniq = []
        for t in facet:
            if t not in uniq:
                uniq.append(t)
        clause_list.append('(' + ' OR '.join(uniq) + ')')

    if not clause_list:
        return ''
    return ' AND '.join(clause_list)


def local_retrieve(query: str, top_k: int = 5):
    """Fallback: use local phase2/phase3 copies inside Query_generator/phases."""
    from pathlib import Path
    phases_path = Path(__file__).resolve().parents[0] / 'phases'
    import sys
    if str(phases_path) not in sys.path:
        sys.path.insert(0, str(phases_path))

    try:
        from phase2.scripts.vector_db import VectorDatabase, Retriever
    except Exception as e:
        raise RuntimeError(f"Local vector DB import failed: {e}")

    vector_db = VectorDatabase()
    if not vector_db.load():
        raise RuntimeError("Local FAISS index could not be loaded")

    retriever = Retriever(vector_db)
    results = retriever.retrieve(query, top_k=top_k)
    return results


def generate_boolean(query: str, top_k: int = 5, expand: bool = False) -> dict:
    """Programmatic entrypoint: returns a dict with suggestions and boolean query."""
    seed_species = extract_species_from_text(query)

    try:
        resp = call_search_api(query, top_k=top_k, expand=expand)
        results = resp.get('results', [])
    except Exception as e:
        # fallback to local retrieval
        try:
            results = local_retrieve(query, top_k=top_k)
        except Exception as e2:
            raise RuntimeError(f"Both API and local retrieval failed: {e2}")

    # build boolean
    boolean = build_boolean(results, seed_organisms=seed_species)

    # prepare suggestions list
    suggestions = []
    for r in results:
        suggestions.append({
            'query_text': r.get('query_text'),
            'query_type': r.get('query_type'),
            'similarity_score': r.get('similarity_score'),
            'kb_data': r.get('kb_data', {})
        })

    return {
        'input': query,
        'boolean_query': boolean,
        'suggestions': suggestions
    }


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    query = ' '.join(sys.argv[1:])
    seed_species = extract_species_from_text(query)

    try:
        resp = call_search_api(query, top_k=5, expand=False)
        results = resp.get('results', [])
    except Exception as e:
        print(f"API not available, falling back to local retrieval: {e}")
        try:
            results = local_retrieve(query, top_k=5)
        except Exception as e2:
            print(f"Local retrieval failed: {e2}")
            sys.exit(1)
    print('\nTop suggestions from semantic search:')
    for i, r in enumerate(results, 1):
        print(f" {i}. [{r.get('similarity_score'):.3f}] {r.get('query_text')} ({r.get('query_type')})")

    boolean = build_boolean(results, seed_organisms=seed_species)
    print('\nGenerated boolean query (NCBI-ready):\n')
    print(boolean)
    print('\n-- End')


if __name__ == '__main__':
    main()
