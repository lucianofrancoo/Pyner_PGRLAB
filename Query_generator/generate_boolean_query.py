#!/usr/bin/env python3
"""
Generate a boolean NCBI query from semantic search results.

Usage:
  python3 Query_generator/generate_boolean_query.py "Arabidopsis root drought"

The script calls the local API `/search`, collects the top results and
maps them to synonym groups to produce a boolean query ready for NCBI/Entrez.

IMPROVED VERSION: Stricter organism extraction with Ollama classification.
"""

import sys
import json
import re
import subprocess
from typing import List, Dict


COMMON_MAP = {
    'mouse': 'Mus musculus',
    'mice': 'Mus musculus',
    'rat': 'Rattus norvegicus',
    'rattus': 'Rattus norvegicus',
    'human': 'Homo sapiens',
    'homo sapiens': 'Homo sapiens',
    'arabidopsis': 'Arabidopsis thaliana',
}

# Organism search variants: {full_name: [variant1, variant2, ...]}
# Includes: full scientific name, genus, and common names
ORGANISM_VARIANTS = {
    'Arabidopsis thaliana': ['Arabidopsis', 'Arabidopsis thaliana'],
    'Homo sapiens': ['Homo', 'Homo sapiens', 'human'],
    'Mus musculus': ['Mus', 'Mus musculus', 'mouse'],
    'Rattus norvegicus': ['Rattus', 'Rattus norvegicus', 'rat'],
    'Danio rerio': ['Danio', 'Danio rerio', 'zebrafish'],
    'Caenorhabditis elegans': ['Caenorhabditis', 'Caenorhabditis elegans'],
    'Drosophila melanogaster': ['Drosophila', 'Drosophila melanogaster'],
    'Escherichia coli': ['Escherichia', 'Escherichia coli', 'E. coli'],
    'Saccharomyces cerevisiae': ['Saccharomyces', 'Saccharomyces cerevisiae'],
    'Oncorhyncus mykiss': ['Oncorhyncus', 'Oncorhyncus mykiss', 'rainbow trout'],
    'Bos taurus': ['Bos', 'Bos taurus', 'cattle', 'cow'],
    'Sus scrofa': ['Sus', 'Sus scrofa', 'pig', 'swine'],
    'Gallus gallus': ['Gallus', 'Gallus gallus', 'chicken'],
    'Canis lupus familiaris': ['Canis', 'Canis lupus familiaris', 'dog'],
    'Equus caballus': ['Equus', 'Equus caballus', 'horse'],
    'Oryza sativa': ['Oryza', 'Oryza sativa', 'rice'],
    'Zea mays': ['Zea', 'Zea mays', 'maize', 'corn'],
    'Triticum aestivum': ['Triticum', 'Triticum aestivum', 'wheat'],
    'Hordeum vulgare': ['Hordeum', 'Hordeum vulgare', 'barley'],
    'Solanum lycopersicum': ['Solanum', 'Solanum lycopersicum', 'tomato'],
    'Plasmodium falciparum': ['Plasmodium', 'Plasmodium falciparum'],
    'Mycobacterium tuberculosis': ['Mycobacterium', 'Mycobacterium tuberculosis'],
    'Streptococcus pneumoniae': ['Streptococcus', 'Streptococcus pneumoniae'],
    'Staphylococcus aureus': ['Staphylococcus', 'Staphylococcus aureus'],
}

STRATEGY_SYNS = [ '"RNA-Seq"', '"RNA sequencing"', 'transcriptome', 'rnaseq' ]
GENE_EXPR_SYNS = [ '"gene expression"', 'transcriptome', '"RNA-Seq"', 'expression' ]
CHIP_SYNS = [ '"ChIP-Seq"', '"ChIP sequencing"', '"chromatin immunoprecipitation"' ]
GUT_SYNS = [ 'gut', 'intestinal', 'fecal' ]
METAGENOME_SYNS = [ 'metagenome', 'metatranscriptome' ]

# Keywords that strongly suggest a term is NOT an organism
NON_ORGANISM_KEYWORDS = {
    'drought', 'stress', 'disease', 'infection', 'treatment', 'condition',
    'root', 'leaf', 'flower', 'tissue', 'organ', 'development', 'growth',
    'sequencing', 'rnaseq', 'chip-seq', 'strategy', 'method', 'analysis',
    'gene expression', 'transcriptome',
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


def load_kb_organisms() -> Dict[str, str]:
    """Load organisms from KB: maps lowercase -> original name."""
    try:
        from pathlib import Path
        kb_path = Path(__file__).resolve().parents[0] / 'phases' / 'phase1' / 'output' / 'stage3_knowledge_base.json'
        if not kb_path.exists():
            return {}
        with open(kb_path, 'r') as f:
            kb = json.load(f)
        orgs = {}
        for org_name in kb.get('organisms', {}).keys():
            orgs[org_name.lower()] = org_name
        return orgs
    except Exception as e:
        print(f"Warning: Could not load KB: {e}", file=sys.stderr)
        return {}


def classify_term_with_ollama(term: str) -> str:
    """Use Ollama to classify: is this an organism, search term, or ambiguous?
    
    Returns: 'organism', 'search_term', or 'ambiguous'
    """
    try:
        prompt = f"""Given the term: "{term}"

Is this primarily an organism name, a search/query term, or ambiguous?

Organism examples: "Arabidopsis thaliana", "Homo sapiens", "mouse gut metagenome"
Search term examples: "drought", "root tissue", "stress response", "gene expression analysis"

Respond with exactly ONE word: organism, search_term, or ambiguous"""
        
        cmd = [
            'curl', '-s', '-X', 'POST',
            'http://localhost:11434/api/generate',
            '-d', json.dumps({
                'model': 'llama2',
                'prompt': prompt,
                'stream': False,
            })
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        if proc.returncode == 0:
            try:
                resp = json.loads(proc.stdout)
                answer = resp.get('response', '').strip().lower()
                if 'organism' in answer:
                    return 'organism'
                elif 'search' in answer or 'term' in answer:
                    return 'search_term'
                else:
                    return 'ambiguous'
            except:
                return 'ambiguous'
    except:
        pass
    
    return 'ambiguous'


def extract_species_from_text(text: str, kb_orgs: Dict[str, str] = None) -> List[str]:
    """Extract organism names from query text — STRICT version.
    
    Strategy:
    1. Check KB for exact organism name matches
    2. Try common name mappings
    3. Use Ollama to validate ambiguous terms
    4. Return only validated organisms
    """
    if kb_orgs is None:
        kb_orgs = load_kb_organisms()
    
    text_l = text.lower()
    found = []
    
    # Strategy 1: Check for exact KB organism names (multi-word, longest first)
    sorted_orgs = sorted(kb_orgs.keys(), key=len, reverse=True)
    for org_lower in sorted_orgs:
        if org_lower in text_l:
            original_name = kb_orgs[org_lower]
            
            # Check if this looks like a false positive (e.g., "root associated fungus metagenome" in "Arabidopsis root drought")
            # by checking if search term keywords dominate
            search_term_score = sum(1 for kw in NON_ORGANISM_KEYWORDS if kw in text_l)
            
            # If this potential organism name contains non-organism keywords, skip it
            if any(kw in org_lower for kw in NON_ORGANISM_KEYWORDS if len(kw) > 3):  # Skip short keywords
                continue
            
            # For ambiguous cases, use Ollama to decide
            class_result = classify_term_with_ollama(original_name)
            if class_result == 'organism':
                if original_name not in found:
                    found.append(original_name)
            elif class_result == 'ambiguous' and search_term_score == 0:
                # If KB-known organism and no conflicting keywords, allow it
                if original_name not in found:
                    found.append(original_name)
    
    # Strategy 2: Check common name map
    for key, latin in COMMON_MAP.items():
        if re.search(rf"\b{re.escape(key)}\b", text_l):
            if latin not in found:
                if latin.lower() in kb_orgs or latin in kb_orgs.values():
                    found.append(latin)
    
    return found


def build_boolean(results: List[dict], seed_organisms: List[str] = None, kb_orgs: Dict[str, str] = None) -> str:
    """Build NCBI-ready boolean query from semantic search results.
    
    For each organism, includes all variants (genus, full name, common names).
    """
    if kb_orgs is None:
        kb_orgs = load_kb_organisms()
    
    organisms = []
    facets = []

    for r in results:
        qtext = r.get('query_text', '')
        qtype = r.get('query_type', '').lower()

        # Always fallback to quoted phrase
        raw_phrase = f'"{qtext}"'

        if qtype == 'organism':
            # Check if extracted species match KB (stricter validation)
            species = extract_species_from_text(qtext, kb_orgs)
            
            if species:
                # Found a validated organism in the text
                for s in species:
                    org_tag = _organism_with_variants(s)
                    if org_tag not in organisms:
                        organisms.append(org_tag)
            elif qtext.lower() in kb_orgs:
                # Direct KB match (exact organism name)
                org_name = kb_orgs[qtext.lower()]
                org_tag = _organism_with_variants(org_name)
                if org_tag not in organisms:
                    organisms.append(org_tag)
            else:
                # Low confidence: this phrase is likely NOT a pure organism
                # Skip it to avoid polluting the boolean query
                pass  # Do NOT add to facets; discard

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
            facets.append([raw_phrase])

    # Build clauses
    clause_list = []
    
    # Include seed organisms detected from the original query
    if seed_organisms:
        for s in seed_organisms:
            org_tag = _organism_with_variants(s)
            if org_tag not in organisms:
                organisms.append(org_tag)

    if organisms:
        clause_list.append('(' + ' OR '.join(organisms) + ')')

    # Deduplicate and merge facets
    uniq_facets = []
    seen = set()
    for facet in facets:
        key = tuple(sorted(facet))
        if key not in seen:
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


def _organism_with_variants(organism_name: str) -> str:
    """Build organism clause with all search variants (genus, full name, common names).
    
    Example: Arabidopsis thaliana -> '("Arabidopsis" OR "Arabidopsis thaliana")[Organism]'
    """
    variants = ORGANISM_VARIANTS.get(organism_name, [organism_name])
    quoted_variants = [f'"{v}"' for v in variants]
    return '(' + ' OR '.join(quoted_variants) + ')[Organism]'


def local_retrieve(query: str, top_k: int = 5):
    """Fallback: use local phase2 vector DB inside Query_generator/phases."""
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
    """Programmatic entrypoint: returns dict with suggestions and boolean query."""
    kb_orgs = load_kb_organisms()
    seed_species = extract_species_from_text(query, kb_orgs)

    try:
        resp = call_search_api(query, top_k=top_k, expand=expand)
        results = resp.get('results', [])
    except Exception as e:
        # Fallback to local retrieval
        try:
            results = local_retrieve(query, top_k=top_k)
        except Exception as e2:
            raise RuntimeError(f"Both API and local retrieval failed: {e2}")

    # Build boolean with KB context
    boolean = build_boolean(results, seed_organisms=seed_species, kb_orgs=kb_orgs)

    # Prepare suggestions list
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
    kb_orgs = load_kb_organisms()
    seed_species = extract_species_from_text(query, kb_orgs)

    try:
        resp = call_search_api(query, top_k=5, expand=False)
        results = resp.get('results', [])
    except Exception as e:
        print(f"API not available, falling back to local retrieval: {e}", file=sys.stderr)
        try:
            results = local_retrieve(query, top_k=5)
        except Exception as e2:
            print(f"Local retrieval failed: {e2}", file=sys.stderr)
            sys.exit(1)
    
    print('\nTop suggestions from semantic search:')
    for i, r in enumerate(results, 1):
        print(f" {i}. [{r.get('similarity_score'):.3f}] {r.get('query_text')} ({r.get('query_type')})")

    boolean = build_boolean(results, seed_organisms=seed_species, kb_orgs=kb_orgs)
    print('\nGenerated boolean query (NCBI-ready):\n')
    print(boolean)
    print('\n-- End')


if __name__ == '__main__':
    main()
