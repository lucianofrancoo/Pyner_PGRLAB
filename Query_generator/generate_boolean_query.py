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

STRATEGY_SYNS = [ '"RNA-Seq"', '"RNA seq"', '"RNA sequencing"', 'transcriptome' ]
RNASEQ_SYNS = [ '"RNA-Seq"', '"RNA seq"', '"RNA sequencing"', '"whole transcriptome"' ]
WGS_SYNS = [ '"whole genome sequencing"', '"WGS"', '"whole genome"' ]
METAGENOMICS_SYNS = [ 'metagenomic', 'metagenomics', '"16S rRNA"', '"16S"' ]
CHIP_SYNS = [ '"ChIP-Seq"', '"ChIP seq"', '"chromatin immunoprecipitation"', '"ChIP sequencing"' ]
ATACSEQ_SYNS = [ '"ATAC-Seq"', '"ATAC seq"', '"assay for transposase-accessible chromatin"' ]
DNASEQ_SYNS = [ '"DNase-Seq"', '"DNase seq"', '"digital genomic footprinting"' ]
WES_SYNS = [ '"whole exome sequencing"', '"WES"', '"exome sequencing"' ]
TARGETSEQ_SYNS = [ '"targeted sequencing"', '"amplicon sequencing"', '"deep sequencing"' ]

# Gene expression specific synonyms - SEPARATE from sequencing strategies
GENE_EXPR_SYNS = [ '"gene expression"', '"expression profiling"', '"expression analysis"', 'microarray' ]
QPCR_SYNS = [ '"qPCR"', '"q-PCR"', '"quantitative PCR"', '"RT-qPCR"', '"real-time PCR"' ]

GUT_SYNS = [ 'gut', 'intestinal', 'fecal' ]
METAGENOME_SYNS = [ 'metagenome', 'metatranscriptome' ]

# Keywords that strongly suggest a term is NOT an organism
NON_ORGANISM_KEYWORDS = {
    'drought', 'stress', 'disease', 'infection', 'treatment', 'condition',
    'root', 'leaf', 'flower', 'tissue', 'organ', 'development', 'growth',
    'sequencing', 'rnaseq', 'chip-seq', 'strategy', 'method', 'analysis',
    'gene expression', 'transcriptome',
}

# Query type relevance mapping: what user search terms should match what result types
# If user searches for these terms, they're probably interested in these query types
RELEVANCE_MAP = {
    # Stress/environment related
    'drought': ['organism', 'strategy'],  # NOT gene_expression generally
    'stress': ['organism', 'strategy'],
    'treatment': ['organism', 'strategy'],
    'condition': ['organism', 'strategy'],
    'disease': ['organism', 'strategy'],
    'infection': ['organism', 'strategy'],
    'immune': ['organism', 'strategy'],
    'pathogen': ['organism', 'strategy'],
    
    # Development/anatomy
    'root': ['organism', 'strategy'],
    'leaf': ['organism', 'strategy'],
    'development': ['organism', 'strategy'],
    'growth': ['organism', 'strategy'],
    'morphogenesis': ['organism', 'strategy'],
    
    # Expression analysis
    'expression': ['gene_expression', 'strategy'],
    'profiling': ['gene_expression', 'strategy'],
    'microarray': ['gene_expression'],
    'transcriptome': ['strategy'],
    'rna': ['strategy'],
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


def extract_query_keywords(query: str) -> set:
    """Extract key search terms from user query.
    
    Example: "Arabidopsis root drought" -> {'arabidopsis', 'root', 'drought'}
    """
    words = query.lower().split()
    keywords = set()
    for word in words:
        # Remove punctuation and add
        clean_word = re.sub(r'[^\w]', '', word)
        if len(clean_word) > 2:  # Skip short words
            keywords.add(clean_word)
    return keywords


def get_relevant_query_types(user_query: str) -> set:
    """Determine which result types are relevant for this user query.
    
    If user searches for 'drought', they probably want organism + strategy results,
    NOT gene_expression results.
    """
    keywords = extract_query_keywords(user_query)
    relevant_types = {'organism', 'strategy'}  # Always include organism and strategy
    
    for keyword in keywords:
        if keyword in RELEVANCE_MAP:
            relevant_types.update(RELEVANCE_MAP[keyword])
    
    return relevant_types


def filter_results_by_relevance(results: List[dict], user_query: str) -> List[dict]:
    """Filter semantic results to keep only those relevant to user query.
    
    If user searches for 'drought' but gets 'gene expression' results, discard them.
    """
    relevant_types = get_relevant_query_types(user_query)
    keywords = extract_query_keywords(user_query)
    
    filtered = []
    for r in results:
        qtype = r.get('query_type', '').lower()
        qtext = r.get('query_text', '').lower()
        
        # Always keep organism type
        if qtype == 'organism':
            filtered.append(r)
            continue
        
        # Keep strategy if it's relevant
        if qtype == 'strategy':
            filtered.append(r)
            continue
        
        # For gene_expression: only keep if user explicitly searched for it
        if qtype == 'gene_expression':
            expr_keywords = {'expression', 'profiling', 'microarray', 'transcriptome'}
            if any(kw in keywords for kw in expr_keywords):
                filtered.append(r)
            # Otherwise skip gene_expression results when user didn't ask for them
            continue
        
        # Keep other types by default
        filtered.append(r)
    
    return filtered


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
            # Detect specific sequencing strategy to use appropriate synonyms
            qtext_lower = qtext.lower()
            
            if 'chip' in qtext_lower or 'chromatin' in qtext_lower or 'immunoprecipitation' in qtext_lower:
                facets.append(CHIP_SYNS)
            elif 'atac' in qtext_lower or 'transposase' in qtext_lower:
                facets.append(ATACSEQ_SYNS)
            elif 'dnase' in qtext_lower or 'dnasei' in qtext_lower:
                facets.append(DNASEQ_SYNS)
            elif 'rna-seq' in qtext_lower or 'rna seq' in qtext_lower or 'rnaseq' in qtext_lower or 'rna sequencing' in qtext_lower:
                facets.append(RNASEQ_SYNS)
            elif 'wgs' in qtext_lower or 'whole genome sequencing' in qtext_lower or ('whole' in qtext_lower and 'genome' in qtext_lower):
                facets.append(WGS_SYNS)
            elif 'wes' in qtext_lower or 'exome' in qtext_lower or ('whole' in qtext_lower and 'exome' in qtext_lower):
                facets.append(WES_SYNS)
            elif 'qpcr' in qtext_lower or 'q-pcr' in qtext_lower or 'real-time pcr' in qtext_lower or 'qt-pcr' in qtext_lower:
                facets.append(QPCR_SYNS)
            elif '16s' in qtext_lower or 'metagenomic' in qtext_lower or '16s rrna' in qtext_lower:
                facets.append(METAGENOMICS_SYNS)
            elif 'amplicon' in qtext_lower or 'targeted' in qtext_lower:
                facets.append(TARGETSEQ_SYNS)
            else:
                # Generic strategy: use raw phrase
                facets.append([raw_phrase])

        elif qtype == 'gene_expression':
            # Gene expression analysis - different from sequencing strategies
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
    
    If organism not in ORGANISM_VARIANTS, asks Ollama to generate variants.
    """
    variants = get_organism_variants(organism_name)
    quoted_variants = [f'"{v}"' for v in variants]
    return '(' + ' OR '.join(quoted_variants) + ')[Organism]'


# Cache for ALL variants (both dict and Ollama-generated)
_VARIANTS_CACHE = {}


def get_organism_variants(organism_name: str) -> List[str]:
    """Get search variants for an organism.
    
    First checks ORGANISM_VARIANTS dict.
    If not found, calls Ollama to generate variants based on pattern.
    Caches all results to avoid repeated lookups.
    """
    # Check cache first (covers both dict-based and Ollama-generated)
    if organism_name in _VARIANTS_CACHE:
        return _VARIANTS_CACHE[organism_name]
    
    # Check if already in hardcoded dict
    if organism_name in ORGANISM_VARIANTS:
        result = ORGANISM_VARIANTS[organism_name]
        _VARIANTS_CACHE[organism_name] = result
        return result
    
    # Call Ollama to generate variants
    try:
        variants = _ask_ollama_for_variants(organism_name)
        if variants:
            _VARIANTS_CACHE[organism_name] = variants
            return variants
    except Exception as e:
        print(f"Warning: Ollama variant generation failed for {organism_name}: {e}", file=sys.stderr)
    
    # Fallback: return just the organism name
    result = [organism_name]
    _VARIANTS_CACHE[organism_name] = result
    return result


def _ask_ollama_for_variants(organism_name: str) -> List[str]:
    """Use Ollama to generate organism search variants (genus, common names, etc.).
    
    Returns: list of variant names, e.g. ["Arabidopsis", "Arabidopsis thaliana"]
    """
    # Build examples from the known dict
    examples = []
    for org, vars_list in list(ORGANISM_VARIANTS.items())[:5]:  # First 5 as examples
        examples.append(f"  {org} -> {vars_list}")
    examples_str = '\n'.join(examples)
    
    prompt = f"""You are an expert in biology and scientific nomenclature.

Given an organism scientific name, generate search variants that would be useful for finding research studies.
Variants should include:
1. The genus name (first word)
2. The full scientific name
3. Common names if they exist (e.g., "mouse" for "Mus musculus", "rainbow trout" for "Oncorhyncus mykiss")

EXAMPLES OF EXPECTED OUTPUT FORMAT:
{examples_str}

Now, for the organism: "{organism_name}"

Generate variants in the SAME format as above. Return ONLY the organism name and its variants list, nothing else.
Format: organism_name -> [variant1, variant2, variant3, ...]"""

    try:
        cmd = [
            'curl', '-s', '-X', 'POST',
            'http://localhost:11434/api/generate',
            '-d', json.dumps({
                'model': 'llama2',
                'prompt': prompt,
                'stream': False,
            })
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        if proc.returncode == 0:
            resp = json.loads(proc.stdout)
            output = resp.get('response', '').strip()
            
            # Parse output: "organism_name -> [variant1, variant2, ...]"
            variants = _parse_ollama_variants(output)
            if variants:
                return variants
    except Exception as e:
        raise RuntimeError(f"Ollama call failed: {e}")
    
    return None


def _parse_ollama_variants(output: str) -> List[str]:
    """Parse Ollama output to extract variants list.
    
    Expected format: "Organism Name -> [variant1, variant2, variant3]"
    Returns: [variant1, variant2, variant3]
    """
    try:
        # Find the arrow and extract the list part
        if '->' not in output:
            return None
        
        list_part = output.split('->', 1)[1].strip()
        
        # Try to parse as Python list
        # Replace single quotes with double quotes for JSON compatibility
        list_part = list_part.replace("'", '"')
        
        # Extract content between brackets
        if list_part.startswith('[') and list_part.endswith(']'):
            variants = json.loads(list_part)
            if isinstance(variants, list) and all(isinstance(v, str) for v in variants):
                return variants
    except Exception as e:
        pass
    
    return None


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
    """Programmatic entrypoint: returns dict with suggestions and boolean query.
    
    Results are filtered for relevance to the user's query keywords.
    """
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

    # Filter results for relevance to the actual user query
    filtered_results = filter_results_by_relevance(results, query)
    
    # Build boolean with filtered, relevant results
    boolean = build_boolean(filtered_results, seed_organisms=seed_species, kb_orgs=kb_orgs)

    # Prepare suggestions list (show what was actually used, not all raw results)
    suggestions = []
    for r in filtered_results:
        suggestions.append({
            'query_text': r.get('query_text'),
            'query_type': r.get('query_type'),
            'similarity_score': r.get('similarity_score'),
            'kb_data': r.get('kb_data', {})
        })

    return {
        'input': query,
        'boolean_query': boolean,
        'suggestions': suggestions,
        'note': f'Filtered {len(results) - len(filtered_results)} irrelevant results'
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
    
    # Filter for relevance
    filtered_results = filter_results_by_relevance(results, query)
    
    print('\nAll semantic search results:')
    for i, r in enumerate(results, 1):
        is_filtered = r not in filtered_results
        status = "❌ FILTERED (irrelevant)" if is_filtered else "✅ KEPT"
        print(f" {i}. [{r.get('similarity_score'):.3f}] {r.get('query_text')} ({r.get('query_type')}) {status}")
    
    if len(filtered_results) < len(results):
        print(f"\n⚠️  Filtered out {len(results) - len(filtered_results)} irrelevant results")

    print(f'\n📊 Using {len(filtered_results)} relevant results to build boolean query')
    boolean = build_boolean(filtered_results, seed_organisms=seed_species, kb_orgs=kb_orgs)
    print('\nGenerated boolean query (NCBI-ready):\n')
    print(boolean)
    print('\n-- End')


if __name__ == '__main__':
    main()
