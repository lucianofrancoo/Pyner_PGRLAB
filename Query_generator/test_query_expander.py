#!/usr/bin/env python3
"""
Test del Query Expander — Flujo completo con LLM + Diccionario MeSH + Fallback PubMedBERT

Flujo:
  1. Input del usuario en lenguaje natural
  2. Ollama LLM extrae los keywords
  3. Lookup en final_synonym_dictionary.json (MeSH + variantes ortográficas)
  4. Si no encuentra → PubMedBERT mapea al MeSH más cercano → relookup en diccionario
  5. Genera boolean query NCBI

Uso:
  python3 test_query_expander.py
  python3 test_query_expander.py --no-llm
  python3 test_query_expander.py --no-fallback
  python3 test_query_expander.py "arabidopsis drought rna-seq"
"""

import sys
import argparse
import json
import time
import subprocess
import re
from pathlib import Path

# Asegurar que podemos importar query_expander
sys.path.insert(0, str(Path(__file__).parent / "phases/phase1/synonym_builder"))
from query_expander import QueryExpander

# ──────────────────────────────────────────────────────────────────────────────
# COLORES (igual que pyner_miner.sh)
# ──────────────────────────────────────────────────────────────────────────────
GREEN  = '\033[0;32m'
BLUE   = '\033[0;34m'
YELLOW = '\033[1;33m'
RED    = '\033[0;31m'
CYAN   = '\033[0;36m'
NC     = '\033[0m'

def c(text, color): return f"{color}{text}{NC}"
def sep(char="━", width=65): print(c(char * width, GREEN))
def header(text): sep(); print(c(text, GREEN)); sep()

# ──────────────────────────────────────────────────────────────────────────────
# PASO LLM: EXTRACCIÓN DE KEYWORDS con Ollama
# ──────────────────────────────────────────────────────────────────────────────
import os
OLLAMA_HOST    = os.getenv('OLLAMA_HOST', 'http://localhost:11434')
DEFAULT_MODEL  = os.getenv('OLLAMA_MODEL', 'qwen2.5:14b')

def ollama_is_available() -> bool:
    try:
        import requests
        r = requests.get(OLLAMA_HOST, timeout=2)
        return r.status_code == 200
    except Exception:
        return False

def extract_keywords_with_llm(user_input: str, model: str = DEFAULT_MODEL) -> list[str] | None:
    """
    Llama a Ollama para extraer los keywords científicos del input del usuario.
    Retorna lista de keywords, o None si Ollama no está disponible.
    """
    try:
        import requests
    except ImportError:
        return None

    prompt = f"""You are a scientific literature search assistant.
Your task: extract the key scientific concepts from the user's query to search PubMed.

Rules:
- Extract biological organisms, experimental conditions, molecular techniques, genes, tissues, diseases
- Return ONLY a JSON array of strings, nothing else
- Keep compound terms together (e.g., "RNA-Seq", "water scarcity", "drought stress")
- DO NOT extract generic generic nouns by themselves (e.g. "development", "expression", "response", "activity", "metabolism", "signaling", "accumulation"). Always keep them attached to their context (e.g. "wheat development", "gene expression", "stress response"). If a generic word lacks clear context, ignore it.
- 3-8 terms maximum
- Use English terms

User query: "{user_input}"

Return format: ["term1", "term2", "term3"]"""

    payload = {
        "model": model,
        "prompt": prompt,
        "stream": False,
        "options": {"temperature": 0.1, "num_predict": 200}
    }

    try:
        r = requests.post(f"{OLLAMA_HOST}/api/generate", json=payload, timeout=120)
        r.raise_for_status()
        response_text = r.json().get("response", "").strip()

        # Extraer JSON del response
        match = re.search(r'\[.*?\]', response_text, re.DOTALL)
        if match:
            keywords = json.loads(match.group())
            if isinstance(keywords, list) and all(isinstance(k, str) for k in keywords):
                clean_kws = [k.strip() for k in keywords if k.strip()]
                return postprocess_keywords(clean_kws, user_input)
    except Exception as e:
        print(c(f"  ⚠ LLM error: {e}", YELLOW))

    return None

def postprocess_keywords(keywords: list[str], user_input: str) -> list[str]:
    """Reconfigura keywords que sean demasiado generales recombinándolos con el contexto"""
    generic_terms = {
        "development", "response", "regulation", "activity", "expression", 
        "metabolism", "signaling", "biosynthesis", "growth", "accumulation",
        "mechanism", "tolerance", "resistance", "adaptation", "quality", "yield"
    }
    
    final_keywords = []
    user_words = user_input.split()
    user_lower = [w.lower() for w in user_words]
    
    for kw in keywords:
        kw_lower = kw.lower()
        if kw_lower in generic_terms:
            # Try to recombine with adjacent word in original query
            try:
                idx = user_lower.index(kw_lower)
                recombined = ""
                # Prefer left neighbor (e.g. "wheat development", "gene expression")
                if idx > 0 and user_lower[idx-1] not in generic_terms:
                    recombined = f"{user_words[idx-1]} {user_words[idx]}"
                # If no left neighbor, try right neighbor
                elif idx < len(user_words) - 1 and user_lower[idx+1] not in generic_terms:
                    recombined = f"{user_words[idx]} {user_words[idx+1]}"
                
                if recombined:
                    recombined = recombined.strip('.,;:')
                    final_keywords.append(recombined)
                    continue
            except ValueError:
                pass
            # Si es muy general y no tiene vecino claro, se descarta (no se añade).
            continue
            
        final_keywords.append(kw)
        
    seen = set()
    return [x for x in final_keywords if not (x.lower() in seen or seen.add(x.lower()))]

def parse_keywords_simple(user_input: str) -> list[str]:
    """Fallback sin LLM: split básico preservando comillas."""
    quoted = re.findall(r'"([^"]+)"', user_input)
    rest = re.sub(r'"[^"]+"', '', user_input)
    words = [w.strip() for w in rest.split() if w.strip() and len(w) > 2]
    return postprocess_keywords(quoted + words, user_input)


# ──────────────────────────────────────────────────────────────────────────────
# DISPLAY DEL RESULTADO
# ──────────────────────────────────────────────────────────────────────────────
def show_result(result: dict):
    for kw, data in result["expanded"].items():
        method = data.get("method", "N/A")
        
        if method == "compound_split":
            print(f"  {c('🧩', GREEN)} '{c(kw, CYAN)}' [{c('FRASE COMPUESTA', GREEN)}]")
            
            # Show org
            o_data = data["org_data"]
            print(f"     ➤ {c('Organismo:', BLUE)} {o_data.get('preferred_term', o_data.get('query'))}")
            if o_data.get("synonyms"):
                print(f"       • {len(o_data['synonyms'])} sinónimos:")
                for s in o_data["synonyms"][:25]:
                    src  = s.get("source", "?")
                    conf = s.get("confidence", "official")
                    print(f"         ◦ {s['term'][:55]:<55} [{src}, conf={conf}]")
                
            # Show proc
            p_data = data["proc_data"]
            if p_data.get("method") == "compound_process_bypassed":
                print(f"     ➤ {c('Proceso:', BLUE)} {p_data.get('query')}")
                print(f"       • {c('Nota:', BLUE)} Integrado en el bloque del organismo (omitiendo cláusula extra para maximizar recall)")
            elif p_data.get("preferred_term"):
                print(f"     ➤ {c('Proceso:', BLUE)} {p_data.get('preferred_term')}")
                print(f"       • {c('Nota:', BLUE)} Integrado en el bloque del organismo (omitiendo cláusula extra para maximizar recall)")
            elif p_data.get("suggestions"):
                print(f"     ➤ {c('Proceso (MeSH fallback):', BLUE)} {p_data['suggestions'][0]['term']}")
            else:
                print(f"     ➤ {c('Proceso:', BLUE)} {p_data.get('query')}")
            print()
            continue

        fallback_remapped = data.get("fallback_mesh_found", False)

        if data.get("mesh_bypassed_for_experimental_condition"):
            icon = c("🚧", YELLOW)
            label = c("MeSH EVADIDO (Condición Experimental No Clínica)", YELLOW)
        elif method == "exact_match":
            icon = c("✅", GREEN)
            label = c("EXACTO", GREEN)
        elif method == "exact_gene_match":
            icon = c("✅", GREEN)
            label = c("EXACTO (Genérico)", GREEN)
        elif method == "gene_bypassed":
            icon = c("🧬", GREEN)  # Green to indicate success of preserving gene/protein
            label = c("GEN DESCONOCIDO (MeSH Evadido)", GREEN)
        elif method == "chemical_bypassed":
            icon = c("🧪", GREEN)
            label = c("ELEMENTO/NUTRIENTE (MeSH Evadido)", GREEN)
        else:
            icon = c("❌", RED)
            label = c("NO ENCONTRADO", RED)

        print(f"  {icon} '{c(kw, CYAN)}' [{label}]")

        if data.get("preferred_term"):
            print(f"     {c('MeSH preferred:', BLUE)} {data['preferred_term']}")

        if data.get("synonyms"):
            n = len(data["synonyms"])
            print(f"     {c('Sinónimos:', BLUE)} {n} ({sum(1 for s in data['synonyms'] if s.get('source')=='mesh')} MeSH + {sum(1 for s in data['synonyms'] if s.get('source')=='text')} texto)")
            for s in data["synonyms"][:25]:
                src  = s.get("source", "?")
                conf = s.get("confidence", "official")
                print(f"       • {s['term'][:55]:<55} [{src}, conf={conf}]")
            if n > 25:
                print(f"       ... y {n - 25} más")
                
        if data.get("lexicon_found"):
            print(f"     {c('Expansiones Experimentales (LLM & Lexicon):', BLUE)} {len(data['lexicon_expansions'])}")
            for s in data["lexicon_expansions"]:
                print(f"       • {s[:55]:<55} [corpus_ngram, conf=high]")
        print()


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Test del Query Expander — Input NL → Keywords LLM → Boolean NCBI"
    )
    parser.add_argument("query", nargs="?", default=None,
                        help="Query en lenguaje natural (opcional, se pide si no se provee)")
    parser.add_argument("--no-llm", action="store_true",
                        help="Desactivar extracción de keywords vía LLM")
    parser.add_argument("--model", type=str, default=DEFAULT_MODEL,
                        help=f"Modelo Ollama a usar (default: {DEFAULT_MODEL})")
    parser.add_argument("--export-json", type=str, default=None,
                        help="Ruta para exportar el diccionario resultante en JSON")
    args = parser.parse_args()

    # Banner
    print()
    print(c("╔═══════════════════════════════════════════════════════════════╗", CYAN))
    print(c("║                                                               ║", CYAN))
    print(c(f"║              {c('🔬 PYNER QUERY EXPANDER TEST', BLUE)}{c('                     ║', CYAN)}", CYAN))
    print(c("║    Input NL → LLM Keywords → MeSH Dict → Boolean Query       ║", CYAN))
    print(c("║                                                               ║", CYAN))
    print(c("╚═══════════════════════════════════════════════════════════════╝", CYAN))
    print()

    # ── PASO 1: Input del usuario ────────────────────────────────────
    header(f"[1/4] {c('Ingrese su consulta de investigación', GREEN)}")
    print(c("  Ejemplos:", CYAN))
    print(f"    - \"Arabidopsis drought stress RNA-Seq\"")
    print(f"    - \"Tomato water scarcity transcriptomics\"")
    print(f"    - \"Wheat heat stress grain quality metabolomics\"")
    print()

    if args.query:
        user_input = args.query
        print(f"  Query: {c(user_input, BLUE)}")
    else:
        user_input = input(c("  Tu consulta: ", YELLOW)).strip()
        if not user_input:
            print(c("ERROR: Input vacío.", RED))
            sys.exit(1)

    print(f"\n  ✓ Query: {c(user_input, BLUE)}")

    # ── PASO 2: LLM extrae keywords ──────────────────────────────────
    print()
    header(f"[2/4] {c('Extracción de keywords', GREEN)}")

    use_llm = not args.no_llm
    keywords = None

    if use_llm:
        if ollama_is_available():
            print(c(f"  Ollama disponible → extrayendo keywords con {args.model}...", CYAN))
            t0 = time.time()
            keywords = extract_keywords_with_llm(user_input, args.model)
            elapsed = time.time() - t0
            if keywords:
                print(c(f"  ✓ Keywords extraídos en {elapsed:.1f}s:", GREEN))
                for kw in keywords:
                    print(f"      • {kw}")
            else:
                print(c("  ⚠ LLM no retornó keywords válidos — usando parsing básico", YELLOW))
        else:
            print(c("  ⚠ Ollama no disponible — usando parsing básico", YELLOW))
            print(c(f"    (Para activar LLM: iniciar Ollama con `ollama serve` y `ollama pull {args.model}`)", YELLOW))

    if not keywords:
        keywords = parse_keywords_simple(user_input)
        print(c("  Parsing básico:", CYAN))
        for kw in keywords:
            print(f"      • {kw}")

    print(f"\n  ✓ Total keywords: {c(str(len(keywords)), BLUE)}")

    # ── PASO 3: Inicializar expander ─────────────────────────────────
    print()
    header(f"[3/4] {c('Expansión de sinónimos', GREEN)}")

    print(c(f"  Diccionario MeSH + variantes ortográficas...", CYAN))
    t0 = time.time()
    expander = QueryExpander()
    print(c(f"  ✓ Listo en {time.time()-t0:.1f}s", GREEN))
    print()

    # ── PASO 4: Expansión y boolean ──────────────────────────────────
    t0 = time.time()
    result = expander.expand_query(keywords)
    elapsed = time.time() - t0

    show_result(result)

    # Boolean query final
    header(f"[4/4] {c('Boolean Query NCBI', GREEN)}")
    print()
    print("NCBI Query:")
    print(f"{result['boolean_query']}")
    print()
    
    # Visual format
    print("Visual representation:")
    parts = result["boolean_query"].split(" AND ")
    if len(parts) > 1:
        print(f"  {c('AND', YELLOW)} ".join(f"  {p}" for p in parts))
    else:
        print(f"  {result['boolean_query']}")

    print()
    s_exact    = sum(1 for d in result["expanded"].values() if d.get("method") == "exact_match")
    s_none     = sum(1 for d in result["expanded"].values() if d.get("method") not in ("exact_match", "exact_gene_match"))
    print(c(f"  ⏱ Tiempo total: {elapsed*1000:.0f} ms", CYAN))
    print(c(f"  📊 Stats: {s_exact} exactos | {s_none} sinónimos base/no encontrados", CYAN))
    print()
    sep("═")
    print()

    if args.export_json:
        try:
            with open(args.export_json, 'w', encoding='utf-8') as f:
                json.dump(result, f, ensure_ascii=False, indent=2)
        except Exception as e:
            print(c(f"Error exportando JSON: {e}", RED))

if __name__ == "__main__":
    main()
