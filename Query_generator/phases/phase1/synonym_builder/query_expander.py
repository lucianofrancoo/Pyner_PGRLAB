#!/usr/bin/env python3
"""
Query Expander with Semantic Fallback

Strategy:
  1. Try exact match in pre-computed dictionary (O(1), fast)
  2. If not found, use PubMedBERT to find closest MeSH term (semantic)
  3. Return expanded synonyms for boolean query

Usage:
    expander = QueryExpander()

    # User query: "drought stress in tomato"
    keywords = ["drought stress", "tomato"]

    for keyword in keywords:
        synonyms = expander.expand(keyword)
        # → Returns MeSH terms + synonyms
"""
import json
import os
import requests
import re
from pathlib import Path
from typing import List, Dict, Tuple, Optional

class QueryExpander:
    def __init__(self, dictionary_path: Optional[str] = None):
        """
        Initialize Query Expander

        Args:
            dictionary_path: Path to final_synonym_dictionary.json
        """

        # Definir sufijos y prefijos que indican un sesgo muy clínico
        self.clinical_terms = {
            "disease", "disorder", "syndrome", "hypersensitivity", "allergy",
            "illness", "deficiency", "pathology", "injury", "neoplasm", "cancer"
        }

        # Paths
        script_dir = Path(__file__).parent
        self.dict_path = Path(dictionary_path) if dictionary_path else \
                        script_dir / "output" / "final_synonym_dictionary.json"

        # Load pre-computed dictionary
        print("Loading synonym dictionary...")
        with open(self.dict_path) as f:
            data = json.load(f)
            self.dictionary = data["dictionary"]
            self.metadata = data["metadata"]

        print(f"✓ Loaded {len(self.dictionary):,} terms")

        # Lazy load context para Gene Dictionary
        self.gene_dict_path = script_dir / "output" / "final_gene_dictionary.json"
        self.gene_dictionary = None
        # Load Experimental Conditions Lexicon
        self.lexicon_path = script_dir / "output" / "experimental_conditions_lexicon.json"
        self.lexicon = {}
        if self.lexicon_path.exists():
            print("Loading experimental conditions lexicon...")
            with open(self.lexicon_path, "r", encoding="utf-8") as f:
                self.lexicon = json.load(f)
            print(f"✓ Loaded {len(self.lexicon):,} n-gram phrases")
        else:
            print("⚠ Warning: Lexicon not found. Experimental condition expansion disabled.")

        # Experimental triggers
        self.condition_cores = {
            "stress", "stresses", "exposure", "treatment", "deficiency", 
            "deprivation", "starvation", "toxicity", "challenge", "infection", 
            "injury", "hypoxia", "anoxia", "irradiation", "salinity",
            "limitation", "depletion", "restriction", "damage", "perturbation", "inflammation",
            "drought", "heat", "cold", "salt", "freezing", "chilling"
        }

        self.condition_modifiers = {
            "heat", "cold", "freezing", "chilling", "temperature", "temperatures",
            "thermal", "radiation", "uv", "ultraviolet", "drug", "drugs", "chemical",
            "salinity", "salt", "water", "drought", "light", "dark", "darkness",
            "osmotic", "oxidative", "metal", "heavy", "biotic", "abiotic", "nutrient", 
            "nutrient deficiency", "starvation", "ozone", "mechanical", "wounding",
            "oxygen", "glucose", "iron", "phosphate", "pressure", "shear",
            "immune", "viral", "bacterial", "toxin", "compound"
        }

        self.compatible_words = {
            "response", "responses", "tolerance", "tolerant", "resistance", "resistant",
            "sensitive", "sensitivity", "adaptation", "adapted", "acclimation", "acclimated",
            "signaling", "pathway", "induced", "induction", "mediated", "regulated", "regulation",
            "damage", "damaged", "perturbation", "inflammatory"
        }

        self.chemical_elements = {
            "zinc", "nitrogen", "iron", "copper", "phosphorus", "potassium", "calcium",
            "magnesium", "sulfur", "oxygen", "carbon", "hydrogen", "sodium", "chlorine",
            "manganese", "boron", "molybdenum", "nickel", "cobalt", "selenium", "silicon",
            "cadmium", "lead", "arsenic", "mercury", "aluminum", "aluminium",
            "nutrient", "nutrients", "metal", "metals", "mineral", "minerals",
            "phosphate", "nitrate", "sulfate", "ammonium"
        }



    def _load_gene_dictionary(self):
        """Lazy load of the gene dictionary to save memory/time on non-gene queries."""
        if self.gene_dictionary is not None:
            return  # Actually already loaded
            
        if self.gene_dict_path.exists():
            print("\n  [LazyLoad] Initializing gene dictionary into memory...")
            with open(self.gene_dict_path, "r", encoding="utf-8") as f:
                self.gene_dictionary = json.load(f)
            print(f"  ✓ Cached {len(self.gene_dictionary):,} genes successfully!")
        else:
            self.gene_dictionary = {}
            print("\n  ⚠ Warning: final_gene_dictionary.json not found. Gene expansion disabled.")

    def expand(self, keyword: str, include_metadata: bool = False) -> List[str]:
        """
        Expand a keyword into synonyms

        Args:
            keyword: User's search term
            include_metadata: Return dict with sources and confidence

        Returns:
            List of synonyms (or dicts if include_metadata=True)
        """
        keyword_lower = keyword.lower()

        # Try exact match first (O(1))
        if keyword_lower in self.dictionary:
            entry = self.dictionary[keyword_lower]

            if include_metadata:
                return {
                    "query": keyword,
                    "found": True,
                    "method": "exact_match",
                    "preferred_term": entry["preferred_term"],
                    "synonyms": entry["synonyms"]
                }
            else:
                # Return unique terms only
                terms = set([entry["preferred_term"]])
                for syn in entry["synonyms"]:
                    terms.add(syn["term"])
                return list(terms)

        # Check if keyword is a gene/protein candidate to avoid catastrophic MeSH remapping
        # Typical patterns: WRKY33, TP53, HSP70, BRCA1, IL-33, NF-kB
        is_gene_candidate = False
        parts_no_space = keyword.replace(" ", "")
        if re.match(r'^[a-zA-Z]+[-_]?\d+[a-zA-Z0-9]*$', parts_no_space) or (keyword.isupper() and 2 <= len(keyword) <= 12):
            is_gene_candidate = True

        if is_gene_candidate:
            # 1. Lazy load gene database if not cached yet
            self._load_gene_dictionary()
            
            # 2. Try exact match in GENE dictionary (O(1))
            if keyword_lower in self.gene_dictionary:
                entry = self.gene_dictionary[keyword_lower]

                if include_metadata:
                    return {
                        "query": keyword,
                        "found": True,
                        "method": "exact_gene_match",
                        "preferred_term": entry["preferred_symbol"],
                        "synonyms": [{"term": s, "source": "gene_db", "confidence": "official"} for s in entry["synonyms"]]
                    }
                else:
                    # Return unique terms only
                    terms = set([entry["preferred_symbol"]])
                    for syn in entry["synonyms"]:
                        terms.add(syn)
                    return list(terms)

            # 3. Not found in main dict nor gene dict: bypass semantic fallback entirely
            if include_metadata:
                # Provide minor spelling variants automatically (e.g. without hyphen)
                minor_syns = []
                cleaned = keyword.replace('-', '').replace('_', '')
                if cleaned != keyword:
                    minor_syns.append(cleaned)
                return {
                    "query": keyword,
                    "found": False,
                    "method": "gene_bypassed",
                    "suggestions": [],
                    "synonyms": [{"term": s, "source": "heuristic", "confidence": "minor_variant"} for s in minor_syns]
                }
            else:
                terms_list = [keyword]
                cleaned = keyword.replace('-', '').replace('_', '')
                if cleaned != keyword:
                    terms_list.append(cleaned)
                return terms_list

        # Check for chemical elements / basic nutrients to prevent specific chemical class remappings
        is_chemical = False
        parts_no_space = keyword_lower.replace(" ", "")
        if keyword_lower in self.chemical_elements or parts_no_space in self.chemical_elements:
            is_chemical = True
            
        if is_chemical:
            # Evadir fallback semántico que ensucia con "Compounds", "Receptors", "Complexes"
            if include_metadata:
                minor_syns = [f"{keyword} ion", f"{keyword} ions"]
                return {
                    "query": keyword,
                    "found": False,
                    "method": "chemical_bypassed",
                    "suggestions": [],
                    "synonyms": [{"term": s, "source": "heuristic", "confidence": "high"} for s in minor_syns]
                }
            else:
                return [keyword, f"{keyword} ion", f"{keyword} ions"]

        # No match found
        if include_metadata:
            return {
                "query": keyword,
                "found": False,
                "method": None,
                "suggestions": []
            }
        else:
            return [keyword]  # Return original term


    def _llm_filter_candidates(self, keyword: str, candidates: List[str]) -> List[str]:
        """
        Use Qwen2.5/Ollama to semantically filter experimental conditions n-grams.
        Returns up to 8 highly relevant candidates.
        """
        prompt = f"""You are a specialized biology terminology assistant.
The user is searching for the experimental condition: "{keyword}".
Below is a pool of text n-grams mined from literature. Many are noisy or describe different kinds of stress/conditions.
Your task: Select up to 8 n-grams from the list that STRICTLY maintain the same conceptual meaning or specify the same condition as "{keyword}". Do NOT include unrelated conditions (e.g. if looking for "heat stress", do NOT include "oxidative stress", "posttraumatic stress", "virus stress", etc).
Return your selection as a JSON format with a single key "candidates" containing a list of strings. DO NOT output any other text or markdown.

List of candidates:
{json.dumps(candidates)}
"""
        host = os.environ.get("OLLAMA_HOST", "http://localhost:11434")
        model = os.environ.get("OLLAMA_MODEL", "qwen2.5:14b")
        
        try:
            response = requests.post(
                f"{host}/api/generate",
                json={
                    "model": model,
                    "prompt": prompt,
                    "format": "json",
                    "stream": False
                },
                timeout=120
            )
            response.raise_for_status()
            content = response.json()["response"]
            result = json.loads(content)
            if isinstance(result, dict) and "candidates" in result:
                return result["candidates"][:8]
            elif isinstance(result, list):
                return result[:8]
            return []
        except Exception as e:
            # Silently fallback to top 5 if LLM fails (e.g. timeout, service down)
            return candidates[:5]

    def expand_query(self, keywords: List[str]) -> Dict:
        """
        Expand a full query (multiple keywords)

        Args:
            keywords: List of user keywords

        Returns:
            Dict with expanded terms and boolean query
        """
        expanded = {}
        generic_processes = {
            "development", "response", "regulation", "activity", "expression", 
            "metabolism", "signaling", "biosynthesis", "growth", "accumulation",
            "mechanism", "tolerance", "resistance", "adaptation", "quality", "yield",
            "transcription", "translation", "degradation", "transport", "homeostasis"
        }

        for keyword in keywords:
            keyword_lower = keyword.lower()
            ordered_tokens = keyword_lower.split()
            
            is_compound = False
            org_part, proc_part = "", ""
            
            if len(ordered_tokens) >= 2:
                for i in range(1, len(ordered_tokens)):
                    left = " ".join(ordered_tokens[:i])
                    right = " ".join(ordered_tokens[i:])
                    
                    if left in self.dictionary and right in generic_processes:
                        entry = self.dictionary[left]
                        if any(s.get("source") == "ncbi_taxonomy" for s in entry.get("synonyms", [])):
                            is_compound = True
                            org_part, proc_part = left, right
                            break
                    elif right in self.dictionary and left in generic_processes:
                        entry = self.dictionary[right]
                        if any(s.get("source") == "ncbi_taxonomy" for s in entry.get("synonyms", [])):
                            is_compound = True
                            org_part, proc_part = right, left
                            break

            if is_compound:
                res_org = self.expand(org_part, include_metadata=True)
                res_proc = self.expand(proc_part, include_metadata=True)
                expanded[keyword] = {
                    "query": keyword,
                    "method": "compound_split",
                    "found": True,
                    "org_data": res_org,
                    "proc_data": res_proc
                }
            else:
                result = self.expand(keyword, include_metadata=True)
                expanded[keyword] = result

        # Build boolean query
        boolean_parts = []
        for keyword, data in expanded.items():
            if data.get("method") == "compound_split":
                terms_org = set([keyword, data["org_data"]["query"]])
                terms_proc = set([keyword, data["proc_data"]["query"]])
                
                # Expand org
                o_data = data["org_data"]
                if o_data.get("found"):
                    terms_org.add(o_data["preferred_term"])
                    for syn in o_data.get("synonyms", []):
                        terms_org.add(syn["term"])
                        
                # Expand proc
                p_data = data["proc_data"]
                
                # REGLA: Blindar el proceso/rasgo genérico dentro de una frase compuesta.
                is_p_bypassed = False
                if p_data.get("method") == "exact_match":
                    pr_term = p_data.get("preferred_term", "").lower()
                    pq_lower = p_data["query"].lower()
                    if pq_lower not in pr_term and pr_term not in pq_lower:
                        is_p_bypassed = True
                        
                if is_p_bypassed or not p_data.get("found"):
                    p_data["method"] = "compound_process_bypassed"
                    p_data["found"] = False
                    p_data.pop("suggestions", None)
                    p_data.pop("preferred_term", None)
                    p_data.pop("synonyms", None)
                else:
                    if p_data.get("found"):
                        terms_proc.add(p_data["preferred_term"])
                        for syn in p_data.get("synonyms", []):
                            terms_proc.add(syn["term"])
                            
                # Save into boolean query
                or_clause_org = " OR ".join([f'"{t}"' for t in sorted(terms_org)])
                
                # REGLA PEDIDA: El proceso ya está implícito en la frase/organismo,
                # omitir la cláusula del proceso para no dañar el recall.
                boolean_parts.append(f"({or_clause_org})")
                continue

            terms = set([keyword])  # ALWAYS include the original keyword

            keyword_lower = keyword.lower()
            keyword_tokens = set(keyword_lower.split())
            
            # Check if it's an experimental condition based on new rules
            is_experimental = False
            if keyword_tokens & self.condition_cores:
                is_experimental = True
            elif (keyword_tokens & self.condition_modifiers) and (keyword_tokens & self.compatible_words):
                is_experimental = True
                
            is_query_clinical = any(c in keyword_lower for c in self.clinical_terms)
            
            # RULE: Evitar remapeos MeSH clínicos ruidosos para condiciones experimentales no clínicas
            bypass_mesh = is_experimental and not is_query_clinical
            
            # Allow MeSH if it's practically equivalent (e.g. contains the same words)
            if bypass_mesh and data.get("preferred_term"):
                if keyword_lower in data["preferred_term"].lower() or data["preferred_term"].lower() in keyword_lower:
                    bypass_mesh = False

            if bypass_mesh:
                # Registramos en el metadata para mostrar que bloqueamos un mapeo clínico
                data["mesh_bypassed_for_experimental_condition"] = True
            else:
                if data["found"] and data["method"] == "exact_match":
                    # Exact match: use preferred term + all synonyms from dictionary
                    terms.add(data["preferred_term"])
                    for syn in data["synonyms"]:
                        terms.add(syn["term"])

                elif data["found"] and data["method"] == "exact_gene_match":
                    # Gene match: use preferred symbol + all synonyms
                    terms.add(data["preferred_term"])
                    for syn in data["synonyms"]:
                        terms.add(syn["term"])

                elif not data["found"] and data.get("method") == "gene_bypassed":
                    # Unknown gene: only original term + minor variants
                    for syn in data.get("synonyms", []):
                        terms.add(syn["term"])

                elif not data["found"] and data.get("method") == "chemical_bypassed":
                    # Unknown element/nutrient: original term + ion variants
                    for syn in data.get("synonyms", []):
                        terms.add(syn["term"])



            # Sort and create OR clause
            # -----------------------------------------------------------------
            # PASO ADICIONAL: EXPANSION LEXICA/SEMÁNTICA PARA CONDICIONES EXP.
            # -----------------------------------------------------------------
            lexicon_expansions = []
            if is_experimental and self.lexicon:
                # 1. Filter candidates lexically (must share at least one INFORMATIVE token with the keyword)
                informative_tokens = keyword_tokens - self.condition_cores
                if not informative_tokens:
                    # If the query is literally just "stress", any token is "informative" to avoid empty sets
                    informative_tokens = keyword_tokens
                    
                candidates_with_freq = []
                for phrase, freq in self.lexicon.items():
                    phrase_tokens = set(phrase.split())
                    if phrase_tokens & informative_tokens: # Shares at least one INFORMATIVE token
                        candidates_with_freq.append((phrase, freq))
                
                # 2. Sort by frequency and get top 20
                if candidates_with_freq:
                    candidates_with_freq.sort(key=lambda x: x[1], reverse=True)
                    top_20_candidates = [c[0] for c in candidates_with_freq[:20]]
                    
                    # 3. Semantic filtering via LLM to get the final 5-8 strictly related variants
                    filtered_candidates = self._llm_filter_candidates(keyword, top_20_candidates)
                    
                    for phrase_cand in filtered_candidates:
                        # Avoid adding exact clinical terms if original was not clinical
                        if any(c in phrase_cand.lower() for c in self.clinical_terms) and not is_query_clinical:
                            continue
                        lexicon_expansions.append(phrase_cand)
                        terms.add(phrase_cand)
                            
                if lexicon_expansions:
                    data["lexicon_expansions"] = lexicon_expansions
                    data["lexicon_found"] = True

            or_clause = " OR ".join([f'"{t}"' for t in sorted(terms)])
            boolean_parts.append(f"({or_clause})")

        # Join with AND
        boolean_query = " AND ".join(boolean_parts)

        return {
            "original_keywords": keywords,
            "expanded": expanded,
            "boolean_query": boolean_query
        }


# CLI interface for testing
def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Query Expander with Semantic Fallback",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Expand single term (in dictionary)
  python3 query_expander.py "drought"

  # Expand term not in dictionary (semantic fallback)
  python3 query_expander.py "water scarcity"

  # Expand multiple terms
  python3 query_expander.py "drought" "arabidopsis" "gene expression"

  # Full query expansion
  python3 query_expander.py --full-query "drought stress" "tomato" "transcriptomics"
        """
    )

    parser.add_argument(
        "terms",
        nargs="+",
        help="Keywords to expand"
    )
    parser.add_argument(
        "--full-query",
        action="store_true",
        help="Generate full boolean query"
    )
    parser.add_argument(
        "--no-fallback",
        action="store_true",
        help="Disable semantic fallback"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.80,
        help="Semantic similarity threshold (default: 0.80)"
    )

    args = parser.parse_args()

    # Initialize expander
    expander = QueryExpander(
        use_semantic_fallback=not args.no_fallback,
        semantic_threshold=args.threshold
    )

    print()
    print("="*80)
    print("QUERY EXPANSION")
    print("="*80)
    print()

    if args.full_query:
        # Expand full query
        result = expander.expand_query(args.terms)

        print("Original keywords:")
        for kw in result["original_keywords"]:
            print(f"  • {kw}")
        print()

        print("Expanded terms:")
        for keyword, data in result["expanded"].items():
            print(f"\n🔍 '{keyword}':")
            print(f"   Method: {data.get('method', 'N/A')}")

            if data["found"]:
                print(f"   Preferred: {data['preferred_term']}")
                print(f"   Synonyms: {len(data['synonyms'])}")
                for syn in data["synonyms"][:5]:
                    print(f"      • {syn['term']} ({syn['source']})")

            elif data.get("suggestions"):
                print(f"   Suggestions (semantic):")
                for sug in data["suggestions"]:
                    print(f"      • {sug['term']} (similarity={sug['similarity']})")
            else:
                print(f"   No matches found")
                
            if data.get("lexicon_found"):
                print(f"   Lexicon Expansions (Experimental Condition):")
                for lex in data["lexicon_expansions"][:10]:
                    print(f"      • {lex} [corpus_ngram, conf=high]")
                if len(data["lexicon_expansions"]) > 10:
                    print(f"      ... y {len(data['lexicon_expansions'])-10} más")

        print()
        print("="*80)
        print("BOOLEAN QUERY")
        print("="*80)
        print(result["boolean_query"])
        print()

    else:
        # Expand individual terms
        for term in args.terms:
            result = expander.expand(term, include_metadata=True)

            print(f"🔍 '{term}':")
            print(f"   Method: {result.get('method', 'N/A')}")

            if result["found"]:
                print(f"   Preferred: {result['preferred_term']}")
                print(f"   Synonyms ({len(result['synonyms'])}):")
                for syn in result["synonyms"][:10]:
                    conf = syn.get('confidence', 'N/A')
                    print(f"      • {syn['term']} ({syn['source']}, conf={conf})")

            elif result.get("suggestions"):
                print(f"   Semantic suggestions:")
                for sug in result["suggestions"]:
                    print(f"      • {sug['term']} (similarity={sug['similarity']})")
            else:
                print(f"   No matches")

            print()


if __name__ == "__main__":
    main()
