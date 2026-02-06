"""
Ollama LLM Integration
======================
Query expansion and result enhancement using Ollama LLMs
"""

import requests
import json
import logging
from typing import List, Dict, Any
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from phase3.config import OLLAMA_HOST, OLLAMA_MODEL, OLLAMA_TIMEOUT, QUERY_EXPANSION_ENABLED

logger = logging.getLogger(__name__)

class OllamaClient:
    """Cliente para interactuar con Ollama LLM"""
    
    def __init__(self):
        self.host = OLLAMA_HOST
        self.model = OLLAMA_MODEL
        self.timeout = OLLAMA_TIMEOUT
        logger.info(f"🤖 OllamaClient initialized: {self.host}/{self.model}")
    
    def health_check(self) -> bool:
        """Verificar que Ollama esté disponible"""
        try:
            response = requests.get(f"{self.host}/api/tags", timeout=5)
            return response.status_code == 200
        except:
            logger.warning(f"⚠️ Ollama not available at {self.host}")
            return False
    
    def expand_query(self, query: str, max_tokens: int = 100) -> List[str]:
        """Expandir query usando LLM"""
        if not QUERY_EXPANSION_ENABLED:
            return [query]
        
        prompt = f"""Given this scientific research query, generate 3 related search variations:
Query: {query}

Respond with ONLY 3 variations, one per line, no numbering or bullets."""
        
        try:
            response = requests.post(
                f"{self.host}/api/generate",
                json={
                    "model": self.model,
                    "prompt": prompt,
                    "stream": False,
                    "temperature": 0.7,
                },
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                result = response.json()
                text = result.get("response", "")
                variations = [v.strip() for v in text.split('\n') if v.strip()]
                logger.info(f"✅ Expanded '{query}' → {len(variations)} variations")
                return [query] + variations[:3]
            
        except Exception as e:
            logger.warning(f"⚠️ Query expansion failed: {e}")
        
        return [query]
    
    def rank_results(self, query: str, results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Re-rankear resultados usando LLM"""
        if len(results) <= 1:
            return results
        
        result_text = "\n".join([f"{i+1}. {r.get('query_text', '')}" for i, r in enumerate(results[:5])])
        prompt = f"""Given this research query and search results, rank them by relevance:

Query: {query}

Results:
{result_text}

Respond with ONLY the ranking numbers (e.g., "3,1,2,4,5"), in order of relevance."""
        
        try:
            response = requests.post(
                f"{self.host}/api/generate",
                json={
                    "model": self.model,
                    "prompt": prompt,
                    "stream": False,
                },
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                ranking_text = response.json().get("response", "")
                try:
                    ranking = [int(x.strip())-1 for x in ranking_text.split(',')[:len(results)]]
                    ranked = [results[i] for i in ranking if i < len(results)]
                    logger.info(f"✅ Re-ranked {len(ranked)} results")
                    return ranked
                except:
                    pass
        except Exception as e:
            logger.warning(f"⚠️ Ranking failed: {e}")
        
        return results


class QueryExpander:
    """Expande queries automáticamente"""
    
    def __init__(self):
        self.ollama = OllamaClient()
        self.fallback_templates = {
            "expression": ["gene expression", "tissue expression", "cell type expression"],
            "disease": ["pathology", "etiology", "molecular mechanisms"],
            "organism": ["genomics", "transcriptomics", "proteomics"],
        }
    
    def expand(self, query: str) -> List[str]:
        """Expandir query con o sin LLM"""
        # Try LLM first
        if self.ollama.health_check():
            return self.ollama.expand_query(query)
        
        # Fallback: heuristic expansion
        variations = [query]
        query_lower = query.lower()
        
        for category, templates in self.fallback_templates.items():
            if any(word in query_lower for word in ["gene", "expression"]):
                variations.extend([f"{query} {t}" for t in templates])
        
        return variations[:4]  # Return max 4 variations


# Test
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    expander = QueryExpander()
    
    test_queries = [
        "COVID-19 research",
        "Gene expression in humans",
        "Bacterial genomics"
    ]
    
    print("\n" + "="*80)
    print("🤖 OLLAMA INTEGRATION TEST")
    print("="*80)
    
    for q in test_queries:
        print(f"\n📝 Original: {q}")
        expanded = expander.expand(q)
        print(f"✓ Expanded ({len(expanded)} variants):")
        for var in expanded[1:]:
            print(f"  - {var}")
