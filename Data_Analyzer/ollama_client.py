"""
Ollama Client for Paper Analysis
=================================
Wrapper for Ollama API calls specific to paper classification
"""

import json
import logging
import requests
from typing import Dict, Optional, List
from config import OLLAMA_BASE_URL, OLLAMA_MODEL, OLLAMA_TIMEOUT

logger = logging.getLogger(__name__)


class OllamaClient:
    """Client for Ollama LLM API"""
    
    def __init__(self, base_url: str = OLLAMA_BASE_URL, model: str = OLLAMA_MODEL):
        self.base_url = base_url.rstrip('/')
        self.model = model
        self.timeout = OLLAMA_TIMEOUT
    
    def is_available(self) -> bool:
        """Check if Ollama server is running"""
        try:
            response = requests.get(f"{self.base_url}/api/tags", timeout=5)
            return response.status_code == 200
        except requests.RequestException:
            return False
    
    def analyze_paper(self, title: str, abstract: str, user_query: str, full_text: Optional[str] = None) -> Dict:
        """
        Analyze a paper and extract structured information
        
        Args:
            title: Paper title
            abstract: Paper abstract
            user_query: Original user query for relevance evaluation
            full_text: Optional full text from PMC (Methods/Results sections)
            
        Returns:
            Dictionary with extracted information
        """
        prompt = self._build_analysis_prompt(title, abstract, user_query, full_text)
        
        try:
            response = self._call_ollama(prompt)
            parsed = self._parse_response(response)
            return parsed
        except Exception as e:
            logger.error(f"Error analyzing paper: {e}")
            return self._empty_result()
    
    def _build_analysis_prompt(self, title: str, abstract: str, user_query: str, full_text: Optional[str] = None) -> str:
        """Build prompt for paper analysis"""
        
        # If full text available, use it for better analysis
        if full_text:
            content_source = f"""
PAPER ABSTRACT: {abstract}

PAPER FULL TEXT (Methods & Results sections):
{full_text}

NOTE: Use both abstract and full text for comprehensive analysis. The full text contains detailed information about organisms, tissues, experimental conditions, and techniques. Look especially in Methods section for experimental techniques.
"""
        else:
            content_source = f"""
PAPER ABSTRACT: {abstract}

NOTE: Only abstract available. Extract what you can from the abstract.
"""
        
        return f"""You are a scientific paper classifier. Analyze the following paper and extract structured information.

USER QUERY (for relevance): {user_query}

PAPER TITLE: {title}

{content_source}

YOUR TASK:
1. Evaluate if this paper is RELEVANT to the user query. Score from 0-10 (0=completely irrelevant, 10=perfectly relevant)
2. Extract ALL organisms mentioned (scientific name preferred, e.g., "Solanum lycopersicum", "Mus musculus")
3. Extract ALL tissues/organs mentioned (e.g., "root", "leaf", "liver", "brain")
4. Extract ALL experimental conditions/treatments (e.g., "drought stress", "heat shock", "nitrogen deficiency")
5. Extract ALL experimental strategies/techniques used in the study. Look for:
   - Molecular techniques: qRT-PCR, RT-PCR, RNA-Seq, DNA microarray, Northern blot, Western blot, ELISA
   - Imaging: confocal microscopy, fluorescence microscopy, histological staining
   - Sequencing: sequencing, Next-Gen sequencing
   - Phenotyping: morphological analysis, gas exchange measurements
   - Biochemical: enzyme assays, antioxidant assays, chromatography
   - Other: bioinformatics, statistical analysis, model plants
   NOTE: If Methods section mentions specific techniques, include them even if not explicitly mentioned in abstract.

RESPONSE FORMAT (JSON only, no explanation):
{{
  "relevance_score": <number 0-10>,
  "organisms": ["organism1", "organism2"],
  "tissues": ["tissue1", "tissue2"],
  "conditions": ["condition1", "condition2"],
  "strategies": ["strategy1", "strategy2"]
}}

IMPORTANT:
- Return ONLY valid JSON
- Use empty arrays [] if category not found (but try harder to find!)
- Use scientific names for organisms when possible
- Be precise and specific
- Include synonyms if mentioned (e.g., "tomato" and "Solanum lycopersicum")
- For strategies: be comprehensive. Include ALL techniques mentioned in Methods/Results
- Never return "N/A" for strategies - use empty array [] if truly no techniques found
- Grammar: "qRT-PCR" (not "qrt-pcr"), "RNA-Seq" (not "RNA-seq")

JSON:"""
    
    def _call_ollama(self, prompt: str) -> str:
        """Call Ollama API"""
        url = f"{self.base_url}/api/generate"
        
        # Dynamically increase timeout for longer texts
        # Longer prompts take more time for LLM to process
        # Estimate: ~1 second per 2000 characters + base timeout
        timeout = max(self.timeout, int(len(prompt) / 2000 * 60))  # ~1 min per 2000 chars
        logger.info(f"Prompt length: {len(prompt)} chars, timeout: {timeout}s")
        
        payload = {
            "model": self.model,
            "prompt": prompt,
            "stream": False,
            "options": {
                "temperature": 0.1,  # Low temperature for consistent extraction
                "top_p": 0.9,
                "num_predict": 500  # Limit response length
            }
        }
        
        logger.info(f"Calling Ollama API: {url}")
        response = requests.post(url, json=payload, timeout=timeout)
        response.raise_for_status()
        
        result = response.json()
        return result.get('response', '')
    
    def _parse_response(self, response: str) -> Dict:
        """Parse Ollama JSON response"""
        try:
            # Try to find JSON in response
            start = response.find('{')
            end = response.rfind('}') + 1
            
            if start == -1 or end == 0:
                logger.warning("No JSON found in response")
                return self._empty_result()
            
            json_str = response[start:end]
            parsed = json.loads(json_str)
            
            # Validate and sanitize
            result = {
                'relevance_score': int(parsed.get('relevance_score', 0)),
                'organisms': parsed.get('organisms', []),
                'tissues': parsed.get('tissues', []),
                'conditions': parsed.get('conditions', []),
                'strategies': parsed.get('strategies', [])
            }
            
            # Ensure all are lists
            for key in ['organisms', 'tissues', 'conditions', 'strategies']:
                if not isinstance(result[key], list):
                    result[key] = []
            
            # Clamp relevance score
            result['relevance_score'] = max(0, min(10, result['relevance_score']))
            
            return result
            
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse JSON: {e}")
            logger.debug(f"Response was: {response}")
            return self._empty_result()
        except Exception as e:
            logger.error(f"Unexpected error parsing response: {e}")
            return self._empty_result()
    
    def _empty_result(self) -> Dict:
        """Return empty result structure"""
        return {
            'relevance_score': 0,
            'organisms': [],
            'tissues': [],
            'conditions': [],
            'strategies': []
        }


def test_ollama_connection():
    """Test Ollama connection"""
    client = OllamaClient()
    
    print("Testing Ollama connection...")
    if client.is_available():
        print(f"✅ Ollama is available at {OLLAMA_BASE_URL}")
        print(f"   Model: {OLLAMA_MODEL}")
        return True
    else:
        print(f"❌ Ollama not available at {OLLAMA_BASE_URL}")
        print("   Make sure Ollama is running: ollama serve")
        return False


if __name__ == "__main__":
    test_ollama_connection()
