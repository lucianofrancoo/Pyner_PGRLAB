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
        
        return f"""You are a scientific paper classifier and metadata extractor. Analyze the following paper and extract structured information about experimental design and biological samples.

USER QUERY (for relevance): {user_query}

PAPER TITLE: {title}

{content_source}

YOUR TASK - Extract comprehensive experimental metadata:

1. RELEVANCE SCORE (0-10): How well does this paper match the user query?
2. ORGANISMS: Scientific names (e.g., "Solanum lycopersicum", "Mus musculus", "Homo sapiens", "Escherichia coli")
3. SPECIES: Primary species studied (singular, e.g., "Solanum lycopersicum")
4. STRAIN_VARIETY: Specific cultivar/strain/cell line (e.g., "Micro-Tom cultivar", "C57BL/6", "HEK293T")
5. GENOTYPE: wild-type, knockout, transgenic, CRISPR-edited, etc. If not specified, say "not described"
6. TISSUES_ORGANS: Specific tissues analyzed (e.g., "root", "leaf", "root meristem", "brain hippocampus")
7. SOURCE_TISSUE_ORIGIN: If cultured cells/tissue culture (e.g., "leaf mesophyll", "embryonic kidney")
8. CELL_TYPE: If applicable (e.g., "HEK293T cells", "primary neurons", "A549 lung cancer cells")

9. DEVELOPMENTAL_STAGE: Specific stage when samples were collected (e.g., "3-leaf stage", "flowering", "seedling", "adult", "exponential phase", "embryonic stage")
10. ORGANISM_AGE: Age at sampling (e.g., "14 days post-germination", "8-week-old", "P21 postnatal", "6-8 hours from inoculation")
11. GROWTH_PHASE: For microorganisms (e.g., "exponential phase", "stationary phase", "lag phase")

12. CONDITIONS: Main experimental conditions/treatments (e.g., "drought stress", "heat shock", "nitrogen deficiency")
13. ENVIRONMENTAL_STRESS: Specific stressors applied (e.g., "50% soil water deficit", "heat 42°C", "UV-B exposure")
14. TEMPERATURE_RANGE: Specific temperatures (e.g., "25°C", "37°C", "22-25°C")
15. LIGHT_CONDITIONS: Photoperiod and intensity if applicable (e.g., "16h light/8h dark, 350 μmol·m−2·s−1")
16. GROWTH_MEDIUM: Culture substrate (e.g., "Murashige and Skoog medium", "LB agar", "soil mixture 3:1")
17. SAMPLE_COLLECTION_CONDITIONS: How/when samples were collected (e.g., "collected at noon", "collected at dark phase", "water not available for 24h")

18. MOLECULES_EXTRACTED: What biological molecules were analyzed (e.g., "RNA", "Protein", "DNA", "Metabolites", "Lipids")
   - Can be semicolon-separated: "RNA ; Protein ; Metabolites"
19. RNA_TYPE: If RNA analyzed, what type (e.g., "total RNA", "mRNA", "small RNA", "miRNA")
20. DNA_TYPE: If DNA analyzed, what type (e.g., "genomic DNA", "mtDNA", "plasmid DNA")
21. PROTEIN_TYPE: If protein analyzed, what type (e.g., "total protein", "soluble protein", "membrane protein")
22. OTHER_MOLECULES: Additional molecules (e.g., "metabolites", "lipids", "secondary metabolites", "hormones")

23. STRATEGIES: All techniques/methods used (e.g., "RNA-Seq", "qRT-PCR", "microarray", "Western blot")
24. MEASUREMENT_TOOLS: Specific instruments/platforms (e.g., "Illumina TrueSeq", "qRT-PCR on ABI 7500", "MALDI-TOF-MS")
25. DETECTION_METHOD: Detection approach (e.g., "next-generation sequencing", "fluorescence detection", "mass spectrometry")

26. TIME_COURSE_DESIGN: Is this a time-course study? "Yes" or "No"
27. TIME_POINTS: Number of time points sampled (e.g., "6", "not specified")
28. TIME_INTERVALS: How often samples collected (e.g., "every 24 hours", "every 6 hours", "0, 2, 4, 6, 12, 24 hours")
29. TIME_DURATION: Total experiment duration (e.g., "7 days", "48 hours", "not specified")

30. SAMPLE_SIZE: Total number of replicates/samples (e.g., "3 biological replicates", "n=6 per group")
31. BIOLOGICAL_REPLICATES: Number of biological replicates (e.g., "3", "not described")
32. TECHNICAL_REPLICATES: Number of technical replicates per sample (e.g., "2", "not described")
33. REPLICATION_DESIGN: Type of experimental design (e.g., "factorial", "randomized block", "completely randomized")

34. TREATMENT_GROUPS: Number of treatment/comparison groups (e.g., "3", "4")
35. CONTROL_TYPE: Type of control (e.g., "untreated control", "mock treatment", "wild-type", "vehicle control")
36. DOSE_RANGE: If dose-response study (e.g., "0-100 μM", "0.5-5 mg/kg")

37. QUALITY_METRICS: Quality control measures reported (e.g., "RNA integrity number > 8", "protein purity > 95%")
38. CONTAMINATION_CHECK: Contamination screening if mentioned (e.g., "mycoplasma screening negative", "bacterial culture negative")

39. PATHWAY_FOCUS: Main biological pathway/process studied (e.g., "stress response", "photosynthesis", "apoptosis", "immune response")
40. BIOMARKERS_MEASURED: Specific biomarkers or measurements (e.g., "gene expression levels", "enzyme activity", "phenotypic traits")
41. DISEASE_MODEL: If disease study, what model (e.g., "Alzheimer's disease", "Type 2 Diabetes", "cancer xenograft")

42. NORMALIZATION_METHOD: Data normalization approach if mentioned (e.g., "quantile normalization", "DESeq2 normalization")
43. STATISTICAL_METHOD: Statistical tests used (e.g., "t-test", "ANOVA", "Kolmogorov-Smirnov test")
44. DIFFERENTIAL_EXPRESSION_THRESHOLD: Cutoff criteria if applicable (e.g., "log2FC > 1.5, p < 0.05")

45. RAW_DATA_AVAILABLE: Are raw data available? "Yes", "No", or "not mentioned"

RESPONSE FORMAT (JSON only, no explanation):
{{
  "relevance_score": <0-10>,
  "organisms": ["organism1", "organism2"],
  "species": "primary species",
  "strain_variety": "strain or cultivar",
  "genotype": "wild-type or transgenic, etc",
  "tissues_organs": ["tissue1", "tissue2"],
  "source_tissue_origin": "origin if cultured",
  "cell_type": "cell type if applicable",
  "developmental_stage": "stage at sampling",
  "organism_age": "age at sampling",
  "growth_phase": "phase for microorganisms",
  "conditions": ["condition1", "condition2"],
  "environmental_stress": "specific stress applied",
  "temperature_range": "temperature value",
  "light_conditions": "photoperiod and intensity",
  "growth_medium": "substrate or medium",
  "sample_collection_conditions": "how samples were collected",
  "molecules_extracted": ["RNA", "Protein"],
  "rna_type": "total RNA or mRNA",
  "dna_type": "genomic DNA or other",
  "protein_type": "total protein or fraction",
  "other_molecules": "metabolites or other",
  "strategies": ["RNA-Seq", "qRT-PCR"],
  "measurement_tools": "specific platforms",
  "detection_method": "sequencing or detection type",
  "time_course_design": "Yes or No",
  "time_points": "number or not specified",
  "time_intervals": "interval description",
  "time_duration": "total duration",
  "sample_size": "number of replicates",
  "biological_replicates": "number",
  "technical_replicates": "number",
  "replication_design": "design type",
  "treatment_groups": "number of groups",
  "control_type": "type of control",
  "dose_range": "range if dose-response",
  "quality_metrics": "QC measures",
  "contamination_check": "contamination screening",
  "pathway_focus": "biological pathway",
  "biomarkers_measured": "specific biomarkers",
  "disease_model": "disease model if applicable",
  "normalization_method": "normalization approach",
  "statistical_method": "statistical tests",
  "differential_expression_threshold": "cutoff criteria",
  "raw_data_available": "Yes or No"
}}

IMPORTANT RULES:
- YOU MUST RETURN ONLY VALID JSON - NO EXPLANATIONS, NO TEXT BEFORE OR AFTER THE JSON
- Start your response with the opening brace
- End your response with the closing brace
- If information NOT FOUND in paper, use "not described" or empty array []
- NEVER INVENT data. If you don't find it, write "not described"
- For arrays: use semicolon separator internally if multiple values found
- Use exact quotes from Methods/Results when possible
- Be precise with units (°C, μmol·m−2·s−1, hours, etc)

CRITICAL: Your entire response must be ONLY the JSON object. Do not add any text before or after the JSON.

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
                "top_p": 0.9
                # num_predict removed - no token limit, let LLM generate complete response
            }
        }
        
        logger.info(f"Calling Ollama API: {url}")
        response = requests.post(url, json=payload, timeout=timeout)
        response.raise_for_status()
        
        result = response.json()
        
        # Log token usage statistics
        if 'prompt_eval_count' in result:
            prompt_tokens = result.get('prompt_eval_count', 0)
            response_tokens = result.get('eval_count', 0)
            total_tokens = prompt_tokens + response_tokens
            eval_duration_ms = result.get('eval_duration', 0) / 1_000_000  # Convert ns to ms
            
            logger.info(f"   📊 Token usage: Prompt={prompt_tokens}, Response={response_tokens}, Total={total_tokens}")
            logger.info(f"   ⏱️  Generation time: {eval_duration_ms:.0f}ms ({response_tokens / (eval_duration_ms / 1000):.1f} tokens/sec)")
        
        return result.get('response', '')
    
    def _parse_response(self, response: str) -> Dict:
        """Parse Ollama JSON response with all experimental metadata fields"""
        try:
            # Try to find JSON in response
            start = response.find('{')
            end = response.rfind('}') + 1
            
            if start == -1 or end == 0:
                logger.warning("No JSON found in response")
                logger.debug(f"Response preview (first 500 chars): {response[:500]}")
                logger.debug(f"Response preview (last 500 chars): {response[-500:]}")
                return self._empty_result()
            
            json_str = response[start:end]
            parsed = json.loads(json_str)
            
            # Validate and sanitize - all fields with defaults for expandedMetadata extraction
            result = {
                # Basic fields (original)
                'relevance_score': int(parsed.get('relevance_score', 0)),
                'organisms': parsed.get('organisms', []),
                
                # Organism Details (NEW)
                'species': parsed.get('species', 'not described'),
                'strain_variety': parsed.get('strain_variety', 'not described'),
                'genotype': parsed.get('genotype', 'not described'),
                'tissues_organs': parsed.get('tissues_organs', []),
                'source_tissue_origin': parsed.get('source_tissue_origin', 'not described'),
                'cell_type': parsed.get('cell_type', 'not described'),
                
                # Developmental Biology (NEW)
                'developmental_stage': parsed.get('developmental_stage', 'not described'),
                'organism_age': parsed.get('organism_age', 'not described'),
                'growth_phase': parsed.get('growth_phase', 'not described'),
                
                # Sample Collection (NEW)
                'conditions': parsed.get('conditions', []),
                'environmental_stress': parsed.get('environmental_stress', 'not described'),
                'temperature_range': parsed.get('temperature_range', 'not described'),
                'light_conditions': parsed.get('light_conditions', 'not described'),
                'growth_medium': parsed.get('growth_medium', 'not described'),
                'sample_collection_conditions': parsed.get('sample_collection_conditions', 'not described'),
                
                # Molecular Analysis (NEW)
                'molecules_extracted': parsed.get('molecules_extracted', []),
                'rna_type': parsed.get('rna_type', 'not described'),
                'dna_type': parsed.get('dna_type', 'not described'),
                'protein_type': parsed.get('protein_type', 'not described'),
                'other_molecules': parsed.get('other_molecules', 'not described'),
                
                # Experimental Techniques (NEW/ENHANCED)
                'strategies': parsed.get('strategies', []),
                'measurement_tools': parsed.get('measurement_tools', 'not described'),
                'detection_method': parsed.get('detection_method', 'not described'),
                
                # Time Course (NEW)
                'time_course_design': parsed.get('time_course_design', 'not described'),
                'time_points': parsed.get('time_points', 'not described'),
                'time_intervals': parsed.get('time_intervals', 'not described'),
                'time_duration': parsed.get('time_duration', 'not described'),
                
                # Replication (NEW)
                'sample_size': parsed.get('sample_size', 'not described'),
                'biological_replicates': parsed.get('biological_replicates', 'not described'),
                'technical_replicates': parsed.get('technical_replicates', 'not described'),
                'replication_design': parsed.get('replication_design', 'not described'),
                
                # Treatment/Comparisons (NEW)
                'treatment_groups': parsed.get('treatment_groups', 'not described'),
                'control_type': parsed.get('control_type', 'not described'),
                'dose_range': parsed.get('dose_range', 'not described'),
                
                # Quality (NEW)
                'quality_metrics': parsed.get('quality_metrics', 'not described'),
                'contamination_check': parsed.get('contamination_check', 'not described'),
                
                # Biological Context (NEW)
                'pathway_focus': parsed.get('pathway_focus', 'not described'),
                'biomarkers_measured': parsed.get('biomarkers_measured', 'not described'),
                'disease_model': parsed.get('disease_model', 'not described'),
                
                # Data Analysis (NEW)
                'normalization_method': parsed.get('normalization_method', 'not described'),
                'statistical_method': parsed.get('statistical_method', 'not described'),
                'differential_expression_threshold': parsed.get('differential_expression_threshold', 'not described'),
                
                # Data Availability (NEW)
                'raw_data_available': parsed.get('raw_data_available', 'not described')
            }
            
            # Ensure all array fields are lists
            array_fields = ['organisms', 'tissues_organs', 'molecules_extracted', 'strategies']
            for key in array_fields:
                if not isinstance(result[key], list):
                    result[key] = []
            
            # Clean up legacy field name
            if 'tissues' in parsed:
                result['tissues_organs'] = parsed.get('tissues', result.get('tissues_organs', []))
            
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
        """Return empty result structure with all metadata fields"""
        return {
            # Relevance
            'relevance_score': 0,
            
            # Organisms
            'organisms': [],
            'species': 'not described',
            'strain_variety': 'not described',
            'genotype': 'not described',
            
            # Tissues
            'tissues_organs': [],
            'source_tissue_origin': 'not described',
            'cell_type': 'not described',
            
            # Developmental Biology
            'developmental_stage': 'not described',
            'organism_age': 'not described',
            'growth_phase': 'not described',
            
            # Sample Collection
            'conditions': [],
            'environmental_stress': 'not described',
            'temperature_range': 'not described',
            'light_conditions': 'not described',
            'growth_medium': 'not described',
            'sample_collection_conditions': 'not described',
            
            # Molecular Analysis
            'molecules_extracted': [],
            'rna_type': 'not described',
            'dna_type': 'not described',
            'protein_type': 'not described',
            'other_molecules': 'not described',
            
            # Techniques
            'strategies': [],
            'measurement_tools': 'not described',
            'detection_method': 'not described',
            
            # Time Course
            'time_course_design': 'not described',
            'time_points': 'not described',
            'time_intervals': 'not described',
            'time_duration': 'not described',
            
            # Replication
            'sample_size': 'not described',
            'biological_replicates': 'not described',
            'technical_replicates': 'not described',
            'replication_design': 'not described',
            
            # Treatment
            'treatment_groups': 'not described',
            'control_type': 'not described',
            'dose_range': 'not described',
            
            # Quality
            'quality_metrics': 'not described',
            'contamination_check': 'not described',
            
            # Biological Context
            'pathway_focus': 'not described',
            'biomarkers_measured': 'not described',
            'disease_model': 'not described',
            
            # Data Analysis
            'normalization_method': 'not described',
            'statistical_method': 'not described',
            'differential_expression_threshold': 'not described',
            
            # Data Availability
            'raw_data_available': 'not described'
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
