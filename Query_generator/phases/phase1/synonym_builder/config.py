"""
Configuration for Synonym Builder
"""
from pathlib import Path

# Base directories
BASE_DIR = Path(__file__).parent
PHASE1_DIR = BASE_DIR.parent
PUBMED_EXPLORER_DIR = PHASE1_DIR / "pubmed_explorer"
PUBMED_BASELINE_DIR = PHASE1_DIR / "pubmed_baseline"

# MeSH data
MESH_DATA_DIR = PUBMED_EXPLORER_DIR / "mesh_data"
MESH_DESC_FILE = MESH_DATA_DIR / "desc2026.xml"

# Output directories
OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# Output files
MESH_SYNONYMS_OUTPUT = OUTPUT_DIR / "step1_mesh_synonyms.json"
ABSTRACT_VARIANTS_OUTPUT = OUTPUT_DIR / "step2_abstract_variants.json"
FINAL_DICTIONARY_OUTPUT = OUTPUT_DIR / "final_synonym_dictionary.json"
STATS_OUTPUT = OUTPUT_DIR / "stats.json"

# Processing parameters
# Etapa 1: MeSH
MESH_INCLUDE_PLURAL = True  # Incluir formas plurales como sinónimos

# Etapa 2: Abstract Mining
ABSTRACT_MIN_FREQUENCY = 5  # Mínimo número de apariciones
ABSTRACT_MIN_PERCENTAGE = 0.05  # Mínimo 5% de artículos con el MeSH term
MAX_NGRAM_SIZE = 4  # Máximo tamaño de n-gramas para buscar variantes

# Merge & Ranking
FINAL_MIN_CONFIDENCE = 0.6  # Mínima confianza para incluir en diccionario final

# Logging
LOG_LEVEL = "INFO"
LOG_PROGRESS_INTERVAL = 1000  # Log cada N artículos procesados
