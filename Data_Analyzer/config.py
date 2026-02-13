"""
Data Analyzer - Configuration
==============================
Configuration for paper analysis with Ollama LLM
"""

from pathlib import Path

# ============================================
# PATHS
# ============================================
DATA_ANALYZER_DIR = Path(__file__).parent
ROOT_DIR = DATA_ANALYZER_DIR.parent
OUTPUT_DIR = DATA_ANALYZER_DIR / "output"
LOGS_DIR = DATA_ANALYZER_DIR / "logs"

# Create directories
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# ============================================
# OLLAMA CONFIGURATION
# ============================================
OLLAMA_BASE_URL = "http://localhost:11434"
OLLAMA_MODEL = "qwen2.5:14b"
OLLAMA_TIMEOUT = 120  # seconds

# ============================================
# ANALYSIS CONFIGURATION
# ============================================
# Relevance threshold (0-10 scale)
RELEVANCE_THRESHOLD = 5

# Maximum abstract length to send to LLM (characters)
MAX_ABSTRACT_LENGTH = 3000

# Use full text from PMC when available (PMCID)
USE_PMC_FULL_TEXT = True

# Maximum full text length to send to LLM
MAX_FULL_TEXT_LENGTH = 15000

#Batch size for processing papers
BATCH_SIZE = 10

# ============================================
# OUTPUT CONFIGURATION
# ============================================
CSV_DELIMITER = ","
MULTIVALUE_SEPARATOR = " ; "  # Separator for multiple organisms/tissues/conditions

# CSV columns order
CSV_COLUMNS = [
    "PMID",
    "PMCID",
    "Title",
    "Relevance_Score",
    "Is_Relevant",
    "Organisms",
    "Tissues",
    "Conditions",
    "Strategies",
    "Year",
    "Journal",
    "DOI",
    "Abstract_Preview"
]

# ============================================
# LOGGING
# ============================================
from datetime import datetime
LOG_FILE = LOGS_DIR / f"analyzer_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

print(f"✅ Data Analyzer Config loaded")
print(f"   Ollama: {OLLAMA_BASE_URL}/{OLLAMA_MODEL}")
print(f"   Output: {OUTPUT_DIR}")
