"""
Pyner Phase 2 - Configuration
=============================
Configuración centralizada para Phase 2: Query Optimizer

Componentes:
- Query generation engine
- Vector database (FAISS)
- LLM integration (Ollama)
- Scoring system
"""

import os
from pathlib import Path

# ============================================
# PATHS
# ============================================
PROJECT_ROOT = Path(__file__).parent.parent
PHASE1_OUTPUT = PROJECT_ROOT / "phase1" / "output"
PHASE2_ROOT = PROJECT_ROOT / "phase2"
DATA_DIR = PHASE2_ROOT / "data"
LOGS_DIR = PHASE2_ROOT / "logs"
SCRIPTS_DIR = PHASE2_ROOT / "scripts"

# Ensure directories exist
DATA_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# Input knowledge bases
KB_STAGE3 = PHASE1_OUTPUT / "stage3_knowledge_base.json"

# Output files
VECTOR_DB_PATH = DATA_DIR / "pyner_vectors.faiss"
VECTOR_INDEX_PATH = DATA_DIR / "pyner_vectors_index.pkl"
QUERY_CACHE_PATH = DATA_DIR / "query_cache.json"
RESULTS_PATH = DATA_DIR / "search_results.json"

# ============================================
# VECTOR DATABASE CONFIG
# ============================================
EMBEDDING_MODEL = "sentence-transformers/all-MiniLM-L6-v2"
EMBEDDING_DIMENSION = 384
VECTOR_DB_INDEX_TYPE = "IVF"  # or "Flat" for exact search
VECTOR_DB_NPROBE = 10  # IVF search candidates
FAISS_FACTORY_STRING = f"IVF{max(1, 50)},Flat"  # Auto cluster count

# ============================================
# QUERY GENERATION CONFIG
# ============================================
MAX_QUERIES_PER_ORGANISM = 3
MAX_ORGANISMS_TO_QUERY = 100  # Limit organism queries
QUERY_TEMPLATES = {
    "organism": "What genes are expressed in {organism}?",
    "strategy": "Studies using {strategy} sequencing strategy",
    "disease": "Genes involved in {disease} pathology in {organism}",
    "tissue": "Tissue-specific gene expression in {tissue}",
    "comparative": "Gene expression differences between {organism1} and {organism2}",
}

# ============================================
# OLLAMA LLM CONFIG
# ============================================
OLLAMA_HOST = "http://localhost:11434"
OLLAMA_MODEL = "llama2"  # or "mistral", "neural-chat", etc
OLLAMA_TIMEOUT = 60
OLLAMA_TEMPERATURE = 0.7
OLLAMA_TOP_P = 0.9
OLLAMA_MAX_TOKENS = 512

# ============================================
# SEARCH CONFIG
# ============================================
TOP_K_RESULTS = 10
MIN_SIMILARITY_SCORE = 0.5
SIMILARITY_METRIC = "cosine"  # or "euclidean", "l2"

# ============================================
# BATCH PROCESSING
# ============================================
BATCH_SIZE_QUERIES = 32
BATCH_SIZE_VECTORS = 128

# ============================================
# LOGGING CONFIG
# ============================================
LOG_LEVEL = "INFO"
LOG_FORMAT = "[%(asctime)s] %(levelname)-8s | %(name)s | %(message)s"
LOG_FILE = LOGS_DIR / "phase2_query_optimizer.log"

# ============================================
# PERFORMANCE CONFIG
# ============================================
USE_GPU = True
GPU_MEMORY_FRACTION = 0.8
NUM_WORKERS = 4
CACHE_EMBEDDINGS = True

# ============================================
# VALIDATION
# ============================================
VALIDATE_KB_ON_STARTUP = True
MIN_ORGANISMS_REQUIRED = 100
REQUIRED_KB_FIELDS = ["organisms", "statistics", "strategies"]

print(f"✅ Phase 2 Config loaded")
print(f"   KB Path: {KB_STAGE3}")
print(f"   Vector DB: {VECTOR_DB_PATH}")
print(f"   Embedding Model: {EMBEDDING_MODEL}")
print(f"   Ollama: {OLLAMA_HOST}/{OLLAMA_MODEL}")
