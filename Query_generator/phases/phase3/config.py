"""
Pyner Phase 3 - Configuration
=============================
Production deployment configuration
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
PHASE2_DATA = PROJECT_ROOT / "phase2" / "data"
PHASE3_ROOT = PROJECT_ROOT / "phase3"
LOGS_DIR = PHASE3_ROOT / "logs"
SUPPORT_DICT_DIR = PHASE3_ROOT / "support_dictionary"

# API Server
API_HOST = "0.0.0.0"
API_PORT = 8000
API_WORKERS = 4
API_RELOAD = True

# Ollama LLM
OLLAMA_HOST = "http://localhost:11434"
OLLAMA_MODEL = "qwen2.5:14b"
OLLAMA_TIMEOUT = 30
QUERY_EXPANSION_ENABLED = True

# Vector DB paths
VECTOR_DB_PATH = PHASE2_DATA / "pyner_vectors.faiss"
QUERY_CACHE_PATH = PHASE2_DATA / "query_cache.json"

# Search config
TOP_K_RESULTS = 10
MIN_SIMILARITY_SCORE = 0.3

# Logging
LOG_LEVEL = "INFO"
LOG_FILE = LOGS_DIR / "phase3_api.log"

print(f"✅ Phase 3 Config loaded")
print(f"   API: {API_HOST}:{API_PORT}")
print(f"   Ollama: {OLLAMA_HOST}/{OLLAMA_MODEL}")
