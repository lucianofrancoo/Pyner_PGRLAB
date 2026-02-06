# Phase 2: Query Optimizer
**Date:** 2026-02-06  
**Status:** 🚀 **READY FOR EXECUTION**

---

## Overview

Phase 2 takes the Knowledge Base from Phase 1 (500K BioProjects) and builds a Query Optimizer system:

1. **Query Generation** - Creates search queries from organisms, strategies, diseases
2. **Vector Embeddings** - Converts queries to embeddings using sentence-transformers
3. **FAISS Index** - Builds fast similarity search database
4. **Retrieval System** - Enables semantic search across all queries

---

## Architecture

```
Phase 1 KB (stage3_knowledge_base.json)
         ↓
    Query Builder (extract organisms, strategies, diseases)
         ↓
    Query Cache (save generated queries)
         ↓
    Vector Embedder (sentence-transformers)
         ↓
    FAISS Index (fast similarity search)
         ↓
    Retriever (semantic search engine)
```

---

## Components

### `config.py`
Centralized configuration:
- Embedding model: `sentence-transformers/all-MiniLM-L6-v2`
- Vector dimension: 384
- FAISS index type: IVF (Inverted File)
- Ollama LLM: llama2 (localhost:11434)

### `scripts/query_builder.py`
Generates queries from KB:
- Organism queries (from 18,413 unique organisms)
- Disease queries (COVID-19, Malaria, TB, etc)
- Strategy queries (RNA-Seq, WGS, 16S, etc)
- Template-based generation

### `scripts/vector_db.py`
Vector database management:
- SentenceTransformer for embeddings
- FAISS for indexing
- Metadata caching
- Batch processing

### `scripts/main.py`
Pipeline orchestrator:
- Validates Phase 1 KB
- Generates queries
- Builds vector database
- Initializes retriever
- Tests retrieval accuracy

---

## Quick Start

### 1. Install Dependencies

```bash
pip install faiss-cpu sentence-transformers numpy
# or for GPU support:
pip install faiss-gpu sentence-transformers numpy
```

### 2. Run Phase 2 Pipeline

```bash
cd /home/lahumada/disco1/Pyner_PGRLAB/phase2
python scripts/main.py
```

**Expected Output:**
```
================================================================================
🚀 PHASE 2: QUERY OPTIMIZER - INITIALIZING
================================================================================

✓ Validating Phase 1 Knowledge Base...
✅ KB validated
   Organisms: 18,413
   Experiments: 365,771,189
   Strategies: 37

📝 Generating queries from KB...
✅ Generated 1,234 queries
   - organism: 450
   - disease: 15
   - strategy: 7
   - comparative: 100

🔨 Building vector database...
📥 Loading embedding model: sentence-transformers/all-MiniLM-L6-v2
✅ Embedder loaded (dims: 384)
📝 Embedding 1,234 queries...
🔨 Creating FAISS index...
✅ Index created with 1,234 vectors

🔎 Initializing retriever...
✅ Retriever initialized

================================================================================
✅ PHASE 2 PIPELINE COMPLETE
================================================================================
```

### 3. Test Retrieval

```python
# Example usage
from scripts.vector_db import VectorDatabase, Retriever

db = VectorDatabase()
db.load()  # Load saved index
retriever = Retriever(db)

# Search
results = retriever.retrieve("COVID-19 gene expression", top_k=5)
for r in results:
    print(f"[{r['similarity_score']:.3f}] {r['query_text']}")
```

---

## Output Files

| File | Purpose |
|------|---------|
| `data/pyner_vectors.faiss` | FAISS index (binary) |
| `data/pyner_vectors_index.pkl` | Metadata cache |
| `data/query_cache.json` | Generated queries |
| `logs/phase2_query_optimizer.log` | Execution log |

---

## Performance Metrics

### Expected Results:
- **Queries Generated:** ~1,200+ diverse search queries
- **Embedding Speed:** 100-200 queries/sec (CPU), 1000+/sec (GPU)
- **Index Size:** ~10-50 MB
- **Search Latency:** <10ms per query (FAISS)
- **Accuracy:** Semantic similarity-based ranking

### Hardware Recommendations:
- CPU: 2+ cores minimum
- RAM: 4GB minimum (16GB+ recommended)
- GPU: Optional (2-3x speedup with CUDA)

---

## Configuration

Edit `config.py` to customize:

```python
# Embedding model
EMBEDDING_MODEL = "sentence-transformers/all-MiniLM-L6-v2"
EMBEDDING_DIMENSION = 384

# Search
TOP_K_RESULTS = 10
MIN_SIMILARITY_SCORE = 0.5

# Ollama LLM (for Phase 2 Advanced)
OLLAMA_HOST = "http://localhost:11434"
OLLAMA_MODEL = "llama2"
```

---

## Monitoring

### View Logs:
```bash
tail -f logs/phase2_query_optimizer.log
```

### Check Database Stats:
```python
from scripts.vector_db import VectorDatabase
db = VectorDatabase()
db.load()
print(db.get_stats())
```

### Query Statistics:
```bash
python -c "import json; d=json.load(open('data/query_cache.json')); print(f'Queries: {len(d[\"queries\"])}')"
```

---

## Advanced Features (Phase 2+)

### Ollama LLM Integration
```python
import requests

response = requests.post(
    "http://localhost:11434/api/generate",
    json={"model": "llama2", "prompt": "What are common COVID-19 genes?"}
)
```

### Semantic Similarity Scoring
- Cosine similarity (default)
- Euclidean distance
- Custom scoring functions

### Batch Retrieval:
```python
queries = ["query1", "query2", "query3"]
results = retriever.batch_retrieve(queries, top_k=5)
```

---

## Troubleshooting

### Import Errors
```bash
# If FAISS import fails:
pip install faiss-cpu  # or faiss-gpu

# If sentence-transformers fails:
pip install sentence-transformers
```

### Memory Issues
- Reduce batch size in config.py
- Use CPU mode instead of GPU
- Process queries in smaller batches

### Slow Embedding Generation
- Enable GPU: `USE_GPU = True`
- Use smaller model: `EMBEDDING_MODEL = "all-MiniLM-L6-v2"`
- Enable batch processing

---

## Next Steps: Phase 3

Once Phase 2 is complete:
1. Deploy retriever as API service
2. Integrate with Ollama LLM
3. Build Web UI for search
4. Add real-time query support
5. Scale to production

---

## Files Structure

```
phase2/
├── config.py                   (Configuration)
├── data/                       (Output data)
│   ├── pyner_vectors.faiss    (FAISS index)
│   ├── pyner_vectors_index.pkl (Metadata)
│   ├── query_cache.json       (Generated queries)
│   └── search_results.json    (Search results)
├── logs/                       (Execution logs)
├── scripts/
│   ├── main.py                (orchestrator)
│   ├── query_builder.py       (query generation)
│   └── vector_db.py           (FAISS management)
└── README.md                  (This file)
```

---

## Support

For issues or questions:
1. Check logs in `logs/phase2_query_optimizer.log`
2. Review [DEVELOPMENT_PLAN.md](../planning/DEVELOPMENT_PLAN.md)
3. Check [KB_SCHEMA.md](../docs/KB_SCHEMA.md)

---

**Status:** ✅ Ready for Execution  
**Last Updated:** 2026-02-06  
**Next Phase:** Phase 3 (Production Deployment)
