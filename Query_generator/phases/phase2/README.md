# Phase 2: Query Optimizer & Vector Database
**Status:** ✅ **COMPLETED - FROZEN ARTIFACTS**  
**Date:** 2026-02-06

---

## Technical Summary

Phase 2 built a **vector database** and **query optimizer** from the Phase 1 Knowledge Base. This pipeline was executed **once** to generate static artifacts now consumed by Phase 3.

### Executed Pipeline (build-time)

```
Phase 1 KB (18,413 unique organisms)
         ↓
  Query Builder → Generate templated queries
         ↓
  Query Cache → Store generated queries
         ↓
  Sentence Transformer → Vector embeddings (384-dim)
         ↓
  FAISS Index → Similarity search index
```

### Purpose

1. **Query generation:** Extract organisms, strategies, and diseases from KB to create search queries
2. **Vectorization:** Convert queries to embeddings using `sentence-transformers/all-MiniLM-L6-v2`
3. **FAISS indexing:** Build IVF index for fast semantic search (<10ms)
4. **Vocabulary cache:** `query_cache.json` used by Phase 3 for term expansion

---

## Preserved Artifacts

| File | Purpose | Used by Phase 3 |
|------|---------|---------------|
| `data/query_cache.json` | Generated queries with KB-extracted terms | ✅ Yes (vocabulary) |
| `data/pyner_vectors.faiss` | FAISS binary index (IVF, 384-dim) | ❌ No (available for future retrieval) |
| `data/pyner_vectors_index.pkl` | Index metadata (pickle) | ❌ No |

### Cleaned Files

- `logs/*.log` - Temporary execution logs
- Raw vectors and intermediate metadata

---

## Technical Components

### `config.py`
Pipeline configuration:
- Embedding model: `all-MiniLM-L6-v2` (384-dim)
- FAISS index: IVF (Inverted File)
- TOP_K: 10 results by default
- MIN_SIMILARITY: 0.5 threshold

### `scripts/query_builder.py`
Templated query generator:
- Extracts unique organisms from Phase 1 KB
- Generates queries by strategy (RNA-Seq, WGS, 16S)
- Generates queries by disease (COVID-19, Malaria, TB)
- Stores in `query_cache.json`

### `scripts/vector_db.py`
Vector database management:
- `VectorDatabase`: FAISS index construction and loading
- `Retriever`: Semantic search via cosine similarity
- Batch processing for efficiency

### `scripts/main.py`
Pipeline orchestrator:
- Validates Phase 1 KB
- Executes query builder
- Generates embeddings
- Builds FAISS index
- Saves artifacts

---

## Reconstruction (Optional)

If you need to regenerate the FAISS index or query cache:

### Dependencies

```bash
pip install faiss-cpu sentence-transformers numpy
# For GPU (faster):
pip install faiss-gpu sentence-transformers numpy
```

### Execution

```bash
cd Query_generator/phases/phase2
python scripts/main.py
```

**Expected output:**
- `data/pyner_vectors.faiss` - FAISS index
- `data/pyner_vectors_index.pkl` - Metadata
- `data/query_cache.json` - Generated query cache

### Retrieval Testing

```python
from scripts.vector_db import VectorDatabase, Retriever

db = VectorDatabase()
db.load()
retriever = Retriever(db)

# Semantic search
results = retriever.retrieve("COVID-19 gene expression", top_k=5)
for r in results:
    print(f"[{r['similarity_score']:.3f}] {r['query_text']}")
```

---

## Performance Metrics

- **Generated queries:** ~1,200+ diverse terms
- **Embedding speed:** 100-200 queries/sec (CPU), 1000+/sec (GPU)
- **Index size:** ~10-50 MB
- **Search latency:** <10ms per query (FAISS)
- **Accuracy:** Cosine similarity ranking

### Recommended Hardware
- CPU: 2+ cores
- RAM: 4GB minimum (16GB+ recommended)
- GPU: Optional (2-3x speedup with CUDA)

---

## Phase 3 Integration

Phase 3 (**Query Generator API**) consumes:
1. `query_cache.json` → To enrich technical vocabulary
2. FAISS index → **NOT currently used** (available for future retrieval)

Phase 3 query generator uses `technical_vocabulary.json` containing synonyms and terms extracted from this cache.

---

## Final Structure

```
phase2/
├── config.py                      # Pipeline configuration
├── README.md                      # This file
├── data/
│   ├── pyner_vectors.faiss       # FAISS index (binary)
│   ├── pyner_vectors_index.pkl   # Metadata (pickle)
│   └── query_cache.json          # Generated queries (used by Phase 3)
├── logs/                          # Empty logs (cleaned)
└── scripts/
    ├── main.py                    # Pipeline orchestrator
    ├── query_builder.py           # Query generator
    └── vector_db.py               # FAISS and embeddings management
```

---

## Troubleshooting

### FAISS import error
```bash
pip install faiss-cpu  # or faiss-gpu for GPU
```

### Memory issues
- Reduce batch size in `config.py`
- Use CPU mode instead of GPU
- Process queries in smaller batches

### Slow embeddings
- Enable GPU: `USE_GPU = True` in config
- Use smaller model (already using most efficient: `all-MiniLM-L6-v2`)

---

**Status:** ✅ Pipeline executed, artifacts frozen  
**Last Updated:** 2026-02-06  
**Next:** Phase 3 consumes `query_cache.json` for vocabulary expansion
