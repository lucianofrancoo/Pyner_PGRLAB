# Pyner Query Generator - 3-Phase Pipeline
**Status:** ✅ **PRODUCTION READY**  
**Last Updated:** 2026-02-11

---

## Overview

Pyner Query Generator is a **3-phase pipeline** that transforms natural language scientific queries into optimized NCBI SRA boolean queries. The system processes 7.2M+ XML files, builds a vector database, and generates accurate queries with organism/strategy synonym expansion.

```
┌─────────────┐      ┌─────────────┐      ┌─────────────┐
│   Phase 1   │─────▶│   Phase 2   │─────▶│   Phase 3   │
│  Knowledge  │      │   Vector    │      │    Query    │
│    Base     │      │  Database   │      │  Generator  │
└─────────────┘      └─────────────┘      └─────────────┘
     (Frozen)            (Frozen)           (Production)
```

---

## Phase 1: Knowledge Base Extraction

**Purpose:** Extract structured metadata from NCBI SRA XML files to build a comprehensive Knowledge Base.

**Input:** 7.2M+ NCBI SRA XML files  
**Output:** JSON Knowledge Base with organism/strategy/tissue/condition mappings

### Key Features
- **GPU-accelerated processing:** 3x NVIDIA RTX 4000 Ada (24GB VRAM each)
- **Parallel execution:** 20 workers, 100-150 files/sec throughput
- **Checkpoint recovery:** Fault-tolerant processing with automatic resume
- **Metadata extraction:** Organism, strategy, tissue, conditions, study details

### Artifacts Generated
| File | Size | Description |
|------|------|-------------|
| `stage3_kb_reduced.json` | 16MB | Reduced KB with 18,413 unique organisms |
| `stage3_knowledge_base.json` | 11KB | Compact KB summary |
| `stage2_knowledge_base.json` | 10KB | Intermediate KB |
| `stage1_indices.json` | 122KB | XML file indices |

### Execution Time
~3-4 hours (full pipeline with GPU acceleration)

**Status:** ✅ Completed, artifacts frozen

---

## Phase 2: Query Optimizer & Vector Database

**Purpose:** Build a vector database and query cache from Phase 1 Knowledge Base for semantic search and vocabulary enrichment.

**Input:** Phase 1 Knowledge Base (`stage3_knowledge_base.json`)  
**Output:** FAISS vector index + query cache

### Key Features
- **Embedding generation:** `sentence-transformers/all-MiniLM-L6-v2` (384-dim)
- **FAISS indexing:** IVF index for fast similarity search (<10ms)
- **Query templates:** Extract organisms, strategies, diseases for search variations
- **Vocabulary cache:** Generate `query_cache.json` for Phase 3 term expansion

### Artifacts Generated
| File | Size | Description | Used By Phase 3 |
|------|------|-------------|----------------|
| `query_cache.json` | varies | Generated queries with KB-extracted terms | ✅ Yes |
| `pyner_vectors.faiss` | 10-50MB | FAISS binary index | ❌ No (future) |
| `pyner_vectors_index.pkl` | varies | Index metadata | ❌ No |

### Performance
- **Queries generated:** ~1,200+ diverse terms
- **Embedding speed:** 100-200 queries/sec (CPU), 1000+/sec (GPU)
- **Search latency:** <10ms per query

**Status:** ✅ Completed, artifacts frozen

---

## Phase 3: Query Generator API

**Purpose:** Production REST API that converts natural language queries into optimized NCBI SRA boolean queries with synonym expansion.

**Input:** User natural language query (e.g., "arabidopsis drought stress RNA-Seq")  
**Output:** NCBI SRA boolean query with organism/strategy synonym expansion

### Key Features
- **Natural language processing:** Parse user queries to extract entities
- **Knowledge Base validation:** Verify organisms/strategies against Phase 1 KB
- **Synonym expansion:** Add common aliases (e.g., "arabidopsis" → "a. thaliana", "thale cress")
- **NCBI query generation:** Build optimized boolean queries with [Organism] and [All Fields] tags
- **FastAPI server:** Production REST API with health checks and statistics
- **Interactive CLI:** User-friendly command-line interface with corrections

### Query Generation Examples

**Input:**
```bash
python main.py "arabidopsis drought rna-seq"
```

**Output:**
```
("Arabidopsis thaliana"[Organism] OR "arabidopsis"[All Fields] OR 
 "a. thaliana"[All Fields] OR "thale cress"[All Fields]) AND 
("RNA-Seq"[Strategy] OR "rna seq"[All Fields] OR 
 "rna-sequencing"[All Fields] OR "transcriptome sequencing"[All Fields]) AND 
("drought"[All Fields] OR "water stress"[All Fields])
```

### Endpoints

**REST API Mode:**
```bash
# Start server
python main.py --server

# Use API
POST /generate
{
  "query": "arabidopsis drought rna-seq",
  "use_llm": true
}

GET /stats
GET /
```

**CLI Mode:**
```bash
# Quick query
python main.py "arabidopsis drought rna-seq"

# Interactive mode with learning
python main.py -i "arabidopsis drought rna-seq"
```

### Technical Components
- **query_generator.py:** Core query generation logic with KnowledgeBaseValidator and NCBIQueryBuilder
- **interactive_cli.py:** Interactive CLI with user corrections and statistics
- **ollama_integration.py:** Optional Ollama LLM integration for query expansion
- **technical_vocabulary.json:** Organism aliases and strategy keyword mappings

**Status:** ✅ Production ready, actively maintained

---

## Data Flow

```
┌──────────────────────────────────────────────────────────────────┐
│                        Phase 1: KB Extraction                    │
│  NCBI SRA XML (7.2M files) → stage3_kb_reduced.json (16MB)      │
└──────────────────┬───────────────────────────────────────────────┘
                   │
                   ▼
┌──────────────────────────────────────────────────────────────────┐
│                    Phase 2: Vector Database                      │
│  stage3_knowledge_base.json → query_cache.json + FAISS index    │
└──────────────────┬───────────────────────────────────────────────┘
                   │
                   ▼
┌──────────────────────────────────────────────────────────────────┐
│                    Phase 3: Query Generator                      │
│  User Query + KB + query_cache → NCBI Boolean Query             │
└──────────────────────────────────────────────────────────────────┘
```

### Artifact Dependencies

- **Phase 3** consumes:
  - `phase1/output/stage3_kb_reduced.json` - Primary KB for validation
  - `phase1/output/stage2_knowledge_base.json` - Fallback KB
  - `phase2/data/query_cache.json` - Vocabulary enrichment
  - `phase3/support_dictionary/technical_vocabulary.json` - Synonym mappings

- **Phase 2** consumes:
  - `phase1/output/stage3_knowledge_base.json` - Source KB for query generation

---

## Quick Start

### Prerequisites

```bash
# Phase 1 (optional - artifacts already frozen)
pip install torch transformers sentencepiece

# Phase 2 (optional - artifacts already frozen)
pip install faiss-cpu sentence-transformers numpy

# Phase 3 (required for production use)
pip install fastapi uvicorn pydantic
```

### Generate a Query (CLI)

```bash
cd Query_generator/phases/phase3

# Quick mode
python api/main.py "arabidopsis drought rna-seq"

# Interactive mode with learning
python api/main.py -i "arabidopsis drought rna-seq"
```

### Start API Server

```bash
cd Query_generator/phases/phase3

# Start FastAPI server
python api/main.py --server

# Server runs at http://localhost:8000
# API docs at http://localhost:8000/docs
```

### Rebuild Pipeline (Optional)

```bash
# Phase 1: Extract Knowledge Base (~3-4 hours)
cd phase1
python scripts/stage1_parse_xml.py
python scripts/stage2_parallel_gpu.py
python scripts/stage3_final.py

# Phase 2: Build Vector Database (~10-20 minutes)
cd ../phase2
python scripts/main.py

# Phase 3: Use frozen artifacts (no rebuild needed)
cd ../phase3
python api/main.py "your query here"
```

---

## Project Structure

```
phases/
├── README.md                          # This file
├── phase1/                            # Knowledge Base Extraction
│   ├── config.py                      # Phase 1 configuration
│   ├── utils.py                       # Helper utilities
│   ├── README.md                      # Phase 1 documentation
│   ├── output/
│   │   ├── stage3_kb_reduced.json    # ⭐ Primary KB (16MB)
│   │   ├── stage3_knowledge_base.json # Compact KB (11KB)
│   │   ├── stage2_knowledge_base.json # Intermediate KB (10KB)
│   │   └── stage1_indices.json        # XML indices (122KB)
│   └── scripts/
│       ├── stage1_parse_xml.py        # XML parsing
│       ├── stage2_parallel_gpu.py     # GPU parallel processing
│       └── stage3_final.py            # Final KB build
│
├── phase2/                            # Vector Database
│   ├── config.py                      # Phase 2 configuration
│   ├── README.md                      # Phase 2 documentation
│   ├── data/
│   │   ├── query_cache.json          # ⭐ Query cache (used by Phase 3)
│   │   ├── pyner_vectors.faiss       # FAISS index (available for retrieval)
│   │   └── pyner_vectors_index.pkl   # Index metadata
│   └── scripts/
│       ├── main.py                    # Pipeline orchestrator
│       ├── query_builder.py           # Query generator
│       └── vector_db.py               # FAISS management
│
└── phase3/                            # Query Generator API
    ├── config.py                      # Phase 3 configuration
    ├── README.md                      # Phase 3 documentation
    ├── INTERACTIVE_MODE_GUIDE.md      # Interactive CLI guide
    ├── support_dictionary/
    │   └── technical_vocabulary.json  # ⭐ Synonym mappings
    └── api/
        ├── main.py                    # ⭐ FastAPI server / CLI entry point
        ├── query_generator.py         # Core query generation logic
        ├── interactive_cli.py         # Interactive CLI
        └── ollama_integration.py      # Optional LLM integration
```

---

## Performance Metrics

### Phase 1
- **Processing speed:** 100-150 files/sec
- **Total execution:** ~13-20 hours (7.2M files)
- **GPU utilization:** 80-95% (3x RTX 4000 Ada)
- **Output size:** 16MB (reduced KB)

### Phase 2
- **Embedding speed:** 100-200 queries/sec (CPU)
- **Index build time:** ~10-20 minutes
- **Search latency:** <10ms per query
- **Index size:** ~10-50MB

### Phase 3
- **Query generation:** <100ms (without LLM)
- **Query generation:** <2s (with Ollama LLM)
- **API latency:** <50ms (typical)
- **Throughput:** 100+ queries/sec

---

## Troubleshooting

### Phase 1 Issues
- **OOM errors:** Reduce `BATCH_SIZE` and `NUM_WORKERS` in config.py
- **GPU not found:** Check CUDA installation with `nvidia-smi`
- **Checkpoint recovery:** Use `--resume` flag to continue from last checkpoint

### Phase 2 Issues
- **FAISS import error:** `pip install faiss-cpu` or `faiss-gpu`
- **Memory issues:** Reduce batch size, use CPU mode
- **Slow embeddings:** Enable GPU with `USE_GPU = True`

### Phase 3 Issues
- **KB not found:** Ensure Phase 1 artifacts exist in `phase1/output/`
- **Query cache missing:** Phase 3 works without it (degraded vocabulary)
- **Ollama unavailable:** Query generation works without LLM (synonym expansion only)

---

## Future Enhancements

- [ ] **Real-time NCBI integration:** Direct query submission to NCBI SRA
- [ ] **FAISS semantic search:** Use Phase 2 vector index for related query suggestions
- [ ] **Multi-language support:** Support queries in Spanish, French, etc.
- [ ] **Query history:** Track and learn from user query patterns
- [ ] **Advanced filters:** Add date ranges, publication types, etc.

---

**Status:** ✅ Pipeline complete, Phase 3 in production  
**Maintainer:** PGR Lab  
**Last Updated:** 2026-02-11

For detailed documentation on each phase, see individual README files in each phase directory.
