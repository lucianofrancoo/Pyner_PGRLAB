# Pyner Query Generator
**Status:** ✅ **PRODUCTION READY**  
**Last Updated:** 2026-02-11

---

## Overview

Pyner Query Generator transforms natural language scientific queries into optimized NCBI SRA boolean queries through a 3-phase pipeline:

1. **Phase 1:** Knowledge Base Extraction (7.2M+ XML files → structured KB)
2. **Phase 2:** Vector Database (KB → FAISS index + query cache)
3. **Phase 3:** Query Generator API (Natural language → NCBI queries)

**Example:**
```bash
Input:  "arabidopsis drought rna-seq"
Output: ("Arabidopsis thaliana"[Organism] OR "arabidopsis"[All Fields] OR 
         "a. thaliana"[All Fields]) AND ("RNA-Seq"[Strategy] OR 
         "rna sequencing"[All Fields]) AND "drought"[All Fields]
```

---

## Quick Start

### Generate a Query (CLI)

```bash
cd phases/phase3
python api/main.py "your query here"

# Interactive mode with learning
python api/main.py -i "arabidopsis drought"
```

### Start API Server

```bash
cd phases/phase3
python api/main.py --server

# Server runs at http://localhost:8000
# API docs at http://localhost:8000/docs
```

### Test API

```bash
# Health check
curl http://localhost:8000/

# Generate query
curl -X POST http://localhost:8000/generate \
  -H "Content-Type: application/json" \
  -d '{"query": "arabidopsis drought rna-seq", "use_llm": false}'

# Or use the test script
./test_api.sh
```

---

## Documentation

📖 **[Complete Documentation](phases/README.md)** - Full 3-phase pipeline guide

### Quick Links
- [Phase 1 README](phases/phase1/README.md) - Knowledge Base Extraction
- [Phase 2 README](phases/phase2/README.md) - Vector Database
- [Phase 3 README](phases/phase3/README.md) - Query Generator API
- [Interactive CLI Guide](phases/phase3/INTERACTIVE_MODE_GUIDE.md) - User-friendly query generation

---

## Project Structure

```
Query_generator/
├── README.md              # This file
├── test_api.sh            # API testing script
└── phases/
    ├── README.md          # 📖 Complete documentation (START HERE)
    │
    ├── phase1/            # Knowledge Base Extraction
    │   ├── README.md
    │   ├── config.py
    │   ├── output/
    │   │   ├── stage3_kb_reduced.json    (16MB - used by Phase 3)
    │   │   └── stage3_knowledge_base.json (11KB - used by Phase 2)
    │   └── scripts/
    │
    ├── phase2/            # Vector Database
    │   ├── README.md
    │   ├── config.py
    │   ├── data/
    │   │   ├── query_cache.json          (used by Phase 3)
    │   │   └── pyner_vectors.faiss       (FAISS index)
    │   └── scripts/
    │
    └── phase3/            # Query Generator API ⭐
        ├── README.md
        ├── INTERACTIVE_MODE_GUIDE.md
        ├── config.py
        ├── api/
        │   ├── main.py                   # FastAPI server / CLI entry point
        │   ├── query_generator.py        # Core query generation logic
        │   ├── interactive_cli.py        # Interactive CLI
        │   └── ollama_integration.py     # Optional LLM integration
        └── support_dictionary/
            └── technical_vocabulary.json # Organism/strategy synonyms
```

---

## Requirements

### Phase 3 (Production - Required)
```bash
pip install fastapi uvicorn pydantic
```

### Phase 1 & 2 (Build - Optional, artifacts already frozen)
```bash
# Phase 1
pip install torch transformers sentencepiece

# Phase 2
pip install faiss-cpu sentence-transformers numpy
```

---

## Examples

### CLI Examples

```bash
# Simple query
python phases/phase3/api/main.py "human cancer rna-seq"

# With organism correction
python phases/phase3/api/main.py "mice liver transcriptome"

# Interactive mode with statistics
python phases/phase3/api/main.py -i "arabidopsis drought"
```

### API Examples

```bash
# Start server
python phases/phase3/api/main.py --server

# Generate query
curl -X POST http://localhost:8000/generate \
  -H "Content-Type: application/json" \
  -d '{
    "query": "arabidopsis thaliana drought stress RNA-Seq",
    "use_llm": false
  }'

# Get statistics
curl http://localhost:8000/stats
```

### Python Integration

```python
from phase3.api.query_generator import QueryGeneratorService
from phase1.config import OUTPUT_DIR

kb_path = OUTPUT_DIR / "stage3_kb_reduced.json"
service = QueryGeneratorService(kb_path)

result = service.generate_query("arabidopsis drought rna-seq")
print(result["ncbi_query"])
```

---

## Performance

- **Query generation:** <100ms (without LLM)
- **API latency:** <50ms (typical)
- **Throughput:** 100+ queries/sec
- **KB size:** 16MB (18,413 organisms)
- **Synonym expansion:** Automatic for organisms and strategies

---

## Troubleshooting

### KB not found
```bash
# Ensure Phase 1 artifacts exist
ls phases/phase1/output/stage3_kb_reduced.json
```

### Module import errors
```bash
# Add to PYTHONPATH
export PYTHONPATH=/home/lahumada/disco1/Pyner_PGRLAB:$PYTHONPATH
```

### Ollama unavailable
Query generation works without LLM (uses synonym expansion only)

---

**For complete documentation, architecture details, and rebuild instructions:**  
👉 **[phases/README.md](phases/README.md)**

**Maintainer:** PGR Lab  
**Repository:** [github.com/lucianofrancoo/Pyner_PGRLAB](https://github.com/lucianofrancoo/Pyner_PGRLAB)
