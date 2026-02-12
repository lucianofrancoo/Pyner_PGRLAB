# Pyner_PGRLAB — Project Overview

Pyner is a local-first NCBI SRA query generator. The current deliverable is a
fully cleaned, production-ready **Query Generator** with a 3-phase pipeline.

---

## What Was Completed

### Query Generator ✅
- **Phase 1:** Knowledge Base extraction from 7.2M+ NCBI SRA XML files.
- **Phase 2:** Vector database and query cache (FAISS + curated vocabulary).
- **Phase 3:** FastAPI + CLI query generator with synonym expansion and
  optional Ollama/Qwen LLM support.
- **Pre-flight checks:** Validates KB availability, Ollama status, query cache,
  and technical vocabulary before starting.
- **Project cleanup:** Obsolete scripts, duplicated docs, logs, and caches were
  removed; documentation consolidated under `Query_generator/`.

### Fetcher_NCBI ✅
- **NCBI SRA fetching:** Query NCBI Sequence Read Archive using boolean queries.
- **Automatic deduplication:** Tracks processed BioProjects to avoid duplicates.
- **Metadata extraction:** Parses XML responses for organism, library strategy, 
  biosample, and more.
- **Rate limiting:** Respects NCBI API guidelines (3 req/sec without key, 10/sec 
  with key).
- **Integration:** Seamless pipeline with Query Generator output.

---

## Key Structure

- `Query_generator/` — Main product folder (cleaned and documented).
  - `Query_generator/README.md` — Main entry point for usage.
  - `Query_generator/phases/README.md` — Full 3-phase pipeline guide.
  - `Query_generator/phases/phase1/` — KB artifacts and scripts.
  - `Query_generator/phases/phase2/` — FAISS + query cache artifacts.
  - `Query_generator/phases/phase3/` — FastAPI + CLI Query Generator.
  - `Query_generator/test_api.sh` — API smoke test.

- `Fetcher_NCBI/` — NCBI SRA data fetcher (NEW).
  - `Fetcher_NCBI/README.md` — Complete usage documentation.
  - `Fetcher_NCBI/main.py` — CLI interface for fetching data.
- `test_fetcher_integrator.sh` — Integration test with Query Generator.

- `scripts_iniciales_beta/` — Archived early scripts (kept for reference).
- `archive_old/` — Archived auxiliary files (kept for reference).
- `docs/` — Legacy docs outside Query Generator (kept as requested).

---

## Complete Workflow

**Query Generator → Fetcher_NCBI**

```bash
# 1. Generate optimized NCBI query from natural language
cd Query_generator/phases/phase3
python api/main.py "arabidopsis drought stress rna-seq"

# 2. Fetch data from NCBI SRA
cd ../../Fetcher_NCBI
python main.py -q "GENERATED_QUERY_HERE" -o results.json

# Or use the integration test (interactive)
cd ../..
bash test_fetcher_integrator.sh
```

**Quick Start (Query Generator only)**

```bash
cd Query_generator/phases/phase3
python api/main.py "arabidopsis drought rna-seq"
```

For full documentation:
- Query Generator: `Query_generator/README.md` and `Query_generator/phases/README.md`
- Fetcher_NCBI: `Fetcher_NCBI/README.md`
