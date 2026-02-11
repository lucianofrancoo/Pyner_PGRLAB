# Pyner_PGRLAB — Project Overview

Pyner is a local-first NCBI SRA query generator. The current deliverable is a
fully cleaned, production-ready **Query Generator** with a 3-phase pipeline.

---

## What Was Completed (Query Generator)

- **Phase 1:** Knowledge Base extraction from 7.2M+ NCBI SRA XML files.
- **Phase 2:** Vector database and query cache (FAISS + curated vocabulary).
- **Phase 3:** FastAPI + CLI query generator with synonym expansion and
  optional Ollama/Qwen LLM support.
- **Pre-flight checks:** Validates KB availability, Ollama status, query cache,
  and technical vocabulary before starting.
- **Project cleanup:** Obsolete scripts, duplicated docs, logs, and caches were
  removed; documentation consolidated under `Query_generator/`.

---

## Key Structure

- `Query_generator/` — Main product folder (cleaned and documented).
  - `Query_generator/README.md` — Main entry point for usage.
  - `Query_generator/phases/README.md` — Full 3-phase pipeline guide.
  - `Query_generator/phases/phase1/` — KB artifacts and scripts.
  - `Query_generator/phases/phase2/` — FAISS + query cache artifacts.
  - `Query_generator/phases/phase3/` — FastAPI + CLI Query Generator.
  - `Query_generator/test_api.sh` — API smoke test.

- `scripts_iniciales_beta/` — Archived early scripts (kept for reference).
- `archive_old/` — Archived auxiliary files (kept for reference).
- `docs/` — Legacy docs outside Query Generator (kept as requested).

---

## Next Phases (Not Described Yet)

- **fetcher_ncbi**
- **data analyzer**

---

## Quick Start

```bash
cd Query_generator/phases/phase3
python api/main.py "arabidopsis drought rna-seq"
```

For full documentation, see `Query_generator/README.md` and
`Query_generator/phases/README.md`.
