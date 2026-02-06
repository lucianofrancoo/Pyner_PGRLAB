# Query_generator — Backend Integration Instructions (for an AI or Frontend)

Purpose: enable another AI agent or a frontend developer to call the backend tools we built and obtain an NCBI-ready boolean query (or use the semantic API). This doc describes how to run the code, the expected inputs and outputs, fallback behavior, and example calls.

1) Requirements
- Python 3.10+ (the repo was developed with CPython 3.12 but 3.10+ should work)
- System packages: `curl`, network access to `localhost` and optionally HF Hub.
- Optional: `ollama` LLM running at `http://localhost:11434/llama2` to enable query expansion.

2) Repository layout (important paths)
- `Query_generator/generate_boolean_query.py` — CLI that produces an NCBI boolean query from a natural language input. Preferred tool for frontend when you need the final boolean string.
- `Query_generator/test_queries.py` — simple CLI that calls the backend `/search` endpoint and prints human-readable results.
- `Query_generator/quick_test.py` — full local quick-test (loads KB + FAISS + Ollama fallback checks).
- `Query_generator/phases/` — copies of `phase1/phase2/phase3` used by the backend (KB, FAISS, API code). Keep these files; they are the authoritative KB and vector index.

3) How the backend works (short)
- Input (natural language) → embeddings → FAISS vector search (top-K hits) → KB enrichment (organisms, strategies, tissues) → optional LLM expansion → boolean query builder.
- The final boolean query builder uses the KB to decide when to tag terms as `[...]` (for example `Mus musculus[Organism]`). If KB does not validate a candidate organism, it leaves the term as an `All Fields` quoted phrase.

4) How to run locally (commands an AI can execute)

- 4.A — If the API is up and running on `localhost:8000` (recommended):

  - Health check (recommended):
    ```bash
    curl -s http://localhost:8000/ | python3 -m json.tool
    ```

  - Generate boolean query by calling the CLI wrapper (no API call needed):
    ```bash
    python3 Query_generator/generate_boolean_query.py "Arabidopsis root drought"
    ```

  - Or call the semantic API directly (returns JSON):
    ```bash
    curl -s -X POST http://localhost:8000/search \
      -H "Content-Type: application/json" \
      -d '{"query":"Arabidopsis root drought","top_k":5,"expand":false}' | python3 -m json.tool
    ```

- 4.B — If API is not running, use the CLI that falls back to local phase copies (this repo contains the fallback):
  ```bash
  python3 Query_generator/generate_boolean_query.py "Mouse RNAseq"
  ```

  The script will attempt to call the API; if it times out it will load the local FAISS index from `Query_generator/phases/phase2/data/pyner_vectors.faiss` and compute the suggestions locally.

5) Input/Output expectations (for the frontend or another AI agent)
- Input: a single natural-language string describing the search intent (e.g. `"Mouse RNAseq"`, `"Arabidopsis root drought"`).
- Output (CLI `generate_boolean_query.py`): prints two sections:
  - Top semantic suggestions (text list with similarity scores).
  - `Generated boolean query (NCBI-ready)`: a single-line boolean string that can be pasted into the NCBI SRA Advanced Search builder or used via Entrez E-utilities.

Example output:
```
(Mus musculus[Organism]) AND ("RNA-Seq" OR "RNA sequencing" OR transcriptome OR rnaseq) AND ("gene expression" OR transcriptome OR "RNA-Seq" OR expression)
```

6) Output format for programmatic consumption
- If you prefer machine-friendly JSON, call the semantic API `/search` to get structured `results` (each result contains `query_text`, `query_type`, `similarity_score`, and `kb_data`). The frontend can then programmatically combine facets.
- Example JSON fields from `/search`:
  - `results`: list of objects with `query_text`, `query_type`, `similarity_score`, `kb_data`.
  - `execution_time`: time in seconds.

7) Integration tips for the frontend developer / AI agent
- Preferred flow for a UI that needs a single boolean to run a query:
  1. Call `/search` (POST) with `expand=false` to receive top-K semantic hits.
  2. Show the top suggestions to the user (optional) and offer toggles (include/exclude strategies, tissues).
  3. When user confirms, call `generate_boolean_query.py` on the backend (or the backend can call the same code path) to build the final boolean string.

- If the frontend wants the backend to return both the boolean and the structured items, expose a small wrapper endpoint that runs the generation function and returns:

  ```json
  {
    "boolean_query": "(Mus musculus[Organism]) AND (...)",
    "suggestions": [ {"query_text":"...","query_type":"...","score":0.52}, ... ]
  }
  ```

8) Environmental variables and optional settings
- `HF_TOKEN` (optional): set it to speed up model downloads if the embedding model needs to fetch weights.
- `OLLAMA_URL` (optional): point to your Ollama instance (default used by code is `http://localhost:11434/llama2`).

9) Security and rate limits
- This backend is local by default. If the frontend is remote, add an authenticated API gateway in front of `localhost:8000` before exposing it publicly.

10) Troubleshooting
- If `generate_boolean_query.py` times out contacting `localhost:8000`, it will attempt a local retrieval from `Query_generator/phases`. Make sure `Query_generator/phases/phase2/data/pyner_vectors.faiss` exists.
- If organism detection is incorrect, check `Query_generator/phases/phase1/output/stage3_knowledge_base.json` — the boolean builder validates organism names against that KB.

11) Quick programmatic example (Python) for another AI agent

```python
import subprocess, json

def get_boolean(query):
    proc = subprocess.run(['python3','Query_generator/generate_boolean_query.py', query], capture_output=True, text=True)
    return proc.stdout

print(get_boolean('Arabidopsis root drought'))
```

12) Notes for frontend design team
- We only provide backend logic and boolean generation. The frontend should:
  - Request the boolean from the backend, display to user, and use it to call NCBI or show in the SRA Advanced Search form.
  - Optionally let users edit the boolean before submission.

Contact: The KB and FAISS index used are under `Query_generator/phases/phase1` and `Query_generator/phases/phase2` respectively.
