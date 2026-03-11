from __future__ import annotations

import copy
import logging
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

ROOT_DIR = Path(__file__).resolve().parents[2]
FETCHER_DIR = ROOT_DIR / "Fetcher_NCBI"
PHASES_DIR = ROOT_DIR / "Query_generator" / "phases"
QUERY_GENERATOR_DIR = ROOT_DIR / "Query_generator"
QUERY_GENERATOR_SCRIPT = QUERY_GENERATOR_DIR / "test_query_expander.py"
DATA_ANALYZER_DIR = ROOT_DIR / "Data_Analyzer"
DATA_VISUALIZATION_DIR = ROOT_DIR / "Data_visualization"

sys.path.insert(0, str(FETCHER_DIR))
sys.path.insert(0, str(PHASES_DIR))
# DATA_ANALYZER_DIR se agrega al final para no interferir con Fetcher_NCBI/config.py
sys.path.append(str(DATA_ANALYZER_DIR))
sys.path.append(str(DATA_VISUALIZATION_DIR))

from boolean_fetcher_integrated import BooleanFetcherIntegrated
from ncbi_linkout import LinkoutFetcher
from phase1.config import OUTPUT_DIR
from phase3.api.ollama_integration import OllamaClient
from phase3.api.query_generator import QueryGeneratorService
from phase3.config import SUPPORT_DICT_DIR

from classification import ResultClassifier

logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(levelname)s | %(message)s")
logger = logging.getLogger("minero_web")

LLM_CLASSIFICATION_LIMITS = {
    "pubmed": 8,
    "bioproject": 5,
}
QUERY_GENERATION_TIMEOUT_SEC = 45
QUERY_GENERATION_SCRIPT_TIMEOUT_SEC = 180
FRONTEND_MAX_QUERY_LIST_ITEMS = 80
FRONTEND_MAX_QUERY_TERM_LENGTH = 180


class SearchRequest(BaseModel):
    natural_query: str = Field(min_length=3, description="Consulta biologica en lenguaje natural")
    source: Literal["pubmed", "pmc", "bioproject"] = "pubmed"
    max_results: int = Field(default=20, ge=1, le=10000)
    use_llm: bool = True
    request_id: Optional[str] = None


class GenerateQueryRequest(BaseModel):
    natural_query: str = Field(min_length=3, description="Consulta biologica en lenguaje natural")
    use_llm: bool = True


class RunSearchRequest(BaseModel):
    source: Literal["pubmed", "pmc", "bioproject"] = "pubmed"
    max_results: int = Field(default=20, ge=1, le=10000)
    use_llm: bool = True
    ncbi_query: str = Field(min_length=3, description="Query booleana NCBI generada")
    query_generation: Dict[str, Any] = Field(default_factory=dict)
    request_id: Optional[str] = None


class AppState:
    query_service: Optional[QueryGeneratorService] = None
    classifier: Optional[ResultClassifier] = None
    llm_client: Optional[OllamaClient] = None
    search_progress: Dict[str, Dict[str, Any]] = {}


state = AppState()

app = FastAPI(
    title="Minero Web API",
    version="1.0.0",
    description="API para app web Minero: query -> fetch -> clasificacion -> visualizacion",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[],
    allow_origin_regex=".*",
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
def startup() -> None:
    logger.info("Inicializando Minero Web API...")

    kb_path = OUTPUT_DIR / "stage3_kb_reduced.json"
    query_cache_path = OUTPUT_DIR.parent.parent / "phase2" / "data" / "query_cache.json"

    if not kb_path.exists():
        raise RuntimeError(f"Knowledge Base no encontrada: {kb_path}")

    llm_client = None
    try:
        probe = OllamaClient()
        if probe.is_available():
            llm_client = probe
            logger.info("Ollama disponible para generacion/clasificacion")
        else:
            logger.warning("Ollama no disponible, se usara fallback heuristico")
    except Exception as exc:
        logger.warning("No se pudo inicializar Ollama: %s", exc)

    state.llm_client = llm_client
    state.query_service = QueryGeneratorService(
        kb_path=kb_path,
        ollama_client=llm_client,
        query_cache_path=query_cache_path,
        cache_dir=SUPPORT_DICT_DIR,
    )
    state.classifier = ResultClassifier(ollama_client=llm_client)


@app.get("/api/minero/health")
def health() -> Dict[str, Any]:
    return {
        "status": "ok",
        "service": "minero-web-api",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "llm_runtime_available": bool(state.classifier and state.classifier.llm_available),
    }


@app.get("/api/minero/progress/{request_id}")
def get_progress(request_id: str) -> Dict[str, Any]:
    progress = state.search_progress.get(request_id)
    if not progress:
        raise HTTPException(status_code=404, detail="Progress not found")
    return progress


def _ensure_services() -> tuple[QueryGeneratorService, ResultClassifier]:
    if not state.query_service or not state.classifier:
        raise HTTPException(status_code=503, detail="Servicio no inicializado")
    return state.query_service, state.classifier


def _set_progress(
    request_id: Optional[str],
    *,
    stage: Optional[str] = None,
    processed: Optional[int] = None,
    target: Optional[int] = None,
    message: Optional[str] = None,
    done: Optional[bool] = None,
) -> None:
    if not request_id:
        return
    current = state.search_progress.get(request_id, {})
    payload = {
        **current,
        "request_id": request_id,
        "stage": stage if stage is not None else current.get("stage", "queued"),
        "processed": int(processed if processed is not None else current.get("processed", 0)),
        "target": int(target if target is not None else current.get("target", 0)),
        "message": message if message is not None else current.get("message", ""),
        "done": bool(done if done is not None else current.get("done", False)),
        "updated_at": datetime.now(timezone.utc).isoformat(),
    }
    state.search_progress[request_id] = payload


def _generate_query_payload(natural_query: str, use_llm_requested: bool) -> Dict[str, Any]:
    # Primary path: use exactly the same script path as pyner_miner.sh.
    phase1_result = _build_phase1_boolean_query(natural_query, use_llm_requested)
    if phase1_result and phase1_result.get("ncbi_query"):
        return {
            "user_input": natural_query,
            "extracted": {
                "organism": None,
                "organism_variants": [],
                "strategies": [],
                "tissues": [],
                "conditions": [],
                "genes": [],
                "free_terms": [],
            },
            "synonyms": {
                "organism": [],
                "strategies": [],
                "tissues": [],
                "conditions": [],
                "genes": [],
            },
            "ncbi_query": phase1_result["ncbi_query"],
            "ready_to_use": True,
            "clarification_needed": False,
            "clarification_message": "",
            "query_mode": "phase1_expander",
        }

    # Fallback: legacy Phase 3 query generator
    query_service, classifier = _ensure_services()
    use_llm = use_llm_requested and classifier.llm_available

    try:
        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(query_service.generate_query, natural_query, use_llm)
            query_payload = future.result(timeout=QUERY_GENERATION_TIMEOUT_SEC)
    except FuturesTimeoutError:
        if use_llm:
            logger.warning(
                "Generate-query timeout with LLM after %ss; retrying with heuristic mode",
                QUERY_GENERATION_TIMEOUT_SEC,
            )
            query_payload = query_service.generate_query(natural_query, use_llm=False)
        else:
            raise HTTPException(
                status_code=504,
                detail=f"Query generation timeout after {QUERY_GENERATION_TIMEOUT_SEC}s",
            )
    except Exception as exc:
        if use_llm:
            logger.warning("Generate-query failed with LLM (%s); retrying with heuristic mode", exc)
            query_payload = query_service.generate_query(natural_query, use_llm=False)
        else:
            raise

    if not query_payload.get("ready_to_use"):
        message = query_payload.get("clarification_message") or "La consulta requiere mas contexto"
        raise HTTPException(status_code=400, detail=message)

    query_payload["query_mode"] = "phase3_legacy"
    return query_payload


def _trim_string_list(values: Any, max_items: int, max_term_length: int) -> tuple[List[str], int]:
    if not isinstance(values, list):
        return [], 0

    trimmed: List[str] = []
    for item in values[:max_items]:
        if isinstance(item, str):
            text = item.strip()
            if text:
                trimmed.append(text[:max_term_length])
    omitted = max(0, len(values) - len(trimmed))
    return trimmed, omitted


def _compact_query_payload_for_frontend(query_payload: Dict[str, Any]) -> Dict[str, Any]:
    """
    Reduce payload size before returning to frontend while preserving backend behavior.
    Upstream query generation can include large synonym/variant arrays.
    """
    compacted = copy.deepcopy(query_payload)
    truncation: Dict[str, int] = {}

    extracted = compacted.get("extracted")
    if isinstance(extracted, dict):
        for key in ("organism_variants", "strategies", "tissues", "conditions", "genes", "free_terms"):
            values = extracted.get(key)
            trimmed, omitted = _trim_string_list(
                values,
                max_items=FRONTEND_MAX_QUERY_LIST_ITEMS,
                max_term_length=FRONTEND_MAX_QUERY_TERM_LENGTH,
            )
            if isinstance(values, list):
                extracted[key] = trimmed
                if omitted:
                    truncation[f"extracted.{key}"] = omitted

    synonyms = compacted.get("synonyms")
    if isinstance(synonyms, dict):
        for key in ("organism", "strategies", "tissues", "conditions", "genes"):
            values = synonyms.get(key)
            trimmed, omitted = _trim_string_list(
                values,
                max_items=FRONTEND_MAX_QUERY_LIST_ITEMS,
                max_term_length=FRONTEND_MAX_QUERY_TERM_LENGTH,
            )
            if isinstance(values, list):
                synonyms[key] = trimmed
                if omitted:
                    truncation[f"synonyms.{key}"] = omitted

    if truncation:
        compacted["frontend_truncation"] = truncation

    return compacted


def _build_phase1_boolean_query(user_input: str, _use_llm_requested: bool) -> Optional[Dict[str, Any]]:
    if not QUERY_GENERATOR_SCRIPT.exists():
        logger.warning("Query generator script not found: %s", QUERY_GENERATOR_SCRIPT)
        return None

    try:
        cmd = ["python3", str(QUERY_GENERATOR_SCRIPT), user_input]
        # Intentionally run with the same invocation as pyner_miner.sh:
        # python3 test_query_expander.py "$USER_INPUT"
        logger.info("[phase1-script] Running: %s", " ".join(cmd))
        proc = subprocess.Popen(
            cmd,
            cwd=str(QUERY_GENERATOR_DIR),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        output_lines: List[str] = []
        assert proc.stdout is not None
        for line in proc.stdout:
            # Mirror script output in backend console for debugging (same as CLI flow).
            print(line, end="")
            output_lines.append(line.rstrip("\n"))

        return_code = proc.wait(timeout=QUERY_GENERATION_SCRIPT_TIMEOUT_SEC)
        if return_code != 0:
            logger.warning("phase1 script returned non-zero exit code: %s", return_code)
            return None

        lines = output_lines
        for idx, line in enumerate(lines):
            if line.strip() == "NCBI Query:":
                for query_line in lines[idx + 1 :]:
                    query = query_line.strip()
                    if query:
                        return {"ncbi_query": query}
                break
        logger.warning("phase1 script did not emit a parsable 'NCBI Query:' block")
    except subprocess.TimeoutExpired:
        logger.warning("Phase1 script timed out after %ss", QUERY_GENERATION_SCRIPT_TIMEOUT_SEC)
    except Exception as exc:
        logger.warning("Phase1 script execution failed; using legacy query. Error: %s", exc)
    return None


def _execute_search(
    source: Literal["pubmed", "pmc", "bioproject"],
    max_results: int,
    use_llm_requested: bool,
    ncbi_query: str,
    query_payload: Dict[str, Any],
    request_id: Optional[str] = None,
) -> Dict[str, Any]:
    _, classifier = _ensure_services()
    use_llm = use_llm_requested and classifier.llm_available

    # PMC se trata igual que pubmed para efectos de límites LLM
    llm_source_key = "pubmed" if source in ("pubmed", "pmc") else "bioproject"
    llm_limit = LLM_CLASSIFICATION_LIMITS[llm_source_key]
    use_llm_classification = use_llm and max_results <= llm_limit

    if use_llm and not use_llm_classification:
        logger.warning(
            "LLM classification disabled for this run: source=%s max_results=%s limit=%s (using heuristic)",
            source,
            max_results,
            llm_limit,
        )

    ncbi_query = query_payload.get("ncbi_query", "") or ncbi_query
    if not ncbi_query:
        raise HTTPException(status_code=500, detail="No se pudo generar query NCBI")

    _set_progress(
        request_id,
        stage="searching",
        processed=0,
        target=max_results,
        message="Submitting query to NCBI...",
        done=False,
    )

    raw_results: List[Dict[str, Any]] = []
    results: List[Dict[str, Any]] = []

    if source in ("pubmed", "pmc"):
        # PMC usa la misma interfaz pero busca en full-text (db='pmc')
        db = "pmc" if source == "pmc" else "pubmed"
        fetcher = LinkoutFetcher()
        last_processed = 0

        def on_pubmed_progress(event: Dict[str, Any]) -> None:
            nonlocal last_processed
            stage = str(event.get("stage") or "searching")
            message = str(event.get("message") or "")
            processed = int(event.get("processed") or 0)
            target = int(event.get("target") or max_results)
            last_processed = max(last_processed, processed)
            _set_progress(
                request_id,
                stage=stage,
                processed=processed,
                target=target,
                message=message,
            )

        raw_results = fetcher.search_publications_by_boolean_query(
            ncbi_query,
            max_results=max_results,
            db=db,
            progress_callback=on_pubmed_progress,
        )
        _set_progress(
            request_id,
            stage="classifying",
            processed=last_processed,
            target=max(len(raw_results), max_results),
            message="Classifying records...",
        )
        results = [
            classifier.classify_pubmed(
                record=item,
                query_payload=query_payload,
                use_llm=use_llm_classification,
            )
            for item in raw_results
        ]
    else:
        fetcher = BooleanFetcherIntegrated()
        _set_progress(
            request_id,
            stage="searching",
            processed=0,
            target=max_results,
            message="Resolving BioProject -> SRA -> PubMed...",
        )

        def on_bioproject_progress(event: Dict[str, Any]) -> None:
            _set_progress(
                request_id,
                stage=str(event.get("stage") or "searching"),
                processed=int(event.get("processed") or 0),
                target=int(event.get("target") or max_results),
                message=str(event.get("message") or ""),
            )

        raw_results = fetcher.run_workflow(
            ncbi_query,
            max_bioproject=max_results,
            progress_callback=on_bioproject_progress,
        )
        _set_progress(
            request_id,
            stage="classifying",
            processed=len(raw_results),
            target=max(len(raw_results), max_results),
            message="Classifying projects...",
        )
        results = [
            classifier.classify_bioproject(
                record=item,
                query_payload=query_payload,
                use_llm=use_llm_classification,
            )
            for item in raw_results
        ]

    partial_success = any(bool(item.get("error")) for item in raw_results)

    if not results:
        status = "empty"
    elif partial_success:
        status = "partial-success"
    else:
        status = "success"

    _set_progress(
        request_id,
        stage="done",
        processed=len(results),
        target=max(len(results), max_results),
        message=f"Completed with {len(results)} records.",
        done=True,
    )

    return {
        "metadata": {
            "status": status,
            "query": ncbi_query,
            "source": source,
            "total_results": len(results),
            "classification_version": classifier.version,
            "classification_timestamp": datetime.now(timezone.utc).isoformat(),
            "llm_runtime_available": classifier.llm_available,
            "model_default": "ollama" if use_llm_classification else "heuristic",
        },
        "query_generation": query_payload,
        "results": results,
    }


@app.post("/api/minero/generate-query")
def generate_query(request: GenerateQueryRequest) -> Dict[str, Any]:
    query_payload = _generate_query_payload(request.natural_query, request.use_llm)
    return {"query_generation": _compact_query_payload_for_frontend(query_payload)}


@app.post("/api/minero/run-search")
def run_search(request: RunSearchRequest) -> Dict[str, Any]:
    _set_progress(
        request.request_id,
        stage="queued",
        processed=0,
        target=request.max_results,
        message="Search queued...",
        done=False,
    )
    query_payload = request.query_generation or {}
    query_payload["ncbi_query"] = request.ncbi_query
    # Si la fuente es PMC, normalizar field tags para PMC antes de buscar
    ncbi_query = request.ncbi_query
    if request.source == "pmc":
        import re as _re
        ncbi_query = _re.sub(r"\[Organism\]", "[all]", ncbi_query)
        ncbi_query = _re.sub(r"\[All Fields\]", "[all]", ncbi_query)
        query_payload["ncbi_query"] = ncbi_query
    response = _execute_search(
        source=request.source,
        max_results=request.max_results,
        use_llm_requested=request.use_llm,
        ncbi_query=ncbi_query,
        query_payload=query_payload,
        request_id=request.request_id,
    )
    response["query_generation"] = _compact_query_payload_for_frontend(response.get("query_generation", {}))
    return response


@app.post("/api/minero/search")
def search(request: SearchRequest) -> Dict[str, Any]:
    _set_progress(
        request.request_id,
        stage="queued",
        processed=0,
        target=request.max_results,
        message="Search queued...",
        done=False,
    )
    query_payload = _generate_query_payload(request.natural_query, request.use_llm)
    ncbi_query = query_payload.get("ncbi_query", "")
    if request.source == "pmc":
        import re as _re
        ncbi_query = _re.sub(r"\[Organism\]", "[all]", ncbi_query)
        ncbi_query = _re.sub(r"\[All Fields\]", "[all]", ncbi_query)
        query_payload["ncbi_query"] = ncbi_query
    response = _execute_search(
        source=request.source,
        max_results=request.max_results,
        use_llm_requested=request.use_llm,
        ncbi_query=ncbi_query,
        query_payload=query_payload,
        request_id=request.request_id,
    )
    response["query_generation"] = _compact_query_payload_for_frontend(response.get("query_generation", {}))
    return response


# ================================================================
# ENDPOINT: MODO PRO — ANÁLISIS PROFUNDO CON PAPER_ANALYZER
# ================================================================

class AnalyzePapersRequest(BaseModel):
    """Recibe la lista de publicaciones PubMed ya obtenidas y una query de contexto."""
    publications: List[Dict[str, Any]] = Field(
        description="Lista de registros PubMed (tal como los devuelve el endpoint de búsqueda)"
    )
    query: str = Field(default="", description="Query original del usuario (para contexto del LLM)")


@app.post("/api/minero/analyze-papers")
def analyze_papers(request: AnalyzePapersRequest) -> Dict[str, Any]:
    """
    Modo Pro: ejecuta PaperAnalyzer sobre una lista de publicaciones PubMed.
    Requiere Ollama con el modelo configurado en Data_Analyzer/config.py.
    Devuelve la lista de resultados con las 53 columnas de metadatos experimentales.
    """
    try:
        # Import dinámico para no bloquear arranque si Data_Analyzer tiene dependencias opcionales
        from paper_analyzer import PaperAnalyzer
        import ollama_client as da_ollama_module
    except ImportError as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Data_Analyzer no disponible: {exc}"
        )

    # Verificar que Ollama está disponible
    try:
        pro_ollama = da_ollama_module.OllamaClient()
        if not pro_ollama.is_available():
            raise HTTPException(
                status_code=503,
                detail="Ollama no está disponible. Verifica que está corriendo y que el modelo está descargado."
            )
    except Exception as exc:
        raise HTTPException(status_code=503, detail=f"Error conectando a Ollama: {exc}")

    if not request.publications:
        raise HTTPException(status_code=400, detail="La lista de publicaciones está vacía.")

    logger.info("[Pro] Iniciando análisis Pro de %d publicaciones...", len(request.publications))

    # Construir el objeto fetcher_data en el formato que espera PaperAnalyzer
    fetcher_data = {
        "metadata": {
            "query": request.query or "(web analysis)",
            "total_results": len(request.publications),
        },
        "publications": request.publications,
    }

    analyzer = PaperAnalyzer(pro_ollama)

    try:
        results = analyzer.analyze_papers(fetcher_data)
    except Exception as exc:
        logger.exception("[Pro] Error durante análisis: %s", exc)
        raise HTTPException(status_code=500, detail=f"Error durante el análisis Pro: {exc}")

    stats = analyzer.stats
    logger.info(
        "[Pro] Análisis completo: %d analizados, %d relevantes, %d errores",
        stats.get("analyzed", 0),
        stats.get("relevant", 0),
        stats.get("errors", 0),
    )

    network_html = ""
    try:
        from paper_visualizer import PaperVisualizer
        visualizer = PaperVisualizer()
        visualizer.set_data(
            query=request.query,
            publications=request.publications,
            classified_results=results
        )
        network_html = visualizer.generate_html_content()
    except Exception as exc:
        logger.exception("[Pro] Error generando visualización: %s", exc)

    return {
        "metadata": {
            "total_analyzed": stats.get("analyzed", 0),
            "total_relevant": stats.get("relevant", 0),
            "total_errors": stats.get("errors", 0),
            "pmc_full_text_used": stats.get("pmc_full_text", 0),
            "abstract_only": stats.get("abstract_only", 0),
            "techniques_enhanced": stats.get("techniques_enhanced", 0),
            "model": pro_ollama.model,
            "analyzed_at": datetime.now(timezone.utc).isoformat(),
        },
        "results": results,
        "network_html": network_html,
    }


if __name__ == "__main__":
    import uvicorn

    uvicorn.run("main:app", host="0.0.0.0", port=8010, reload=True)
