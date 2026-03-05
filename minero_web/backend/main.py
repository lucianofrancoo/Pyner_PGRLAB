from __future__ import annotations

import logging
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field

ROOT_DIR = Path(__file__).resolve().parents[2]
FETCHER_DIR = ROOT_DIR / "Fetcher_NCBI"
PHASES_DIR = ROOT_DIR / "Query_generator" / "phases"
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


class SearchRequest(BaseModel):
    natural_query: str = Field(min_length=3, description="Consulta biologica en lenguaje natural")
    source: Literal["pubmed", "pmc", "bioproject"] = "pubmed"
    max_results: int = Field(default=20, ge=1, le=10000)
    use_llm: bool = True


class GenerateQueryRequest(BaseModel):
    natural_query: str = Field(min_length=3, description="Consulta biologica en lenguaje natural")
    use_llm: bool = True


class RunSearchRequest(BaseModel):
    source: Literal["pubmed", "pmc", "bioproject"] = "pubmed"
    max_results: int = Field(default=20, ge=1, le=10000)
    use_llm: bool = True
    ncbi_query: str = Field(min_length=3, description="Query booleana NCBI generada")
    query_generation: Dict[str, Any] = Field(default_factory=dict)


class AppState:
    query_service: Optional[QueryGeneratorService] = None
    classifier: Optional[ResultClassifier] = None
    llm_client: Optional[OllamaClient] = None


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


def _ensure_services() -> tuple[QueryGeneratorService, ResultClassifier]:
    if not state.query_service or not state.classifier:
        raise HTTPException(status_code=503, detail="Servicio no inicializado")
    return state.query_service, state.classifier


def _generate_query_payload(natural_query: str, use_llm_requested: bool) -> Dict[str, Any]:
    query_service, classifier = _ensure_services()
    use_llm = use_llm_requested and classifier.llm_available

    query_payload = query_service.generate_query(natural_query, use_llm=use_llm)
    if not query_payload.get("ready_to_use"):
        message = query_payload.get("clarification_message") or "La consulta requiere mas contexto"
        raise HTTPException(status_code=400, detail=message)

    return query_payload


def _execute_search(
    source: Literal["pubmed", "pmc", "bioproject"],
    max_results: int,
    use_llm_requested: bool,
    ncbi_query: str,
    query_payload: Dict[str, Any],
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

    raw_results: List[Dict[str, Any]] = []
    results: List[Dict[str, Any]] = []

    if source in ("pubmed", "pmc"):
        # PMC usa la misma interfaz pero busca en full-text (db='pmc')
        db = "pmc" if source == "pmc" else "pubmed"
        fetcher = LinkoutFetcher()
        raw_results = fetcher.search_publications_by_boolean_query(
            ncbi_query, max_results=max_results, db=db
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
        raw_results = fetcher.run_workflow(ncbi_query, max_bioproject=max_results)
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
    return {"query_generation": query_payload}


@app.post("/api/minero/run-search")
def run_search(request: RunSearchRequest) -> Dict[str, Any]:
    query_payload = request.query_generation or {}
    query_payload["ncbi_query"] = request.ncbi_query
    # Si la fuente es PMC, normalizar field tags para PMC antes de buscar
    ncbi_query = request.ncbi_query
    if request.source == "pmc":
        import re as _re
        ncbi_query = _re.sub(r"\[Organism\]", "[all]", ncbi_query)
        ncbi_query = _re.sub(r"\[All Fields\]", "[all]", ncbi_query)
        query_payload["ncbi_query"] = ncbi_query
    return _execute_search(
        source=request.source,
        max_results=request.max_results,
        use_llm_requested=request.use_llm,
        ncbi_query=ncbi_query,
        query_payload=query_payload,
    )


@app.post("/api/minero/search")
def search(request: SearchRequest) -> Dict[str, Any]:
    query_payload = _generate_query_payload(request.natural_query, request.use_llm)
    ncbi_query = query_payload.get("ncbi_query", "")
    if request.source == "pmc":
        import re as _re
        ncbi_query = _re.sub(r"\[Organism\]", "[all]", ncbi_query)
        ncbi_query = _re.sub(r"\[All Fields\]", "[all]", ncbi_query)
        query_payload["ncbi_query"] = ncbi_query
    return _execute_search(
        source=request.source,
        max_results=request.max_results,
        use_llm_requested=request.use_llm,
        ncbi_query=ncbi_query,
        query_payload=query_payload,
    )


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
