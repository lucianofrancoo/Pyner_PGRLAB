"""
Pyner Phase 3 - FastAPI Server
===============================
Production REST API for Pyner

Endpoints:
- GET  / - Health check
- POST /search - Semantic search
- POST /expand - Query expansion
- GET  /stats - System statistics
"""

import logging
import sys
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime
import json

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

try:
    from fastapi import FastAPI, HTTPException
    from fastapi.responses import JSONResponse
    from pydantic import BaseModel
    FASTAPI_AVAILABLE = True
except ImportError:
    FASTAPI_AVAILABLE = False
    logging.error("⚠️ FastAPI not installed")

from phase3.config import (
    API_HOST, API_PORT, VECTOR_DB_PATH, QUERY_CACHE_PATH,
    TOP_K_RESULTS, LOG_FILE, LOGS_DIR
)
from phase3.api.ollama_integration import QueryExpander, OllamaClient

# ============================================
# LOGGING
# ============================================
LOGS_DIR.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)-8s | %(name)s | %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# ============================================
# REQUEST/RESPONSE MODELS
# ============================================
class SearchRequest(BaseModel):
    query: str
    top_k: int = TOP_K_RESULTS
    expand: bool = True

class SearchResult(BaseModel):
    query_text: str
    query_type: str
    similarity_score: float
    rank: int

class SearchResponse(BaseModel):
    query: str
    expanded_queries: List[str]
    results: List[SearchResult]
    total_results: int
    execution_time: float

class StatsResponse(BaseModel):
    status: str
    vector_db_ready: bool
    ollama_available: bool
    queries_cached: int
    timestamp: str

# ============================================
# FASTAPI APP
# ============================================
if FASTAPI_AVAILABLE:
    app = FastAPI(
        title="Pyner Semantic Search",
        description="Production API for NCBI BioProject Research",
        version="3.0.0"
    )
    
    # Global state
    vector_db = None
    retriever = None
    expander = None
    queries_cache = []
    
    @app.on_event("startup")
    async def startup():
        """Inicializar al arrancar"""
        global vector_db, retriever, expander, queries_cache
        
        logger.info("="*80)
        logger.info("🚀 PHASE 3: PRODUCTION API - STARTUP")
        logger.info("="*80)
        
        # Load Vector DB
        try:
            from phase2.scripts.vector_db import VectorDatabase, Retriever
            vector_db = VectorDatabase()
            if vector_db.load():
                retriever = Retriever(vector_db)
                logger.info(f"✅ Vector DB loaded: {VECTOR_DB_PATH}")
            else:
                logger.error("❌ Failed to load Vector DB")
        except Exception as e:
            logger.error(f"❌ Error loading Vector DB: {e}")
        
        # Load Query Cache
        try:
            with open(QUERY_CACHE_PATH, 'r') as f:
                cache_data = json.load(f)
                queries_cache = cache_data.get('queries', [])
                logger.info(f"✅ Query cache loaded: {len(queries_cache)} queries")
        except Exception as e:
            logger.warning(f"⚠️ Query cache not available: {e}")
        
        # Initialize Expander
        expander = QueryExpander()
        logger.info(f"✅ Query Expander initialized")
        
        logger.info("✅ PHASE 3 STARTUP COMPLETE\n")
    
    # ============================================
    # ENDPOINTS
    # ============================================
    
    @app.get("/", response_class=JSONResponse)
    async def health_check():
        """Health check endpoint"""
        return {
            "status": "ok",
            "service": "Pyner Semantic Search Phase 3",
            "version": "3.0.0",
            "timestamp": str(datetime.now())
        }
    
    @app.post("/search", response_model=SearchResponse)
    async def search(request: SearchRequest):
        """Semantic search endpoint"""
        import time
        start_time = time.time()
        
        if not retriever:
            raise HTTPException(status_code=503, detail="Vector DB not initialized")
        
        try:
            # Expand query if requested
            expanded_queries = [request.query]
            if request.expand and expander:
                expanded_queries = expander.expand(request.query)
            
            # Search with all variations
            all_results = []
            for exp_query in expanded_queries:
                results = retriever.retrieve(exp_query, top_k=request.top_k)
                all_results.extend(results)
            
            # Deduplicate and sort by score
            unique_results = {}
            for r in all_results:
                key = r['query_text']
                if key not in unique_results or r['similarity_score'] > unique_results[key]['similarity_score']:
                    unique_results[key] = r
            
            sorted_results = sorted(unique_results.values(), 
                                   key=lambda x: x['similarity_score'], 
                                   reverse=True)[:request.top_k]
            
            # Format response
            search_results = [
                SearchResult(
                    query_text=r['query_text'],
                    query_type=r['query_type'],
                    similarity_score=round(r['similarity_score'], 3),
                    rank=i+1
                )
                for i, r in enumerate(sorted_results)
            ]
            
            execution_time = round(time.time() - start_time, 3)
            
            return SearchResponse(
                query=request.query,
                expanded_queries=expanded_queries,
                results=search_results,
                total_results=len(search_results),
                execution_time=execution_time
            )
            
        except Exception as e:
            logger.error(f"❌ Search error: {e}")
            raise HTTPException(status_code=500, detail=str(e))
    
    @app.post("/expand")
    async def expand_query(query: str):
        """Query expansion endpoint"""
        if not expander:
            return {"query": query, "variations": [query]}
        
        variations = expander.expand(query)
        return {
            "query": query,
            "variations": variations,
            "count": len(variations)
        }
    
    @app.get("/stats", response_model=StatsResponse)
    async def get_stats():
        """System statistics endpoint"""
        ollama_client = OllamaClient()
        
        return StatsResponse(
            status="ready",
            vector_db_ready=vector_db is not None and vector_db.index is not None,
            ollama_available=ollama_client.health_check(),
            queries_cached=len(queries_cache),
            timestamp=str(datetime.now())
        )


# ============================================
# MAIN
# ============================================
if __name__ == "__main__":
    if not FASTAPI_AVAILABLE:
        print("❌ FastAPI not available. Install: pip install fastapi uvicorn")
        sys.exit(1)
    
    import uvicorn
    
    logger.info(f"\n🚀 Starting Pyner Phase 3 API Server")
    logger.info(f"   Host: {API_HOST}")
    logger.info(f"   Port: {API_PORT}")
    logger.info(f"   URL: http://{API_HOST}:{API_PORT}\n")
    
    uvicorn.run(
        app,
        host=API_HOST,
        port=API_PORT,
        log_level="info",
        access_log=True
    )
