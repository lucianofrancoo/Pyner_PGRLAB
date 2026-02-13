"""
Pyner Phase 3 - FastAPI Server (Reimplemented)
===============================================
API REST para generación de queries NCBI

Endpoints:
- GET  / - Health check
- POST /generate - Generar query NCBI desde lenguaje natural
- GET  /stats - Estadísticas del sistema
"""

import logging
import sys
from pathlib import Path
from typing import Dict, List
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

try:
    from fastapi import FastAPI, HTTPException
    from fastapi.responses import JSONResponse
    from pydantic import BaseModel
    FASTAPI_AVAILABLE = True
except ImportError:
    FASTAPI_AVAILABLE = False
    logging.error("⚠️ FastAPI not installed")

from phase3.config import API_HOST, API_PORT, LOG_FILE, LOGS_DIR, SUPPORT_DICT_DIR
from phase3.api.query_generator import QueryGeneratorService
from phase3.api.ollama_integration import OllamaClient
from phase1.config import OUTPUT_DIR

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
# DEPENDENCY CHECKS
# ============================================
def check_dependencies() -> Dict[str, bool]:
    """
    Verify all required dependencies and services before startup.
    Returns dict with status of each dependency.
    """
    status = {
        "kb_file": False,
        "ollama": False,
        "query_cache": False,
        "technical_vocab": False,
        "all_critical": False
    }
    
    print("\n" + "="*80)
    print("🔍 PYNER QUERY GENERATOR - PRE-FLIGHT CHECK")
    print("="*80)
    
    # 1. Check Knowledge Base (CRITICAL)
    kb_path = OUTPUT_DIR / "stage3_kb_reduced.json"
    if kb_path.exists():
        kb_size = kb_path.stat().st_size / (1024*1024)  # MB
        print(f"✅ Knowledge Base found: {kb_path.name} ({kb_size:.1f} MB)")
        status["kb_file"] = True
    else:
        print(f"❌ CRITICAL: Knowledge Base NOT found: {kb_path}")
        print(f"   Expected location: {kb_path}")
        print(f"   → Phase 3 cannot start without KB file")
        print(f"   → Run Phase 1 to generate KB or check file path")
    
    # 2. Check Ollama LLM (OPTIONAL)
    try:
        from phase3.api.ollama_integration import OllamaClient
        ollama = OllamaClient()
        if ollama.is_available():
            print(f"✅ Ollama LLM available at {ollama.host}")
            status["ollama"] = True
        else:
            print(f"⚠️  Ollama NOT available (checked {ollama.host})")
            print(f"   → Query generation will work WITHOUT LLM")
            print(f"   → Using synonym expansion only")
            print(f"   → To enable LLM: Install Ollama (https://ollama.ai)")
    except Exception as e:
        print(f"⚠️  Ollama check failed: {e}")
        print(f"   → Continuing without LLM support")
    
    # 3. Check Query Cache (OPTIONAL - enriches vocabulary)
    query_cache_path = OUTPUT_DIR.parent.parent / "phase2" / "data" / "query_cache.json"
    if query_cache_path.exists():
        cache_size = query_cache_path.stat().st_size / 1024  # KB
        print(f"✅ Query cache found: {query_cache_path.name} ({cache_size:.1f} KB)")
        status["query_cache"] = True
    else:
        print(f"ℹ️  Query cache not found (optional)")
        print(f"   → Vocabulary will work with technical_vocabulary.json only")
    
    # 4. Check Technical Vocabulary (IMPORTANT)
    vocab_path = SUPPORT_DICT_DIR / "technical_vocabulary.json"
    if vocab_path.exists():
        vocab_size = vocab_path.stat().st_size / 1024  # KB
        print(f"✅ Technical vocabulary found: {vocab_path.name} ({vocab_size:.1f} KB)")
        status["technical_vocab"] = True
    else:
        print(f"⚠️  Technical vocabulary not found: {vocab_path}")
        print(f"   → Synonym expansion will be limited")
    
    # 5. Summary
    print("\n" + "-"*80)
    status["all_critical"] = status["kb_file"]
    
    if status["all_critical"]:
        print("✅ ALL CRITICAL DEPENDENCIES SATISFIED")
        if not status["ollama"]:
            print("ℹ️  Running in BASIC mode (no LLM - synonym expansion only)")
        if not status["query_cache"]:
            print("ℹ️  Running without Phase 2 query cache")
        print("="*80 + "\n")
        return status
    else:
        print("❌ CRITICAL DEPENDENCIES MISSING - CANNOT START")
        print("="*80 + "\n")
        return status

# ============================================
# FASTAPI APP
# ============================================
if FASTAPI_AVAILABLE:
    # REQUEST/RESPONSE MODELS (only needed for FastAPI)
    class GenerateQueryRequest(BaseModel):
        query: str
        use_llm: bool = True

    class ExtractedTerms(BaseModel):
        organism: str = None
        strategies: List[str] = []
        tissues: List[str] = []
        conditions: List[str] = []
        free_terms: List[str] = []

    class GenerateQueryResponse(BaseModel):
        user_input: str
        extracted: ExtractedTerms
        ncbi_query: str
        ready_to_use: bool
        clarification_needed: bool = False
        clarification_message: str = ""

    class StatsResponse(BaseModel):
        status: str
        kb_loaded: bool
        organisms_count: int
        strategies_count: int
        ollama_available: bool
        timestamp: str
    
    # FastAPI application
    app = FastAPI(
        title="Pyner Query Generator",
        description="Generador de queries booleanas NCBI desde lenguaje natural",
        version="3.0.0"
    )
    
    # Global state
    query_service = None
    
    @app.on_event("startup")
    async def startup():
        """Inicializar al arrancar"""
        global query_service
        
        logger.info("="*80)
        logger.info("🚀 PHASE 3: QUERY GENERATOR - STARTUP")
        logger.info("="*80)
        
        # Path al KB
        kb_path = OUTPUT_DIR / "stage3_kb_reduced.json"
        
        if not kb_path.exists():
            logger.error(f"❌ KB not found: {kb_path}")
            return
        
        # Inicializar Ollama (opcional)
        ollama_client = None
        try:
            ollama_client = OllamaClient()
            if ollama_client.is_available():
                logger.info("✅ Ollama LLM available")
            else:
                logger.warning("⚠️ Ollama not available - will use basic query generation")
        except Exception as e:
            logger.warning(f"⚠️ Ollama initialization failed: {e}")
        
        # Query cache (para vocabularios adicionales)
        query_cache_path = OUTPUT_DIR.parent.parent / "phase2" / "data" / "query_cache.json"
        
        # Support dictionaries directory for technical vocabulary
        cache_dir = SUPPORT_DICT_DIR

        # Inicializar servicio
        query_service = QueryGeneratorService(
            kb_path,
            ollama_client,
            query_cache_path=query_cache_path,
            cache_dir=cache_dir
        )
        logger.info("✅ Query Generator Service initialized")
        logger.info("✅ PHASE 3 STARTUP COMPLETE\n")
    
    # ============================================
    # ENDPOINTS
    # ============================================
    
    @app.get("/", response_class=JSONResponse)
    async def health_check():
        """Health check endpoint"""
        return {
            "status": "ok",
            "service": "Pyner Query Generator Phase 3",
            "version": "3.0.0",
            "timestamp": str(datetime.now())
        }
    
    @app.post("/generate", response_model=GenerateQueryResponse)
    async def generate_query(request: GenerateQueryRequest):
        """
        Generar query NCBI desde lenguaje natural
        
        Ejemplo:
        {
            "query": "Arabidopsis thaliana drought stress RNA-Seq",
            "use_llm": true
        }
        
        Returns:
        {
            "user_input": "...",
            "extracted": {
                "organism": "Arabidopsis thaliana",
                "strategies": ["RNA-Seq"],
                "free_terms": ["drought", "stress"]
            },
            "ncbi_query": "(\"Arabidopsis thaliana\"[Organism] AND \"drought\"[All Fields] AND \"stress\"[All Fields] AND \"RNA-Seq\"[Strategy])",
            "ready_to_use": true
        }
        """
        
        if not query_service:
            raise HTTPException(status_code=503, detail="Query service not initialized")
        
        try:
            result = query_service.generate_query(request.query, request.use_llm)
            return result
        except Exception as e:
            logger.error(f"❌ Error generating query: {e}")
            raise HTTPException(status_code=500, detail=str(e))
    
    @app.get("/stats", response_model=StatsResponse)
    async def get_stats():
        """Obtener estadísticas del sistema"""
        
        if not query_service:
            return {
                "status": "not_ready",
                "kb_loaded": False,
                "organisms_count": 0,
                "strategies_count": 0,
                "ollama_available": False,
                "timestamp": str(datetime.now())
            }
        
        return {
            "status": "ready",
            "kb_loaded": True,
            "organisms_count": len(query_service.validator.organisms),
            "strategies_count": len(query_service.validator.strategies),
            "ollama_available": query_service.query_builder.ollama_client is not None,
            "timestamp": str(datetime.now())
        }

else:
    logger.error("❌ FastAPI not available - cannot start server")


# ============================================
# CLI INTERFACE (sin servidor)
# ============================================
def generate_query_cli(user_input: str, use_llm: bool = True, interactive: bool = True):
    """
    Generate query from command line
    
    Usage:
        python main.py "Arabidopsis thaliana drought stress RNA-Seq"
        python main.py -i "Arabidopsis thaliana drought stress RNA-Seq"  # Interactive
    """
    from phase1.config import OUTPUT_DIR
    from phase3.api.ollama_integration import OllamaClient
    from phase3.api.interactive_cli import InteractiveQueryGenerator
    
    kb_path = OUTPUT_DIR / "stage3_kb_reduced.json"
    cache_dir = SUPPORT_DICT_DIR
    
    # Initialize
    ollama_client = OllamaClient() if use_llm else None
    
    if interactive:
        # Interactive mode with learning
        generator = InteractiveQueryGenerator(kb_path, ollama_client, cache_dir=cache_dir)
        generator.run_interactive(user_input)
        generator.show_statistics()
    else:
        # Quick mode without interaction
        query_cache_path = OUTPUT_DIR.parent.parent / "phase2" / "data" / "query_cache.json"
        service = QueryGeneratorService(kb_path, ollama_client, query_cache_path=query_cache_path, cache_dir=cache_dir)
        result = service.generate_query(user_input, use_llm)
        
        # Display results
        print("\n" + "="*80)
        print("🔍 PYNER QUERY GENERATOR")
        print("="*80)
        print(f"\n📝 Input: {result['user_input']}")
        print(f"\n🔎 Extracted Terms:")
        
        # Show organism variants
        org_variants = result['extracted'].get('organism_variants', [])
        if org_variants:
            if len(org_variants) == 1:
                print(f"   Organism:    {org_variants[0]}")
            else:
                print(f"   Organisms:   {len(org_variants)} variants found")
                for i, var in enumerate(org_variants, 1):
                    print(f"                {i}. {var}")
        else:
            print(f"   Organism:    (none)")
        
        print(f"   Strategies:  {', '.join(result['extracted']['strategies']) or 'None'}")
        print(f"   Tissues:     {', '.join(result['extracted']['tissues']) or 'None'}")
        print(f"   Conditions:  {', '.join(result['extracted']['conditions']) or 'None'}")
        print(f"   Keywords:    {', '.join(result['extracted']['free_terms']) or 'None'}")

        syn = result.get('synonyms', {})
        org_syn = ', '.join(syn.get('organism', []) or []) or 'None'
        strat_syn = ', '.join(syn.get('strategies', []) or []) or 'None'
        tissue_syn = ', '.join(syn.get('tissues', []) or []) or 'None'
        cond_syn = ', '.join(syn.get('conditions', []) or []) or 'None'
        print(f"\n🔗 Synonyms:")
        print(f"   Organism:    {org_syn}")
        print(f"   Strategies:  {strat_syn}")
        print(f"   Tissues:     {tissue_syn}")
        print(f"   Conditions:  {cond_syn}")
        if result.get('clarification_needed'):
            print("\n⚠️ Need more details:")
            print(f"   {result.get('clarification_message', '')}")
            print("\n" + "="*80)
            print("❓ Try adding a condition, tissue, or strategy")
            print("="*80 + "\n")
        else:
            print(f"\n✅ NCBI Query:")
            print(f"   {result['ncbi_query']}")
            print("\n" + "="*80)
            print("📋 Copy the query above and paste it into NCBI SRA search")
            print("="*80 + "\n")


if __name__ == "__main__":
    import sys
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Pyner Phase 3 - Query Generator (NCBI SRA)",
        usage="python main.py [OPTIONS] QUERY or python main.py --server"
    )
    
    parser.add_argument("query", nargs="?", help="Search query (supports Spanish and English)")
    parser.add_argument("-s", "--server", action="store_true", help="Start API server mode")
    parser.add_argument("-i", "--interactive", action="store_true", default=True, 
                        help="Interactive mode with learning (default)")
    parser.add_argument("-q", "--quick", action="store_true", help="Quick mode without interaction")
    parser.add_argument("--no-llm", action="store_true", help="Disable LLM for faster processing")
    parser.add_argument("--stats", action="store_true", help="Show dictionary statistics")
    
    args = parser.parse_args()
    
    # ============================================
    # PRE-FLIGHT CHECK
    # ============================================
    dep_status = check_dependencies()
    
    if not dep_status["all_critical"]:
        print("\n❌ Cannot continue - critical dependencies missing")
        print("Please ensure Phase 1 has been executed and KB file exists")
        sys.exit(1)
    
    # ============================================
    # START REQUESTED MODE
    # ============================================
    if args.server:
        # Start API server
        import uvicorn
        print("\n🚀 Starting API server...")
        print(f"   Host: {API_HOST}:{API_PORT}")
        print("   URL: http://localhost:8000")
        print("   Docs: http://localhost:8000/docs\n")
        uvicorn.run(app, host=API_HOST, port=API_PORT)
    
    elif args.stats:
        # Show statistics
        from phase3.api.interactive_cli import InteractiveQueryGenerator
        from phase1.config import OUTPUT_DIR
        kb_path = OUTPUT_DIR / "stage3_kb_reduced.json"
        cache_dir = SUPPORT_DICT_DIR
        generator = InteractiveQueryGenerator(kb_path, cache_dir=cache_dir)
        generator.show_statistics()
    
    elif args.query:
        # CLI mode
        use_llm = not args.no_llm
        interactive_mode = not args.quick
        generate_query_cli(args.query, use_llm=use_llm, interactive=interactive_mode)
    
    else:
        # Interactive REPL mode
        from phase3.api.interactive_cli import InteractiveQueryGenerator
        from phase1.config import OUTPUT_DIR
        kb_path = OUTPUT_DIR / "stage3_kb_reduced.json"
        cache_dir = SUPPORT_DICT_DIR
        generator = InteractiveQueryGenerator(kb_path, cache_dir=cache_dir)

        print("\n" + "="*80)
        print("🔍 PYNER QUERY GENERATOR")
        print("="*80)
        print("Type a search phrase and press Enter.")
        print("Type 'exit', 'quit', 'salir', or 'cancel' to stop.")
        print("Examples:")
        print("  - Arabidopsis thaliana drought stress RNA-Seq")
        print("  - Mus musculus liver transcriptome")
        print("  - quiero saber todo sobre trucha arcoiris")
        print("="*80)

        while True:
            try:
                user_input = input("\n> ").strip()
            except (EOFError, KeyboardInterrupt):
                print("\nExiting.")
                break

            if not user_input:
                continue

            if user_input.lower() in {"exit", "quit", "salir", "cancel"}:
                print("\nExiting.")
                break

            generator.run_interactive(user_input)
