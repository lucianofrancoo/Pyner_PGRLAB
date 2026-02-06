"""
Pyner Phase 2 - Main Orchestrator
==================================
Orquesta el Query Optimizer: Generation → Embedding → Retrieval

Pipeline:
1. Load Phase 1 Knowledge Base
2. Generate queries from KB
3. Create vector embeddings
4. Build FAISS index
5. Ready for search
"""

import json
import logging
import sys
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))

from config import (
    KB_STAGE3, LOG_LEVEL, LOG_FILE, LOGS_DIR,
    VALIDATE_KB_ON_STARTUP, MIN_ORGANISMS_REQUIRED,
    REQUIRED_KB_FIELDS
)
from scripts.query_builder import QueryBuilder
from scripts.vector_db import VectorDatabase, Retriever

# ============================================
# LOGGING SETUP
# ============================================
def setup_logging():
    """Configurar logging"""
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=getattr(logging, LOG_LEVEL),
        format="[%(asctime)s] %(levelname)-8s | %(name)s | %(message)s",
        handlers=[
            logging.FileHandler(LOG_FILE),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


logger = setup_logging()

# ============================================
# PHASE 2 ORCHESTRATOR
# ============================================
class Phase2Orchestrator:
    """Orquesta Phase 2: Query Optimizer"""
    
    def __init__(self):
        """Inicializar orchestrator"""
        self.kb_data = None
        self.query_builder = None
        self.vector_db = None
        self.retriever = None
        self.queries = []
        
        logger.info("="*80)
        logger.info("🚀 PHASE 2: QUERY OPTIMIZER - INITIALIZING")
        logger.info("="*80)
    
    def validate_kb(self) -> bool:
        """Validar Knowledge Base del Phase 1"""
        logger.info("✓ Validating Phase 1 Knowledge Base...")
        
        if not KB_STAGE3.exists():
            logger.error(f"❌ KB not found: {KB_STAGE3}")
            return False
        
        try:
            with open(KB_STAGE3, 'r') as f:
                self.kb_data = json.load(f)
            
            # Check required fields
            for field in REQUIRED_KB_FIELDS:
                if field not in self.kb_data:
                    logger.error(f"❌ Missing required KB field: {field}")
                    return False
            
            # Check organism count
            organisms = self.kb_data.get('organisms', {})
            if len(organisms) < MIN_ORGANISMS_REQUIRED:
                logger.error(f"❌ Insufficient organisms: {len(organisms)} < {MIN_ORGANISMS_REQUIRED}")
                return False
            
            stats = self.kb_data.get('statistics', {})
            logger.info(f"✅ KB validated")
            logger.info(f"   Organisms: {len(organisms)}")
            logger.info(f"   Experiments: {stats.get('total_experiments', 0):,}")
            logger.info(f"   Strategies: {stats.get('unique_strategies', 0)}")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ KB validation failed: {e}")
            return False
    
    def generate_queries(self) -> bool:
        """Generar queries a partir del KB"""
        logger.info("\n📝 Generating queries from KB...")
        
        try:
            self.query_builder = QueryBuilder(KB_STAGE3)
            self.queries = self.query_builder.build_all_queries()
            
            logger.info(f"✅ Generated {len(self.queries)} queries")
            
            # Show breakdown by type
            by_type = {}
            for q in self.queries:
                qtype = q.get('type', 'unknown')
                by_type[qtype] = by_type.get(qtype, 0) + 1
            
            for qtype, count in sorted(by_type.items()):
                logger.info(f"   - {qtype}: {count}")
            
            # Save cache
            self.query_builder.save_query_cache()
            
            return True
            
        except Exception as e:
            logger.error(f"❌ Query generation failed: {e}")
            return False
    
    def build_vector_db(self) -> bool:
        """Construir base de datos de vectores"""
        logger.info("\n🔨 Building vector database...")
        
        try:
            self.vector_db = VectorDatabase()
            
            # Create index
            if not self.vector_db.create_index(self.queries):
                logger.error("❌ Failed to create vector index")
                return False
            
            # Save index
            if not self.vector_db.save():
                logger.error("❌ Failed to save vector index")
                return False
            
            logger.info(f"✅ Vector database created and saved")
            
            # Get stats
            stats = self.vector_db.get_stats()
            logger.info(f"   Total vectors: {stats['total_vectors']}")
            logger.info(f"   Embedding dimension: {stats['embedding_dimension']}")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ Vector DB build failed: {e}")
            return False
    
    def initialize_retriever(self) -> bool:
        """Inicializar retriever"""
        logger.info("\n🔎 Initializing retriever...")
        
        try:
            if not self.vector_db:
                logger.error("❌ Vector DB not initialized")
                return False
            
            self.retriever = Retriever(self.vector_db)
            logger.info(f"✅ Retriever initialized")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ Retriever initialization failed: {e}")
            return False
    
    def run_pipeline(self) -> bool:
        """Ejecutar pipeline completo"""
        logger.info("\n🔄 Starting Phase 2 pipeline...\n")
        
        # Step 1: Validate KB
        if not self.validate_kb():
            return False
        
        # Step 2: Generate queries
        if not self.generate_queries():
            return False
        
        # Step 3: Build vector DB
        if not self.build_vector_db():
            return False
        
        # Step 4: Initialize retriever
        if not self.initialize_retriever():
            return False
        
        logger.info("\n" + "="*80)
        logger.info("✅ PHASE 2 PIPELINE COMPLETE")
        logger.info("="*80)
        
        return True
    
    def test_retrieval(self, test_queries: List[str] = None) -> None:
        """Test retrieval functionality"""
        if not self.retriever:
            logger.warning("❌ Retriever not initialized")
            return
        
        if not test_queries:
            test_queries = [
                "COVID-19 research",
                "Gene expression in humans",
                "Bacterial genomics",
                "Metagenomic studies",
                "RNA sequencing techniques"
            ]
        
        logger.info("\n" + "="*80)
        logger.info("🔎 TEST RETRIEVAL")
        logger.info("="*80)
        
        for test_query in test_queries[:3]:  # Test first 3
            logger.info(f"\n🔍 Query: {test_query}")
            results = self.retriever.retrieve(test_query, top_k=3)
            
            for i, result in enumerate(results, 1):
                score = result['similarity_score']
                text = result['query_text'][:60]
                logger.info(f"   {i}. [{score:.3f}] {text}...")
        
        logger.info("\n" + "="*80)
    
    def get_status(self) -> Dict[str, Any]:
        """Obtener estado actual"""
        status = {
            "kb_loaded": self.kb_data is not None,
            "queries_generated": len(self.queries),
            "vector_db_ready": self.vector_db is not None and self.vector_db.index is not None,
            "retriever_ready": self.retriever is not None,
            "timestamp": str(datetime.now())
        }
        return status


# ============================================
# MAIN EXECUTION
# ============================================
def main():
    """Main entry point"""
    orchestrator = Phase2Orchestrator()
    
    # Run pipeline
    success = orchestrator.run_pipeline()
    
    if success:
        # Test retrieval
        orchestrator.test_retrieval()
        
        # Print status
        status = orchestrator.get_status()
        logger.info(f"\n📊 Phase 2 Status:")
        for key, val in status.items():
            logger.info(f"   {key}: {val}")
    
    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
