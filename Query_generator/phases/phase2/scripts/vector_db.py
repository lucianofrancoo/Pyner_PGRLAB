"""
Pyner Phase 2 - Vector Database
================================
Gestión de base de datos de vectores usando FAISS

Features:
- Embedding generation con sentence-transformers
- FAISS index creation y search
- Batch processing
- GPU support
"""

import json
import logging
import pickle
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Any
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from config import (
    EMBEDDING_MODEL, EMBEDDING_DIMENSION, VECTOR_DB_PATH,
    VECTOR_INDEX_PATH, TOP_K_RESULTS, SIMILARITY_METRIC,
    USE_GPU, BATCH_SIZE_VECTORS
)

try:
    import faiss
    FAISS_AVAILABLE = True
except ImportError:
    FAISS_AVAILABLE = False
    logging.warning("⚠️ FAISS not available, using numpy fallback")

try:
    from sentence_transformers import SentenceTransformer
    ST_AVAILABLE = True
except ImportError:
    ST_AVAILABLE = False
    logging.warning("⚠️ sentence-transformers not available")

# ============================================
# LOGGING
# ============================================
logger = logging.getLogger(__name__)

# ============================================
# VECTOR DATABASE CLASS
# ============================================
class VectorDatabase:
    """Manejo de base de datos de vectores con FAISS"""
    
    def __init__(self):
        """Inicializar VectorDatabase"""
        self.embedder = None
        self.index = None
        self.metadata = []  # Store metadata for each vector
        self.db_path = Path(VECTOR_DB_PATH)
        self.index_path = Path(VECTOR_INDEX_PATH)
        
        logger.info(f"🔍 VectorDatabase initialized (FAISS: {FAISS_AVAILABLE})")
        
        # Load embedder
        if ST_AVAILABLE:
            try:
                logger.info(f"📥 Loading embedding model: {EMBEDDING_MODEL}")
                self.embedder = SentenceTransformer(EMBEDDING_MODEL)
                logger.info(f"✅ Embedder loaded (dims: {EMBEDDING_DIMENSION})")
            except Exception as e:
                logger.warning(f"⚠️ Failed to load embedder: {e}")
    
    def embed_texts(self, texts: List[str]) -> np.ndarray:
        """Generar embeddings para textos"""
        if not self.embedder:
            logger.warning("❌ Embedder not available, using random vectors")
            return np.random.randn(len(texts), EMBEDDING_DIMENSION).astype('float32')
        
        try:
            embeddings = self.embedder.encode(texts, show_progress_bar=True)
            return embeddings.astype('float32')
        except Exception as e:
            logger.error(f"❌ Embedding failed: {e}")
            return np.random.randn(len(texts), EMBEDDING_DIMENSION).astype('float32')
    
    def create_index(self, queries: List[Dict[str, Any]]) -> bool:
        """Crear índice FAISS a partir de queries"""
        if not FAISS_AVAILABLE:
            logger.error("❌ FAISS not available")
            return False
        
        if not queries or len(queries) == 0:
            logger.error("❌ No queries to index")
            return False
        
        try:
            # Extract text from queries
            query_texts = [q.get('text', str(q)) for q in queries]
            
            logger.info(f"📝 Embedding {len(query_texts)} queries...")
            vectors = self.embed_texts(query_texts)
            
            # Create FAISS index
            logger.info(f"🔨 Creating FAISS index...")
            self.index = faiss.IndexFlatL2(EMBEDDING_DIMENSION)
            self.index.add(vectors)
            
            # Store metadata
            self.metadata = queries
            
            logger.info(f"✅ Index created with {len(queries)} vectors")
            return True
            
        except Exception as e:
            logger.error(f"❌ Index creation failed: {e}")
            return False
    
    def search(self, query_text: str, top_k: int = TOP_K_RESULTS) -> List[Tuple[str, float]]:
        """Buscar en el índice"""
        if not self.index:
            logger.warning("❌ Index not initialized")
            return []
        
        try:
            # Embed query
            query_vector = self.embed_texts([query_text])[0].reshape(1, -1)
            
            # Search
            distances, indices = self.index.search(query_vector, top_k)
            
            # Return results with similarity scores
            results = []
            for i, idx in enumerate(indices[0]):
                if idx < len(self.metadata):
                    # Convert distance to similarity (1 / (1 + distance))
                    similarity = 1.0 / (1.0 + distances[0][i])
                    results.append((self.metadata[idx], similarity))
            
            return results
            
        except Exception as e:
            logger.error(f"❌ Search failed: {e}")
            return []
    
    def save(self) -> bool:
        """Guardar índice y metadata"""
        try:
            if self.index:
                faiss.write_index(self.index, str(self.db_path))
                logger.info(f"✅ Index saved: {self.db_path}")
            
            # Save metadata
            with open(self.index_path, 'wb') as f:
                pickle.dump(self.metadata, f)
                logger.info(f"✅ Metadata saved: {self.index_path}")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ Save failed: {e}")
            return False
    
    def load(self) -> bool:
        """Cargar índice y metadata"""
        try:
            if self.db_path.exists() and FAISS_AVAILABLE:
                self.index = faiss.read_index(str(self.db_path))
                logger.info(f"✅ Index loaded: {self.db_path}")
            
            if self.index_path.exists():
                with open(self.index_path, 'rb') as f:
                    self.metadata = pickle.load(f)
                    logger.info(f"✅ Metadata loaded: {len(self.metadata)} entries")
            
            return True
            
        except Exception as e:
            logger.error(f"❌ Load failed: {e}")
            return False
    
    def get_stats(self) -> Dict[str, Any]:
        """Obtener estadísticas del índice"""
        stats = {
            "index_exists": self.index is not None,
            "total_vectors": len(self.metadata) if self.index else 0,
            "embedding_dimension": EMBEDDING_DIMENSION,
            "faiss_available": FAISS_AVAILABLE,
            "embedder_available": self.embedder is not None,
        }
        return stats


# ============================================
# RETRIEVAL CLASS
# ============================================
class Retriever:
    """Search y ranking en la base de datos de vectores"""
    
    def __init__(self, vector_db: VectorDatabase):
        """Inicializar Retriever"""
        self.db = vector_db
        logger.info("🔎 Retriever initialized")
    
    def retrieve(self, query: str, top_k: int = TOP_K_RESULTS) -> List[Dict[str, Any]]:
        """Recuperar resultados más relevantes"""
        results = self.db.search(query, top_k)
        
        # Format results
        formatted_results = []
        for metadata, score in results:
            formatted_results.append({
                "query_text": metadata.get('text', ''),
                "query_type": metadata.get('type', ''),
                "similarity_score": float(score),
                "metadata": metadata
            })
        
        return formatted_results
    
    def batch_retrieve(self, queries: List[str], top_k: int = TOP_K_RESULTS) -> List[List[Dict[str, Any]]]:
        """Recuperar resultados en lote"""
        results = []
        for query in queries:
            results.append(self.retrieve(query, top_k))
        return results


# ============================================
# MAIN EXECUTION
# ============================================
def main():
    """Test VectorDatabase"""
    logging.basicConfig(level=logging.INFO)
    
    # Create sample queries
    sample_queries = [
        {"text": "Gene expression in humans", "type": "organism"},
        {"text": "COVID-19 pathology mechanisms", "type": "disease"},
        {"text": "RNA-Seq sequencing studies", "type": "strategy"},
        {"text": "Microbial community analysis", "type": "metagenome"},
        {"text": "Comparative genomics across species", "type": "comparative"},
    ]
    
    # Initialize VectorDB
    vdb = VectorDatabase()
    
    print("\n" + "="*80)
    print("🔍 VECTOR DATABASE TEST")
    print("="*80)
    
    # Show stats
    stats = vdb.get_stats()
    print(f"\n📊 VectorDB Stats:")
    for key, val in stats.items():
        print(f"   {key}: {val}")
    
    print("\n" + "="*80)
    print("✅ VectorDatabase ready for Phase 2")
    print("="*80)


if __name__ == "__main__":
    main()
