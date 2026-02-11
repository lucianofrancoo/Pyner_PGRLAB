"""
Pyner Phase 2 - FAISS Vector Indexing
======================================
Construye índice vectorial FAISS para búsqueda semántica rápida

Objetivo:
- Cargar metadatos de Phase 1
- Generar embeddings con sentence-transformers
- Construir y guardar índice FAISS para retrieval O(log n)

Ejecución:
    python scripts/phase2_faiss_indexing.py --model all-minilm-l6-v2 --index-type ivf

Performance:
- Vectorization: ~5-10 min (GPU)
- FAISS building: ~2-5 min (CPU)
- Total: ~10-15 min

Output:
- phase2_vectors.faiss        (2.1 GB FAISS index)
- vectors.npy                 (2.0 GB raw vectors)
- index_metadata.json         (config + mapping)
"""

import json
import numpy as np
import logging
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple
import sys
from collections import defaultdict

try:
    import faiss
    FAISS_AVAILABLE = True
except ImportError:
    FAISS_AVAILABLE = False
    print("⚠️  FAISS not installed. Install: pip install faiss-gpu or faiss-cpu")

try:
    from sentence_transformers import SentenceTransformer
    SENTENCE_TRANSFORMERS_AVAILABLE = True
except ImportError:
    SENTENCE_TRANSFORMERS_AVAILABLE = False
    print("⚠️  sentence-transformers not installed. Install: pip install sentence-transformers")

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import KB_STAGE3, DATA_DIR, LOGS_DIR, EMBEDDING_MODEL, USE_GPU


# ============================================
# UTILITY FUNCTIONS
# ============================================

def setup_logging(name: str):
    """Setup basic logging"""
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger(name)

def print_section_header(logger, title: str):
    """Print section header"""
    logger.info("=" * 80)
    logger.info(f"  {title}")
    logger.info("=" * 80)

def print_stage_summary(logger, title: str, stats: dict, duration: float):
    """Print stage summary"""
    logger.info("=" * 80)
    logger.info(f"✅ {title} - Summary")
    logger.info(f"   Duration: {duration:.1f} seconds")
    for key, value in stats.items():
        logger.info(f"   {key}: {value}")
    logger.info("=" * 80)


# ============================================
# CONFIGURATION
# ============================================

EMBEDDING_MODELS = {
    'all-minilm-l6-v2': {
        'name': 'sentence-transformers/all-MiniLM-L6-v2',
        'dimensions': 384,
        'size_mb': 44,
        'speed': 'fast',
        'description': 'General-purpose, lightweight'
    },
    'all-mpnet-base-v2': {
        'name': 'sentence-transformers/all-mpnet-base-v2',
        'dimensions': 768,
        'size_mb': 420,
        'speed': 'medium',
        'description': 'Higher quality but slower'
    },
    'specter': {
        'name': 'allenai/specter',
        'dimensions': 768,
        'size_mb': 500,
        'speed': 'slow',
        'description': 'Scientific papers (best if we have abstracts)'
    }
}

INDEX_TYPES = {
    'flat': {
        'description': 'Exact search O(n), slow but perfect recall',
        'search_time_ms': 100
    },
    'ivf': {
        'description': 'Approximate O(log n), balanced, recommended',
        'search_time_ms': 30
    },
    'ivf-pq': {
        'description': 'Ultra-compressed, very fast, lower recall',
        'search_time_ms': 15
    },
    'hnsw': {
        'description': 'Best recall, medium speed',
        'search_time_ms': 40
    }
}


# ============================================
# VECTORIZATION
# ============================================

def load_phase1_data(logger) -> Tuple[List[str], Dict]:
    """Cargar datos de Phase 1"""
    
    logger.info("📂 Cargando datos de Phase 1...")
    kb_file = KB_STAGE3
    
    if not kb_file.exists():
        logger.error(f"❌ No encontrado: {kb_file}")
        logger.error("   Asegúrate que Phase 1 haya completado exitosamente")
        return [], {}
    
    with open(kb_file, 'r', encoding='utf-8') as f:
        kb_data = json.load(f)
    
    logger.info(f"✅ Metadatos cargados: {len(kb_data)} records")
    
    # Crear descripciones concatenadas para cada record
    # Estas descripciones se vectorizarán
    descriptions = []
    accessions = []
    metadata_list = []
    
    for org_name, org_data in kb_data.get('organisms', {}).items():
        for strategy, count in kb_data.get('strategies', {}).items():
            # Descripción sintética: "organism strategy tissue"
            desc = f"{org_name} {strategy} sequencing"
            descriptions.append(desc)
            accessions.append(f"{org_name}_{strategy}")
            metadata_list.append({
                'organism': org_name,
                'strategy': strategy,
                'count': count if isinstance(count, int) else count.get('count', 0)
            })
    
    logger.info(f"📊 Total descripciones para vectorizar: {len(descriptions):,}")
    
    return descriptions, {
        'descriptions': descriptions,
        'accessions': accessions,
        'metadata': metadata_list,
        'kb_data': kb_data
    }


def vectorize_descriptions(descriptions: List[str], 
                          model_name: str,
                          batch_size: int = 1000,
                          logger: logging.Logger = None) -> np.ndarray:
    """Vectorizar descripciones con sentence-transformers"""
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info(f"🤖 Cargando modelo: {model_name}")
    
    # Cargar modelo
    if not SENTENCE_TRANSFORMERS_AVAILABLE:
        logger.error("❌ sentence-transformers no disponible")
        return None
    
    model = SentenceTransformer(model_name)
    
    # GPU acceleration
    if torch.cuda.is_available() and USE_GPU:
        device = 'cuda:0'
        model = model.to(device)
        logger.info(f"🎮 GPU: {device}")
    else:
        logger.warning("⚠️  GPU no disponible, usando CPU (lento)")
    
    # Vectorizar en batches
    logger.info(f"⚡ Vectorizando {len(descriptions):,} descripciones (batch_size={batch_size})...")
    
    vectors = []
    start_time = time.time()
    
    for i in range(0, len(descriptions), batch_size):
        batch = descriptions[i:i+batch_size]
        batch_vecs = model.encode(batch, batch_size=batch_size, 
                                  show_progress_bar=(i % (10*batch_size) == 0),
                                  convert_to_numpy=True)
        vectors.extend(batch_vecs)
        
        if (i + batch_size) % (10*batch_size) == 0:
            elapsed = time.time() - start_time
            processed = i + batch_size
            rate = processed / elapsed
            remaining_sec = (len(descriptions) - processed) / rate
            logger.info(
                f"  ✓ {processed:,}/{len(descriptions):,} "
                f"({100*processed/len(descriptions):.1f}%) "
                f"| Rate: {rate:.0f} vecs/sec "
                f"| ETA: {remaining_sec//60:.0f}m"
            )
    
    vectors = np.array(vectors, dtype=np.float32)
    
    elapsed = time.time() - start_time
    logger.info(f"✅ Vectorización completada en {elapsed:.1f}s")
    logger.info(f"   Dimensiones: {vectors.shape}")
    logger.info(f"   Throughput: {len(descriptions)/elapsed:.0f} vecs/sec")
    
    return vectors


# ============================================
# FAISS INDEXING
# ============================================

def build_faiss_index(vectors: np.ndarray, 
                     index_type: str = 'ivf',
                     logger: logging.Logger = None) -> faiss.Index:
    """Construir índice FAISS"""
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    if not FAISS_AVAILABLE:
        logger.error("❌ FAISS no disponible")
        return None
    
    n_samples, dim = vectors.shape
    
    logger.info(f"🔨 Normalizando vectores...")
    faiss.normalize_L2(vectors)
    
    logger.info(f"📊 Datos de entrada:")
    logger.info(f"   Vectores: {n_samples:,}")
    logger.info(f"   Dimensiones: {dim}")
    logger.info(f"   Tipo índice: {index_type}")
    
    start_time = time.time()
    
    if index_type == 'flat':
        logger.info("🔍 Construyendo índice: Flat (búsqueda exacta)")
        index = faiss.IndexFlatL2(dim)
        index.add(vectors)
    
    elif index_type == 'ivf':
        logger.info("🔍 Construyendo índice: IVF-Flat (aproximado, balanceado)")
        n_centroids = 512
        quantizer = faiss.IndexFlatL2(dim)
        index = faiss.IndexIVFFlat(quantizer, dim, n_centroids)
        
        logger.info(f"   Training on {min(100000, n_samples):,} vectors...")
        train_vectors = vectors[:min(100000, n_samples)]
        index.train(train_vectors)
        
        logger.info(f"   Adding {n_samples:,} vectors to index...")
        index.add(vectors)
    
    elif index_type == 'hnsw':
        logger.info("🔍 Construyendo índice: HNSW (máxima recall)")
        index_hnsw = faiss.IndexHNSWFlat(dim, 16)
        index = index_hnsw
        index.add(vectors)
    
    else:
        logger.error(f"❌ Tipo de índice desconocido: {index_type}")
        return None
    
    elapsed = time.time() - start_time
    logger.info(f"✅ Índice construido en {elapsed:.1f}s")
    
    return index


def save_faiss_index(index: faiss.Index,
                    vectors: np.ndarray,
                    accessions: List[str],
                    kb_data: Dict,
                    model_info: Dict,
                    index_type: str,
                    logger: logging.Logger = None):
    """Guardar índice y metadatos"""
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info(f"💾 Guardando índice FAISS...")
    
    phase2_dir = DATA_DIR
    phase2_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Índice FAISS
    faiss_path = phase2_dir / "pyner_vectors.faiss"
    faiss.write_index(index, str(faiss_path))
    faiss_size_mb = faiss_path.stat().st_size / (1024*1024)
    logger.info(f"   ✅ {faiss_path.name} ({faiss_size_mb:.1f} MB)")
    
    # 2. Vectores crudos (para debugging)
    vectors_path = phase2_dir / "vectors.npy"
    np.save(vectors_path, vectors)
    vec_size_mb = vectors_path.stat().st_size / (1024*1024)
    logger.info(f"   ✅ {vectors_path.name} ({vec_size_mb:.1f} MB)")
    
    # 3. Metadatos
    metadata = {
        "execution_date": datetime.now().isoformat(),
        "index_size": len(vectors),
        "dimensions": vectors.shape[1],
        "index_type": index_type,
        "embedding_model": model_info['name'],
        "model_dimensions": model_info['dimensions'],
        "total_accessions": len(accessions),
        "kb_stats": {
            "total_experiments": kb_data.get('metadata', {}).get('total_experiments'),
            "unique_organisms": kb_data.get('metadata', {}).get('unique_organisms'),
            "unique_strategies": kb_data.get('metadata', {}).get('unique_strategies'),
        },
        "search_parameters": {
            "index_type": index_type,
            "estimated_search_time_ms": INDEX_TYPES.get(index_type, {}).get('search_time_ms'),
            "estimated_recall": 0.95 if index_type == 'ivf' else 0.99 if index_type == 'hnsw' else 1.0
        }
    }
    
    metadata_path = phase2_dir / "index_metadata.json"
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    logger.info(f"   ✅ {metadata_path.name}")
    
    # 4. Mapeo accesiones
    accession_map = {i: acc for i, acc in enumerate(accessions)}
    map_path = phase2_dir / "accession_mapping.json"
    with open(map_path, 'w') as f:
        json.dump(accession_map, f, indent=2)
    logger.info(f"   ✅ {map_path.name}")
    
    logger.info(f"\n📂 Todos los archivos guardados en: {phase2_dir}")
    
    return {
        'index_path': str(faiss_path),
        'vectors_path': str(vectors_path),
        'metadata_path': str(metadata_path),
        'mapping_path': str(map_path)
    }


# ============================================
# VALIDATION & TESTING
# ============================================

def test_index(index: faiss.Index, vectors: np.ndarray, logger: logging.Logger = None):
    """Probar índice con búsquedas de ejemplo"""
    
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info("\n🧪 Probando índice FAISS con búsquedas de ejemplo...")
    
    # Buscar los KNN del primer vector
    n_test = min(5, len(vectors))
    
    for i in range(n_test):
        query_vec = vectors[i].reshape(1, -1)
        distances, indices = index.search(query_vec, k=5)
        
        logger.debug(f"  Query {i}: Top-5 nearest")
        for idx, dist in zip(indices[0], distances[0]):
            logger.debug(f"    - Index {idx}: distance {dist:.4f}")
    
    logger.info("✅ Test completado")


# ============================================
# MAIN EXECUTION
# ============================================

def main():
    """Pipeline completo Phase 2"""
    
    logger = setup_logging("phase2_faiss")
    
    print_section_header(logger, "PHASE 2: FAISS VECTOR INDEXING")
    
    # Parse arguments
    import argparse
    parser = argparse.ArgumentParser(description="Pyner Phase 2 - FAISS Indexing")
    parser.add_argument('--model', default='all-minilm-l6-v2', 
                       choices=list(EMBEDDING_MODELS.keys()),
                       help='Embedding model')
    parser.add_argument('--index-type', default='ivf',
                       choices=list(INDEX_TYPES.keys()),
                       help='FAISS index type')
    parser.add_argument('--batch-size', type=int, default=1000,
                       help='Batch size for vectorization')
    args = parser.parse_args()
    
    # ============ SETUP ============
    print_section_header(logger, "SETUP: Configuration")
    
    model_info = EMBEDDING_MODELS[args.model]
    logger.info(f"📦 Model: {model_info['name']}")
    logger.info(f"   Dimensions: {model_info['dimensions']}")
    logger.info(f"   Size: {model_info['size_mb']} MB")
    logger.info(f"   Speed: {model_info['speed']}")
    
    logger.info(f"\n🔍 Index type: {args.index_type}")
    logger.info(f"   {INDEX_TYPES[args.index_type]['description']}")
    logger.info(f"   Est. search time: {INDEX_TYPES[args.index_type]['search_time_ms']} ms")
    
    # ============ PHASE 1 DATA ============
    print_section_header(logger, "LOAD: Phase 1 Data")
    
    descriptions, data = load_phase1_data(logger)
    
    if not descriptions:
        logger.error("❌ Failed to load Phase 1 data")
        return 1
    
    logger.info(f"✅ Loaded {len(descriptions):,} descriptions")
    
    # ============ VECTORIZATION ============
    print_section_header(logger, "VECTORIZE: Generate Embeddings")
    
    vectors = vectorize_descriptions(
        descriptions,
        model_info['name'],
        batch_size=args.batch_size,
        logger=logger
    )
    
    if vectors is None:
        logger.error("❌ Vectorization failed")
        return 1
    
    # ============ INDEX BUILDING ============
    print_section_header(logger, "INDEX: Build FAISS")
    
    index = build_faiss_index(vectors, index_type=args.index_type, logger=logger)
    
    if index is None:
        logger.error("❌ Index building failed")
        return 1
    
    # ============ VALIDATION ============
    print_section_header(logger, "TEST: Validate Index")
    
    test_index(index, vectors, logger=logger)
    
    # ============ SAVING ============
    print_section_header(logger, "SAVE: Persist Index")
    
    paths = save_faiss_index(
        index, vectors, data['accessions'], data['kb_data'],
        model_info, args.index_type, logger=logger
    )
    
    # ============ SUMMARY ============
    stats = {
        'Total vectors': f"{len(vectors):,}",
        'Dimensions': f"{vectors.shape[1]}",
        'Model': model_info['name'].split('/')[-1],
        'Index type': args.index_type,
        'Estimated search time': f"{INDEX_TYPES[args.index_type]['search_time_ms']}ms",
        'Output size (FAISS)': f"{Path(paths['index_path']).stat().st_size / (1024*1024):.1f} MB",
    }
    
    print_stage_summary(logger, "PHASE 2: FAISS INDEXING", stats, 0)
    
    logger.info("✅ Phase 2 completado exitosamente")
    logger.info(f"\n📂 Índice guardado en: {Path(paths['index_path']).parent}")
    
    return 0


if __name__ == "__main__":
    exit(main())
