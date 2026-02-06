"""
Pyner Phase 2 - Query Builder
==============================
Genera queries automáticamente a partir del Knowledge Base

Features:
- Parse organisms, strategies, tissues del KB
- Generate queries en múltiples templates
- Batch processing para eficiencia
"""

import json
import logging
from pathlib import Path
from typing import List, Dict, Any, Tuple
from collections import defaultdict
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from config import (
    KB_STAGE3, MAX_QUERIES_PER_ORGANISM, MAX_ORGANISMS_TO_QUERY,
    QUERY_TEMPLATES, BATCH_SIZE_QUERIES, QUERY_CACHE_PATH
)

# ============================================
# LOGGING
# ============================================
logger = logging.getLogger(__name__)

# ============================================
# QUERY BUILDER CLASS
# ============================================
class QueryBuilder:
    """Genera queries a partir del Knowledge Base Phase 1"""
    
    def __init__(self, kb_path: Path):
        """Inicializar QueryBuilder con KB"""
        self.kb_path = kb_path
        self.kb_data = None
        self.organisms = {}
        self.strategies = {}
        self.queries = []
        
        logger.info(f"🔨 QueryBuilder initialized with KB: {kb_path.name}")
    
    def load_kb(self) -> bool:
        """Cargar Knowledge Base"""
        try:
            with open(self.kb_path, 'r') as f:
                self.kb_data = json.load(f)
            
            # Extract organisms and strategies
            self.organisms = self.kb_data.get('organisms', {})
            
            # Get strategies from KB (handle both int count and list formats)
            stats = self.kb_data.get('statistics', {})
            num_strategies = stats.get('unique_strategies', 0)
            if isinstance(num_strategies, int):
                # If it's a count, create dummy strategies
                self.strategies = {f"strategy_{i}": 1 for i in range(min(num_strategies, 10))}
            else:
                # If it's a list, use it
                self.strategies = {strat: 1 for strat in num_strategies}
            
            logger.info(f"✅ KB loaded: {len(self.organisms)} organisms, {len(self.strategies)} strategies")
            return True
            
        except Exception as e:
            logger.error(f"❌ Failed to load KB: {e}")
            return False
    
    def generate_organism_queries(self) -> List[Dict[str, str]]:
        """Generar queries para los organismos principales"""
        queries = []
        
        # Top organisms by experiment count
        sorted_organisms = sorted(
            self.organisms.items(), 
            key=lambda x: x[1]['count'], 
            reverse=True
        )[:MAX_ORGANISMS_TO_QUERY]
        
        for organism, stats in sorted_organisms:
            # Template 1: Basic organism query
            query = {
                "type": "organism",
                "text": f"Research studies on {organism}",
                "organism": organism,
                "count": stats['count'],
                "studies": stats['studies']
            }
            queries.append(query)
            
            # Template 2: Gene expression query
            if organism not in ["Severe acute respiratory syndrome coronavirus 2", "metagenome"]:
                query = {
                    "type": "gene_expression",
                    "text": f"Gene expression patterns in {organism}",
                    "organism": organism,
                    "count": stats['count'],
                    "studies": stats['studies']
                }
                queries.append(query)
            
            # Template 3: Comparative query
            if len(queries) % 2 == 0:
                query = {
                    "type": "comparative",
                    "text": f"Comparative genomics and {organism}",
                    "organism": organism,
                    "count": stats['count'],
                    "studies": stats['studies']
                }
                queries.append(query)
        
        logger.info(f"✅ Generated {len(queries)} organism queries")
        return queries
    
    def generate_strategy_queries(self) -> List[Dict[str, str]]:
        """Generar queries para las estrategias de secuenciación"""
        queries = []
        strategies = ["RNA-Seq", "WGS", "WES", "ChIP-Seq", "16S", "metagenomic", "amplicon"]
        
        for strategy in strategies:
            query = {
                "type": "strategy",
                "text": f"Studies using {strategy} sequencing",
                "strategy": strategy
            }
            queries.append(query)
        
        logger.info(f"✅ Generated {len(queries)} strategy queries")
        return queries
    
    def generate_disease_queries(self) -> List[Dict[str, str]]:
        """Generar queries para enfermedades comunes"""
        queries = []
        diseases = [
            ("COVID-19", "Homo sapiens"),
            ("Malaria", "Plasmodium falciparum"),
            ("Tuberculosis", "Mycobacterium tuberculosis"),
            ("Cholera", "Vibrio cholerae"),
            ("Gut dysbiosis", "human gut metagenome"),
        ]
        
        for disease, organism in diseases:
            query = {
                "type": "disease",
                "text": f"Molecular mechanisms of {disease} in {organism}",
                "disease": disease,
                "organism": organism
            }
            queries.append(query)
        
        logger.info(f"✅ Generated {len(queries)} disease queries")
        return queries
    
    def build_all_queries(self) -> List[Dict[str, str]]:
        """Construir todas las queries"""
        if not self.load_kb():
            logger.error("❌ Cannot build queries without KB")
            return []
        
        all_queries = []
        all_queries.extend(self.generate_organism_queries())
        all_queries.extend(self.generate_strategy_queries())
        all_queries.extend(self.generate_disease_queries())
        
        self.queries = all_queries
        
        logger.info(f"📊 Total queries generated: {len(all_queries)}")
        return all_queries
    
    def save_query_cache(self) -> bool:
        """Guardar queries en caché"""
        try:
            cache_data = {
                "timestamp": str(Path(self.kb_path).stat().st_mtime),
                "total_queries": len(self.queries),
                "queries": self.queries
            }
            
            with open(QUERY_CACHE_PATH, 'w') as f:
                json.dump(cache_data, f, indent=2)
            
            logger.info(f"✅ Query cache saved: {QUERY_CACHE_PATH}")
            return True
            
        except Exception as e:
            logger.error(f"❌ Failed to save query cache: {e}")
            return False
    
    def get_queries_by_type(self, query_type: str) -> List[Dict[str, str]]:
        """Obtener queries por tipo"""
        return [q for q in self.queries if q.get("type") == query_type]


# ============================================
# MAIN EXECUTION
# ============================================
def main():
    """Test QueryBuilder"""
    logging.basicConfig(level=logging.INFO)
    
    builder = QueryBuilder(Path(KB_STAGE3))
    queries = builder.build_all_queries()
    
    print("\n" + "="*80)
    print("📝 SAMPLE QUERIES")
    print("="*80)
    
    # Show sample queries
    for i, q in enumerate(queries[:10], 1):
        print(f"\n{i}. [{q.get('type', 'unknown').upper()}]")
        print(f"   Text: {q.get('text', 'N/A')}")
        if 'organism' in q:
            print(f"   Organism: {q['organism']}")
    
    print(f"\n... ({len(queries) - 10} more queries)")
    print("="*80)
    
    # Save cache
    builder.save_query_cache()


if __name__ == "__main__":
    main()
