#!/usr/bin/env python3
"""
🧪 DEMO: Query Generator Prototype (moved)
====================================

Propósito: Demostrar generación de queries desde lenguaje natural
Usa: KB local (500K) como proof-of-concept

Nota: Esto es sesgado (solo 500K) pero demuestra el concepto.
      Phase 3b añadirá NCBI API integration.
"""

import sys
import json
import time
from pathlib import Path

# Ensure repo root is on sys.path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

try:
    from phase3.api.ollama_integration import QueryExpander
    from phase2.scripts.vector_db import VectorDatabase, Retriever
    IMPORTS_OK = True
except ImportError as e:
    print(f"⚠️ Import error: {e}")
    IMPORTS_OK = False


class DemoQueryGenerator:
    """Demostración interactiva de generación de queries"""
    
    def __init__(self):
        """Inicializar componentes"""
        print("\n" + "="*80)
        print("🧪 PYNER DEMO: Query Generator (Sesgado - 500K KB)")
        print("="*80)
        
        # Cargar FAISS index
        print("\n📚 Cargando Vector DB...")
        self.vector_db = VectorDatabase()
        if self.vector_db.load():
            self.retriever = Retriever(self.vector_db)
            print(f"   ✅ FAISS Index cargado: {self.vector_db.index.ntotal} vectores")
        else:
            print("   ❌ Error loading FAISS index")
            self.retriever = None
        
        # Cargar KB
        print("\n📖 Cargando Knowledge Base (500K)...")
        self.kb = self._load_kb()
        print(f"   ✅ KB cargado: {len(self.kb.get('organisms', {}))} organismos")
        
        # Inicializar expander
        print("\n🤖 Inicializando Query Expander (Ollama LLM)...")
        self.expander = QueryExpander()
        print("   ✅ Query Expander listo")
        
        self.ready = True
        
    def _load_kb(self):
        """Cargar Knowledge Base"""
        kb_path = Path("phase1/output/stage3_knowledge_base.json")
        if not kb_path.exists():
            print("   ❌ KB no encontrado")
            return {}
        
        with open(kb_path, 'r') as f:
            return json.load(f)
    
    def generate_query(self, natural_language: str, top_k: int = 5) -> dict:
        """
        Generar query desde lenguaje natural
        """
        
        if not self.ready:
            return {"error": "Generator not initialized"}
        
        start_time = time.time()
        
        print(f"\n{'='*80}")
        print(f"❓ INPUT: \"{natural_language}\"")
        print(f"{'='*80}")
        
        # Step 1: Vector search
        print(f"\n⚡ Paso 1: Búsqueda vectorial (FAISS)...")
        if self.retriever:
            results = self.retriever.retrieve(natural_language, top_k=top_k)
            print(f"   ✅ Encontrados {len(results)} resultados similares")
        else:
            print("   ❌ FAISS not available")
            return {"error": "FAISS not available"}
        
        # Step 2: Enrich with KB
        print(f"\n📊 Paso 2: Enriquecimiento con KB...")
        enriched = []
        for i, result in enumerate(results, 1):
            query_text = result['query_text']
            query_type = result['query_type']
            score = result['similarity_score']
            
            # Lookup en KB
            kb_data = self._lookup_kb(query_text, query_type)
            
            enriched_result = {
                'rank': i,
                'query_text': query_text,
                'query_type': query_type,
                'similarity_score': score,
                'kb_data': kb_data
            }
            enriched.append(enriched_result)
            
            print(f"   {i}. [{score:.3f}] {query_type:15s} | {query_text[:50]}...")
            if kb_data:
                print(f"      └─ KB: {kb_data.get('label', 'N/A')}")
        
        execution_time = time.time() - start_time
        
        # Step 3: Format response
        print(f"\n✅ Paso 3: Respuesta formateada")
        
        response = {
            "input": natural_language,
            "generated_queries": results,
            "enriched_results": enriched,
            "kb_coverage": {
                "total_organisms": len(self.kb.get('organisms', {})),
                "total_experiments": self.kb.get('statistics', {}).get('total_experiments', 0),
                "total_samples": self.kb.get('statistics', {}).get('total_samples', 0),
                "note": "⚠️ Sesgado: solo 500K de 3.6M archivos NCBI"
            },
            "execution_time": round(execution_time, 3)
        }
        
        return response
    
    def _lookup_kb(self, query_text: str, query_type: str) -> dict:
        """Buscar en KB para enriquecer resultado"""
        
        if query_type == 'organism':
            for org_name, org_data in self.kb.get('organisms', {}).items():
                if org_name.lower() in query_text.lower() or query_text.lower() in org_name.lower():
                    return {
                        'label': org_name,
                        'experiments': org_data.get('count', 0),
                        'studies': org_data.get('studies', 0),
                        'percentage': round((org_data.get('count', 0) / 
                                           self.kb.get('statistics', {}).get('total_experiments', 1)) * 100, 1)
                    }
        
        elif query_type == 'strategy':
            for strat_name, strat_count in self.kb.get('strategies', {}).items():
                if strat_name.lower() in query_text.lower() or query_text.lower() in strat_name.lower():
                    return {
                        'label': strat_name,
                        'experiments': strat_count,
                        'percentage': round((strat_count / 
                                           self.kb.get('statistics', {}).get('total_experiments', 1)) * 100, 1)
                    }
        
        elif query_type == 'gene_expression':
            return {
                'label': 'Gene Expression',
                'note': 'Múltiples organismos disponibles',
                'cached_patterns': len(self.kb.get('gene_expression', {}))
            }
        
        elif query_type == 'disease':
            return {
                'label': 'Disease Query',
                'cached_diseases': len(self.kb.get('diseases', {}))
            }
        
        return {}
    
    def print_results(self, response: dict):
        """Mostrar resultados de forma presentable"""
        
        print(f"\n{'='*80}")
        print("📊 RESULTADOS")
        print(f"{'='*80}")
        
        # Inputs
        print(f"\n🔍 Consulta: \"{response['input']}\"")
        print(f"⏱️  Tiempo: {response['execution_time']}ms")
        
        # Coverage warning
        coverage = response.get('kb_coverage', {})
        print(f"\n⚠️  COBERTURA DEL KB:")
        print(f"   Organismos: {coverage.get('total_organisms', 0):,}")
        print(f"   Experiments: {coverage.get('total_experiments', 0):,}")
        print(f"   {coverage.get('note', '')}")
        
        # Results
        print(f"\n🎯 TOP RESULTADOS (desde KB 500K):")
        for result in response.get('enriched_results', []):
            print(f"\n   Rank {result['rank']}: [{result['similarity_score']:.3f}]")
            print(f"   Query: {result['query_text']}")
            print(f"   Type: {result['query_type']}")
            
            kb_data = result.get('kb_data', {})
            if kb_data:
                print(f"   KB Data:")
                for key, val in kb_data.items():
                    if isinstance(val, (int, float)):
                        print(f"      • {key}: {val:,}")
                    else:
                        print(f"      • {key}: {val}")
        
        print(f"\n{'='*80}\n")


def interactive_demo():
    """Modo interactivo"""
    
    if not IMPORTS_OK:
        print("❌ Imports faltan. Instala: pip install biopython")
        return
    
    generator = DemoQueryGenerator()
    
    if not generator.ready:
        print("❌ No se pudo inicializar")
        return
    
    print("\n\n🎓 MODO INTERACTIVO")
    print("="*80)
    print("Escribe preguntas en lenguaje natural.")
    print("Ejemplo: 'Virus que infectan humanos'")
    print("Escribe 'exit' para salir.\n")
    
    while True:
        try:
            query = input("\n❓ Pregunta: ").strip()
            
            if query.lower() == 'exit':
                print("\n👋 Adiós!\n")
                break
            
            if not query:
                continue
            
            # Generate
            response = generator.generate_query(query, top_k=5)
            
            if 'error' in response:
                print(f"❌ Error: {response['error']}")
            else:
                generator.print_results(response)
        
        except KeyboardInterrupt:
            print("\n\n👋 Adiós!\n")
            break
        except Exception as e:
            print(f"❌ Error: {e}")


def demo_with_examples():
    """Demostración con ejemplos precargados"""
    
    if not IMPORTS_OK:
        print("❌ Imports faltan. Instala: pip install biopython")
        return
    
    examples = [
        "¿Qué virus infectan a humanos?",
        "Secuenciación de bacteria del suelo",
        "Expresión génica en plantas",
        "Antibiótesistencia en patógenos",
    ]
    
    print("\n\n🎬 DEMO CON EJEMPLOS")
    print("="*80)
    
    generator = DemoQueryGenerator()
    
    if not generator.ready:
        print("❌ No se pudo inicializar")
        return
    
    for example in examples:
        response = generator.generate_query(example, top_k=3)
        
        if 'error' not in response:
            generator.print_results(response)
        
        # Pequeña pausa entre queries
        time.sleep(1)


if __name__ == "__main__":
    
    if len(sys.argv) > 1 and sys.argv[1] == "interactive":
        # Modo interactivo
        interactive_demo()
    else:
        # Demostración con ejemplos
        demo_with_examples()
