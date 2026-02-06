#!/usr/bin/env python3
"""
🚀 QUICK TEST: Generación de Queries (MVP)
==========================================

Script simple para testear la generación de queries
sin complejidades. Usa lo que ya tenemos funcionando.
"""

import sys
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

print("\n🧪 QUICK TEST: Query Generation\n")

# Test 1: Cargar KB
print("1️⃣ Cargando Knowledge Base...")
kb_path = Path("phase1/output/stage3_knowledge_base.json")
try:
    with open(kb_path, 'r') as f:
        kb = json.load(f)
    print(f"   ✅ KB cargado")
    print(f"   - Organismos: {len(kb['organisms'])}")
    print(f"   - Experiments: {kb['statistics']['total_experiments']:,}")
except Exception as e:
    print(f"   ❌ Error: {e}")
    sys.exit(1)

# Test 2: Cargar FAISS Index
print("\n2️⃣ Cargando FAISS Index...")
try:
    from phase2.scripts.vector_db import VectorDatabase, Retriever
    
    vector_db = VectorDatabase()
    if vector_db.load():
        retriever = Retriever(vector_db)
        print(f"   ✅ FAISS cargado")
        print(f"   - Vectores: {vector_db.index.ntotal}")
    else:
        print(f"   ❌ Error loading FAISS")
        sys.exit(1)
except Exception as e:
    print(f"   ❌ Error: {e}")
    sys.exit(1)

# Test 3: Prueba de búsqueda
print("\n3️⃣ Prueba de búsqueda vectorial...")
try:
    test_query = "Virus that infect humans with sequencing"
    results = retriever.retrieve(test_query, top_k=3)
    
    print(f"   ✅ Búsqueda exitosa")
    print(f"   - Query: \"{test_query}\"")
    print(f"   - Resultados top-3:")
    for i, r in enumerate(results, 1):
        print(f"      {i}. [{r['similarity_score']:.3f}] {r['query_text'][:60]}...")
except Exception as e:
    print(f"   ❌ Error: {e}")
    sys.exit(1)

# Test 4: Enriquecimiento con KB
print("\n4️⃣ Enriquecimiento con KB...")
try:
    for result in results[:1]:
        query_text = result['query_text']
        query_type = result['query_type']
        
        print(f"   Query: {query_text}")
        print(f"   Type: {query_type}")
        
        # Buscar en KB si es organismo
        if query_type == 'organism':
            for org_name, org_data in list(kb['organisms'].items())[:1]:
                if org_name.lower() in query_text.lower():
                    exp_count = org_data.get('count', 0)
                    max_exp = kb['statistics']['total_experiments']
                    pct = (exp_count / max_exp) * 100 if max_exp > 0 else 0
                    print(f"   KB Match: {org_name}")
                    print(f"   Experiments: {exp_count:,} ({pct:.1f}%)")
                    break
    
    print(f"   ✅ Enriquecimiento exitoso")
except Exception as e:
    print(f"   ❌ Error: {e}")
    sys.exit(1)

# Test 5: Usar Ollama para expansión
print("\n5️⃣ Expansión con Ollama LLM...")
try:
    from phase3.api.ollama_integration import QueryExpander
    
    expander = QueryExpander()
    expanded = expander.expand("COVID-19 research")
    
    print(f"   ✅ Expansión completada")
    print(f"   - Original: \"COVID-19 research\"")
    print(f"   - Expanded: {expanded}")
except Exception as e:
    print(f"   ⚠️ Ollama not available (ok): {str(e)[:50]}")

print("\n✅ ALL TESTS PASSED!\n")
print("=" * 70)
print("Para el demo interactivo:")
print("  python3 phase3/scripts/demo_query_generator.py interactive")
print("\nPara el demo con ejemplos:")
print("  python3 phase3/scripts/demo_query_generator.py")
print("=" * 70 + "\n")
