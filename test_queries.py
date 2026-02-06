#!/usr/bin/env python3
"""
🎯 PYNER: Test Query Generation (Simple)
=========================================

Modo rápido sin dependencias complejas.
Genera queries y muestra resultados desde KB (500K).

Uso:
    # Ejemplo específico
    python3 test_queries.py "virus humanos"
    
    # Modo interactivo
    python3 test_queries.py --interactive
"""

import sys
import json
import subprocess
from pathlib import Path

# Simpler approach - use curl to test API
def test_via_api(query: str):
    """Test usando el API que ya está corriendo"""
    
    print(f"\n{'='*70}")
    print(f"❓ Query: \"{query}\"")
    print(f"{'='*70}\n")
    
    # Check if API is running
    import subprocess
    import time
    
    result = subprocess.run(
        ['curl', '-s', 'http://localhost:8000/'],
        capture_output=True,
        text=True,
        timeout=5
    )
    
    if result.returncode != 0:
        print("❌ API no está corriendo en localhost:8000")
        print("   Incia con: uvicorn phase3.api.main:app --host 0.0.0.0 --port 8000")
        return
    
    print("✅ API está disponible\n")
    
    # Test search endpoint
    print("🔍 Búsqueda en KB (FAISS)...\n")
    
    cmd = [
        'curl', '-s', '-X', 'POST',
        'http://localhost:8000/search',
        '-H', 'Content-Type: application/json',
        '-d', json.dumps({
            'query': query,
            'top_k': 5,
            'expand': False
        })
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
    
    if result.returncode != 0:
        print(f"❌ Error: {result.stderr}")
        return
    
    try:
        data = json.loads(result.stdout)
    except:
        print(f"❌ Error parsing response: {result.stdout[:200]}")
        return
    
    # Display results
    print(f"📊 Resultados (Top-{len(data['results'])}):\n")
    
    for i, res in enumerate(data['results'], 1):
        score = res['similarity_score']
        bar = "█" * int(score * 20) + "░" * (20 - int(score * 20))
        
        print(f"  {i}. [{score:.3f}] {bar}")
        print(f"     Type: {res['query_type']}")
        print(f"     Query: {res['query_text']}")
        print()
    
    print(f"⏱️  Tiempo: {data['execution_time']*1000:.1f}ms\n")
    
    # Warning about bias
    print(f"{'='*70}")
    print("⚠️  NOTA: Estos resultados son SESGADOS")
    print("   • Solo usa 500K de 3.6M archivos NCBI")
    print("   • Falta 85.7% de la cobertura")
    print("   • Phase 3b añadirá consulta real a NCBI")
    print(f"{'='*70}\n")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nEjemplos:\n")
        examples = [
            "virus humanos",
            "CRISPR plantas",
            "bacteria suelo",
            "expresión génica"
        ]
        for ex in examples:
            print(f"  python3 test_queries.py \"{ex}\"")
        return
    
    if sys.argv[1] == "--interactive":
        print("\n🎓 MODO INTERACTIVO")
        print("=" * 70)
        print("Escribe queries. Exit para salir.\n")
        
        while True:
            try:
                query = input("❓ Query: ").strip()
                if query.lower() == 'exit':
                    print("👋 Adiós!\n")
                    break
                if query:
                    test_via_api(query)
            except KeyboardInterrupt:
                print("\n\n👋 Adiós!\n")
                break
    else:
        # Ejecutar query específica
        query = " ".join(sys.argv[1:])
        test_via_api(query)


if __name__ == "__main__":
    main()
