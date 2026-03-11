#!/usr/bin/env python3
"""
Test Simple: LLM + Query Expander → NCBI Boolean Query

Flujo:
1. Usuario da query en lenguaje natural
2. LLM (Ollama) extrae keywords
3. query_expander.py expande sinónimos
4. Genera boolean query para NCBI

NO usa datos de SRA (Phase 1/2), solo el diccionario de sinónimos.
"""
import sys
import json
import requests
from pathlib import Path
from query_expander import QueryExpander

# Configuración Ollama
OLLAMA_URL = "http://localhost:11434/api/generate"
OLLAMA_MODEL = "qwen2.5:14b"  # Cambiar según tu modelo instalado

def call_ollama(prompt, model=OLLAMA_MODEL):
    """Llamar a Ollama API local"""
    try:
        response = requests.post(
            OLLAMA_URL,
            json={
                "model": model,
                "prompt": prompt,
                "stream": False
            },
            timeout=30
        )

        if response.status_code == 200:
            return response.json().get("response", "").strip()
        else:
            return None
    except Exception as e:
        print(f"❌ Error conectando con Ollama: {e}")
        print(f"   Asegúrate de que Ollama esté corriendo: ollama serve")
        return None

def extract_keywords_with_llm(user_query):
    """Extraer keywords científicos de la query del usuario con LLM"""

    prompt = f"""Eres un experto en biología molecular y genómica.

Analiza la siguiente consulta científica y extrae SOLO los términos clave científicos más importantes.

Consulta: "{user_query}"

Instrucciones:
- Extrae términos relacionados con: organismos, técnicas, condiciones, fenotipos, genes
- Devuelve SOLO una lista de términos separados por comas
- NO incluyas palabras vacías (de, en, sobre, etc.)
- Normaliza términos (ej: "sequencing de RNA" → "rna sequencing")
- Máximo 5 términos

Ejemplo:
Consulta: "Estudios de estrés por sequía en arabidopsis usando rna-seq"
Respuesta: drought stress, arabidopsis, rna-seq

Respuesta (solo términos separados por comas):"""

    print("🤖 Llamando a LLM para extraer keywords...")
    response = call_ollama(prompt)

    if not response:
        return None

    # Parsear respuesta
    keywords = [k.strip() for k in response.split(",")]
    keywords = [k for k in keywords if k]  # Remover vacíos

    return keywords

def generate_ncbi_query(user_query, use_llm=True):
    """
    Genera una query NCBI desde lenguaje natural

    Args:
        user_query: Query del usuario en lenguaje natural
        use_llm: Usar LLM para extraer keywords (True) o manual (False)

    Returns:
        Dict con keywords, expansiones y boolean query final
    """
    print("="*80)
    print("QUERY GENERATOR - LLM + Synonym Expander")
    print("="*80)
    print()

    # Paso 1: Extraer keywords
    print(f"📝 Query del usuario: \"{user_query}\"")
    print()

    if use_llm:
        keywords = extract_keywords_with_llm(user_query)
        if not keywords:
            print("❌ No se pudieron extraer keywords con LLM")
            return None
    else:
        # Modo manual: extraer palabras simples
        keywords = [w.strip() for w in user_query.lower().split()
                   if len(w.strip()) > 2]

    print(f"🔑 Keywords extraídos: {keywords}")
    print()

    # Paso 2: Inicializar Query Expander
    print("📚 Cargando diccionario de sinónimos...")
    try:
        expander = QueryExpander(
            use_semantic_fallback=True,
            semantic_threshold=0.80,
            max_fallback_results=3
        )
    except Exception as e:
        print(f"❌ Error cargando expander: {e}")
        return None

    print("✓ Diccionario cargado")
    print()

    # Paso 3: Expandir query
    print("🔍 Expandiendo sinónimos...")
    result = expander.expand_query(keywords)

    # Mostrar expansiones
    print()
    print("="*80)
    print("EXPANSIONES DE SINÓNIMOS")
    print("="*80)
    print()

    for keyword, data in result["expanded"].items():
        print(f"📌 '{keyword}':")
        print(f"   Método: {data.get('method', 'N/A')}")

        if data["found"]:
            print(f"   Preferred: {data['preferred_term']}")
            print(f"   Total sinónimos: {len(data['synonyms'])}")

            # Mostrar top 5 sinónimos
            for syn in data["synonyms"][:5]:
                source = syn.get("source", "N/A")
                conf = syn.get("confidence", "N/A")
                print(f"      • {syn['term']} ({source}, conf={conf})")

        elif data.get("suggestions"):
            print(f"   ⚠️ No encontrado en diccionario, usando fallback semántico:")
            for sug in data["suggestions"]:
                print(f"      • {sug['term']} (similarity={sug['similarity']})")
        else:
            print(f"   ⚠️ Sin sinónimos (se usará término original)")

        print()

    # Paso 4: Generar boolean query
    print("="*80)
    print("QUERY NCBI FINAL")
    print("="*80)
    print()
    print(result["boolean_query"])
    print()

    return result

def interactive_mode():
    """Modo interactivo para probar múltiples queries"""
    print("="*80)
    print("QUERY GENERATOR - MODO INTERACTIVO")
    print("="*80)
    print()
    print("Escribe queries en lenguaje natural para generar boolean queries NCBI")
    print("Comandos:")
    print("  - 'exit' o 'quit': Salir")
    print("  - 'manual': Cambiar a extracción manual de keywords")
    print("  - 'llm': Cambiar a extracción con LLM")
    print()

    use_llm = True

    while True:
        try:
            query = input("\n🔬 Query: ").strip()

            if not query:
                continue

            if query.lower() in ['exit', 'quit', 'q']:
                print("\n👋 ¡Hasta luego!")
                break

            if query.lower() == 'manual':
                use_llm = False
                print("✓ Modo manual activado")
                continue

            if query.lower() == 'llm':
                use_llm = True
                print("✓ Modo LLM activado")
                continue

            # Generar query
            result = generate_ncbi_query(query, use_llm=use_llm)

            if result:
                # Preguntar si quiere guardar
                save = input("\n💾 ¿Guardar esta query? (y/n): ").strip().lower()
                if save == 'y':
                    filename = f"query_{len(list(Path('.').glob('query_*.json')))}.json"
                    with open(filename, 'w') as f:
                        json.dump({
                            "user_query": query,
                            "keywords": result["original_keywords"],
                            "boolean_query": result["boolean_query"],
                            "expanded": result["expanded"]
                        }, f, indent=2)
                    print(f"✓ Guardado en {filename}")

        except KeyboardInterrupt:
            print("\n\n👋 ¡Hasta luego!")
            break
        except Exception as e:
            print(f"\n❌ Error: {e}")

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Test del Query Generator con LLM + Synonym Expander",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:

# Query única con LLM
python3 test_llm_query_generator.py "estudios de sequía en arabidopsis con rnaseq"

# Query única sin LLM (manual)
python3 test_llm_query_generator.py --no-llm "drought arabidopsis rnaseq"

# Modo interactivo
python3 test_llm_query_generator.py -i

Queries de ejemplo:
  - "estudios de estrés hídrico en plantas usando transcriptómica"
  - "drought tolerance in wheat using rna sequencing"
  - "arabidopsis heat stress proteomics"
  - "rice salt stress gene expression"
        """
    )

    parser.add_argument(
        "query",
        nargs="?",
        help="Query del usuario en lenguaje natural"
    )
    parser.add_argument(
        "-i", "--interactive",
        action="store_true",
        help="Modo interactivo"
    )
    parser.add_argument(
        "--no-llm",
        action="store_true",
        help="No usar LLM, extracción manual de keywords"
    )
    parser.add_argument(
        "--model",
        default=OLLAMA_MODEL,
        help=f"Modelo Ollama a usar (default: {OLLAMA_MODEL})"
    )

    args = parser.parse_args()

    # Actualizar modelo global
    global OLLAMA_MODEL
    OLLAMA_MODEL = args.model

    if args.interactive:
        interactive_mode()
    elif args.query:
        result = generate_ncbi_query(args.query, use_llm=not args.no_llm)
        if not result:
            sys.exit(1)
    else:
        parser.print_help()
        print("\n" + "="*80)
        print("QUICK START")
        print("="*80)
        print()
        print("1. Asegúrate de que Ollama esté corriendo:")
        print("   ollama serve")
        print()
        print("2. Ejecuta una query de prueba:")
        print('   python3 test_llm_query_generator.py "drought stress in arabidopsis"')
        print()
        print("3. O usa modo interactivo:")
        print("   python3 test_llm_query_generator.py -i")
        print()

if __name__ == "__main__":
    main()
