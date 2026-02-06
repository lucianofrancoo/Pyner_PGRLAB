#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pyner_search_v0.2.py
--------------------
Generador de términos de búsqueda optimizados para NCBI SRA (no PubMed ni PubChem).
Ejemplo:
    python3 Pyner_search_v0.2.py solanum lycopersicum nitrogen
"""

import sys
import json
import ollama

# ============================================
# 1️⃣ ENTRADA DE KEYWORDS
# ============================================

if len(sys.argv) < 2:
    print(" Uso: python3 Pyner_search_v0.2.py <keywords separados por espacio>")
    sys.exit(1)

keywords = sys.argv[1:]
print("\n Generando término de búsqueda a partir de keywords:")
print("   ", ", ".join(keywords))
print("--------------------------------------------------")

# ============================================
# 2️⃣ CLIENTE DE OLLAMA
# ============================================

client = ollama.Client()

# ============================================
# 3️⃣ PROMPT ACTUALIZADO (ENFOQUE SRA)
# ============================================

prompt = f"""
Eres un asistente experto en bioinformática especializado en el uso de NCBI SRA (Sequence Read Archive).

Tu tarea es generar un término de búsqueda booleano válido para NCBI SRA usando E-utilities.

IMPORTANTE:
- No uses campos como [substance], [mesh], ni otros específicos de PubMed.
- Usa solo campos que existen en SRA, por ejemplo:
  [Organism], [All Fields], [Title], [Strategy], [Selection], [Platform], [Source], [Study]
- Si una palabra no encaja en un campo específico, usa [All Fields].
- Devuelve la búsqueda en formato booleano con AND y OR según corresponda.
- No incluyas explicaciones, solo el JSON.

Las palabras clave son:
{', '.join(keywords)}

Devuelve un JSON con este formato:
{{
  "natural_query": "frase natural que describe la búsqueda",
  "esearch_query": "consulta compatible con NCBI SRA"
}}
Considera incluir al menos 3 sinónimos relevantes y términos relacionados para maximizar la recuperación de datos, si no existen entonces no.
"""

# ============================================
# 4️⃣ LLAMADA AL MODELO LOCAL
# ============================================

print("\n Llamando al LLM (qwen2.5:14b)...\n")

try:
    response = client.generate(
        model='qwen2.5:14b',
        prompt=prompt,
        options={'temperature': 0.2}
    )
    respuesta_llm = response['response'].strip()
except Exception as e:
    print(f" Error al llamar al modelo: {e}")
    sys.exit(1)

# ============================================
# 5️⃣ INTERPRETAR LA RESPUESTA
# ============================================

print(" Respuesta cruda del modelo:\n")
print(respuesta_llm)
print("\n--------------------------------------------------")

try:
    data = json.loads(respuesta_llm)
    natural_query = data.get("natural_query", "N/A")
    esearch_query = data.get("esearch_query", "N/A")
except json.JSONDecodeError:
    print(" El modelo no devolvió un JSON válido.")
    natural_query = respuesta_llm
    esearch_query = "No disponible"

# ============================================
# 6️⃣ MOSTRAR RESULTADOS
# ============================================

print(" Término de búsqueda generado:\n")
print(f" Natural query : {natural_query}")
print(f" ESearch query : {esearch_query}\n")

print("--------------------------------------------------")
print(" Copia el valor de 'esearch_query' para usarlo directamente en Pyner_V0.2.py")
print("--------------------------------------------------")
