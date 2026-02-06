#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pyner_search v0.1 - Generador de términos de búsqueda para NCBI SRA
=================================================================
Propósito: Crear términos de búsqueda booleanos optimizados para NCBI SRA
Uso: python3 Pyner_search_v0.1.py <palabras_clave>

Ejemplo de uso:
    python3 Pyner_search_v0.1.py arabidopsis drought
    python3 Pyner_search_v0.1.py solanum lycopersicum nitrogen

El script:
1. Recibe palabras clave del usuario como argumentos de línea de comandos
2. Envía esas palabras a un LLM local que:
   - Entiende la sintaxis de búsqueda de NCBI SRA
   - Genera sinónimos biológicos relevantes
   - Construye una consulta booleana óptima
3. Devuelve una frase natural y una consulta para usar en NCBI

Output: JSON con dos campos:
- natural_query: descripción legible de la búsqueda
- esearch_query: término booleano para copiar en Pyner_v0.2.py o NCBI directamente
"""

import sys          # Para leer argumentos de línea de comandos
import json         # Para parsear JSON del LLM
import ollama       # Para usar modelo local

# ============================================
# 1️⃣ VALIDAR ENTRADA Y PROCESAR KEYWORDS
# ============================================

# Verificar que el usuario proporcionó palabras clave
if len(sys.argv) < 2:
    print(" Uso: python3 Pyner_search_v0.1.py <keywords separados por espacio>")
    print(" Ejemplo: python3 Pyner_search_v0.1.py arabidopsis drought")
    sys.exit(1)

# Obtener todas las palabras clave proporcionadas
keywords = sys.argv[1:]
print("\n Generando término de búsqueda a partir de keywords:")
print("   ", ", ".join(keywords))
print("--------------------------------------------------")

# ============================================
# 2️⃣ INICIALIZAR CLIENTE DE OLLAMA
# ============================================

# Cliente para llamar al modelo local (debe estar corriendo localmente)
client = ollama.Client()

# ============================================
# 3️⃣ CONSTRUIR PROMPT PARA EL LLM
# ============================================
# El prompt es muy importante: instrucciones claras = mejores resultados

prompt = f"""
Eres un asistente experto en bioinformática especializado en el uso de NCBI SRA (Sequence Read Archive).

Tu tarea es generar un término de búsqueda booleano válido para NCBI SRA usando E-utilities.

IMPORTANTE - RESTRICCIONES DE SINTAXIS:
- No uses campos como [substance], [mesh], ni otros específicos de PubMed.
- Usa solo campos que existen en SRA:
  [Organism], [All Fields], [Title], [Strategy], [Selection], [Platform], [Source], [Study]
- Si una palabra no encaja en un campo específico, usa [All Fields].
- Devuelve la búsqueda en formato booleano con AND y OR según corresponda.
- No incluyas explicaciones, solo el JSON.

Las palabras clave objetivo son:
{', '.join(keywords)}

Devuelve SOLO un JSON con este formato exacto:
{{
  "natural_query": "frase natural que describe la búsqueda",
  "esearch_query": "consulta compatible con NCBI SRA"
}}

Considera incluir al menos 3 sinónimos relevantes y términos relacionados para maximizar la recuperación de datos.
Si no existen sinónimos relevantes, no los fuerces.
"""

# ============================================
# 4️⃣ LLAMAR AL MODELO LOCAL
# ============================================

print("\n Llamando al LLM (qwen2.5:14b - local)...\n")

try:
    # Llamar al modelo con baja temperatura para respuestas más consistentes
    response = client.generate(
        model='qwen2.5:14b',
        prompt=prompt,
        options={'temperature': 0.2}  # 0.2 = determinístico pero flexible
    )
    # Extraer la respuesta del diccionario
    respuesta_llm = response['response'].strip()
except Exception as e:
    print(f" Error al llamar al modelo: {e}")
    print(" Asegúrate de que ollama está corriendo: ollama serve")
    sys.exit(1)

# ============================================
# 5️⃣ PARSEAR LA RESPUESTA JSON
# ============================================

print(" Respuesta cruda del modelo:\n")
print(respuesta_llm)
print("\n--------------------------------------------------")

# Intentar extraer JSON de la respuesta
try:
    data = json.loads(respuesta_llm)
    natural_query = data.get("natural_query", "N/A")
    esearch_query = data.get("esearch_query", "N/A")
except json.JSONDecodeError:
    # Si el LLM no devuelve JSON válido, usar la respuesta como está
    print(" Advertencia: El modelo no devolvió JSON válido")
    natural_query = respuesta_llm
    esearch_query = "No disponible"

# ============================================
# 6️⃣ MOSTRAR RESULTADOS
# ============================================

print(" Término de búsqueda generado:\n")
print(f" Natural query : {natural_query}")
print(f" ESearch query : {esearch_query}\n")

print("--------------------------------------------------")
print(" 📋 Instrucciones:")
print(f"    1. Copia el valor de 'esearch_query'")
print(f"    2. Pégalo en la variable 'search_term' de Pyner_v0.2.py")
print(f"    3. Ejecuta: python3 Pyner_v0.2.py")
print("--------------------------------------------------")
