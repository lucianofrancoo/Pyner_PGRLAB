#!/usr/bin/env python3
"""
Pyner v0.1 - Script ultra simple para buscar en NCBI
Busca: arabidopsis + sequía
"""

import time
from Bio import Entrez
import ollama
import json
import pandas as pd

# ============================================
# CONFIGURACIÓN BÁSICA
# ============================================

# Tu email es OBLIGATORIO para NCBI
Entrez.email = ""  # CAMBIAR ESTO!

# Tu API key (opcional pero recomendado)
Entrez.api_key = ""  # CAMBIAR ESTO!

# Cliente de Ollama
client = ollama.Client()

# ============================================
# PASO 1: BUSCAR EN NCBI
# ============================================

print("🔍 Buscando en NCBI: arabidopsis + sequía...")

# Construir la búsqueda - empezamos con GEO que tiene mucha metadata
search_term = 'arabidopsis[Organism] AND (drought OR "water stress" OR "water deficit" OR dehydration OR "sequía")'

# Buscar en GEO primero
handle = Entrez.esearch(
    db="gds",  # GEO DataSets
    term=search_term,
    retmax=10  # Solo 10 resultados para este test
)

# Leer los resultados
search_results = Entrez.read(handle)
handle.close()

print(f"✅ Encontrados {search_results['Count']} estudios totales")
print(f"📊 Procesando los primeros {len(search_results['IdList'])} estudios...\n")

# ============================================
# PASO 2: VER EL PRIMER ESTUDIO
# ============================================

print("="*50)
print("VEAMOS EL PRIMER ESTUDIO EN DETALLE:")
print("="*50)

# Tomar solo el primer ID
primer_id = search_results['IdList'][0]
print(f"\nID del estudio: {primer_id}")

# Obtener detalles del primer estudio
handle = Entrez.esummary(db="gds", id=primer_id)
estudio = Entrez.read(handle)[0]
handle.close()

# Imprimir TODOS los campos que vienen
print("\nCAMPOS DISPONIBLES:")
print("-"*30)
for campo, valor in estudio.items():
    print(f"{campo}: {valor}")

print("\n" + "="*50)
print("AHORA CON EL LLM:")
print("="*50)

# Preparar texto para el LLM
texto_para_llm = f"""
Title: {estudio.get('title', '')}
Summary: {estudio.get('summary', '')}
"""

print(f"\nTexto que le pasaremos al LLM:\n{texto_para_llm[:500]}...")

# Prompt simple para el LLM
prompt = f"""
Analyze this biological study and extract information.

{texto_para_llm}

Return a JSON with these fields:
- experimental_conditions: list all treatments or conditions studied
- is_time_series: true or false
- tissues_studied: list of tissues or cell types

Only return the JSON, no explanations.
"""

print("\nLlamando al LLM...")
response = client.generate(
    model='qwen2.5:14b',
    prompt=prompt,
    options={'temperature': 0.1}
)

print("\nRespuesta del LLM:")
print(response['response'])
