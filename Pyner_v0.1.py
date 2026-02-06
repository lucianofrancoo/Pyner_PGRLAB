#!/usr/bin/env python3
"""
Pyner v0.1 - Script simple para buscar en NCBI GEO Dataset
================================================
Propósito: Buscar estudios de expresión génica en arabidopsis relacionados con sequía
Fuente: NCBI Gene Expression Omnibus (GEO) - Base de datos de microarrays y RNA-seq curados

Flujo del script:
1. Conecta a NCBI usando API de Entrez (BioPython)
2. Busca datasets de arabidopsis con palabras clave relacionadas a sequía
3. Extrae información del primer resultado
4. Envía los datos a un LLM local (ollama) para análisis automático
5. Imprime el análisis (condiciones, series temporales, tejidos)
"""

import time
from Bio import Entrez  # Para conectarse a NCBI
import ollama          # Para usar modelo local (qwen2.5)
import json
import pandas as pd

# ============================================
# CONFIGURACIÓN BÁSICA
# ============================================

# Tu email es OBLIGATORIO para NCBI (ética: NCBI quiere saber quién eres)
Entrez.email = ""  # CAMBIAR ESTO!

# Tu API key (opcional pero recomendado para aumentar límites de request)
Entrez.api_key = ""  # CAMBIAR ESTO!

# Cliente de Ollama para análisis con LLM
client = ollama.Client()

# ============================================
# PASO 1: BUSCAR EN NCBI
# ============================================

print("🔍 Buscando en NCBI: arabidopsis + sequía...")

# Construir la búsqueda usando sintaxis de NCBI
# [Organism] = campo específico para organismo
# AND = los resultados deben cumplir todas las condiciones
# OR = sinónimos de "sequía" para mayor cobertura
search_term = 'arabidopsis[Organism] AND (drought OR "water stress" OR "water deficit" OR dehydration OR "sequía")'

# Realizar búsqueda en GEO (Gene Expression Omnibus - estudios curados de expresión génica)
handle = Entrez.esearch(
    db="gds",                    # Base de datos: GEO DataSets
    term=search_term,
    retmax=10                    # Limitar a 10 resultados para pruebas
)

# Procesar respuesta XML de NCBI
search_results = Entrez.read(handle)
handle.close()

# Mostrar estadísticas
print(f"✅ Encontrados {search_results['Count']} estudios totales")
print(f"📊 Procesando los primeros {len(search_results['IdList'])} estudios...\n")

# ============================================
# PASO 2: EXTRAER DETALLES DEL PRIMER ESTUDIO
# ============================================

print("="*50)
print("DETALLES DEL PRIMER ESTUDIO:")
print("="*50)

# Tomar el primer ID del resultado
primer_id = search_results['IdList'][0]
print(f"\nID del estudio: {primer_id}")

# Obtener resumen detallado del primer estudio
handle = Entrez.esummary(db="gds", id=primer_id)
estudio = Entrez.read(handle)[0]
handle.close()

# Mostrar todos los campos disponibles (para inspeccionar estructura)
print("\nCAMPOS DISPONIBLES:")
print("-"*30)
for campo, valor in estudio.items():
    print(f"{campo}: {valor}")

print("\n" + "="*50)
print("ANÁLISIS CON LLM LOCAL:")
print("="*50)

# Preparar datos del estudio para enviar al LLM
# Extraemos título y resumen (los campos más importantes)
texto_para_llm = f"""
Title: {estudio.get('title', '')}
Summary: {estudio.get('summary', '')}
"""

print(f"\nDatos enviados al LLM:\n{texto_para_llm[:500]}...")

# Crear prompt para que el LLM extraiga información estructurada
prompt = f"""
Analyze this biological study and extract information.

{texto_para_llm}

Return a JSON with these fields:
- experimental_conditions: list all treatments or conditions studied
- is_time_series: true or false
- tissues_studied: list of tissues or cell types

Only return the JSON, no explanations.
"""

print("\nLlamando al LLM (qwen2.5:14b)...")
response = client.generate(
    model='qwen2.5:14b',
    prompt=prompt,
    options={'temperature': 0.1}  # Baja temperatura = respuestas más determinísticas
)

print("\nRespuesta del LLM:")
print(response['response'])
