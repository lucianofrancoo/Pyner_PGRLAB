#!/usr/bin/env python3
"""
Pyner v0.2 - Búsqueda avanzada en NCBI SRA con deduplicación automática
==================================================================
Propósito: Extraer datos de secuenciación (RNA-seq, genomics) de arabidopsis con sequía
Fuente: NCBI Sequence Read Archive (SRA) - Acceso público a datos crudos de secuenciación

Características principales:
1. Busca en NCBI SRA (millones de estudios de secuenciación)
2. Evita procesar datos duplicados del mismo BioProject
3. Extrae metadatos (estrategia de biblioteca, tejido, organismo, etc.)
4. Envía cada estudio a un LLM para análisis automático
5. Exporta resultados como CSV para análisis posterior

Flujo:
1. Buscar en SRA
2. Para cada resultado:
   - Parsear XML de metadatos
   - Extraer BioProject ID
   - Si es nuevo BioProject: procesar y analizar con LLM
   - Guardar en lista de resultados
3. Exportar resultados a CSV
"""

import time
from Bio import Entrez        # Para conectarse a NCBI
import ollama                 # Para usar modelo local
import json
import pandas as pd           # Para tabular resultados
import xml.etree.ElementTree as ET  # Para parsear respuestas XML de NCBI

# ============================================
# CONFIGURACIÓN BÁSICA - CREDENCIALES NCBI
# ============================================

# Credenciales para acceso a NCBI (obligatorias)
Entrez.email = "lucianofranco.a@gmail.com"
# API key aumenta límites: 10 req/seg sin key, 3 req/sec con key
Entrez.api_key = "4579ad4ab4c2144aa84a340b2e4374111308"

# Cliente para usar modelo local con ollama
client = ollama.Client()

# ============================================
# PASO 1: BÚSQUEDA INICIAL EN NCBI SRA
# ============================================

print(" Buscando en NCBI SRA: arabidopsis + sequía...")

# Términos de búsqueda booleana para SRA
search_term = 'arabidopsis[Organism] AND (drought OR "water stress" OR "water deficit" OR dehydration OR "sequía")'

# Realizar búsqueda en SRA
handle = Entrez.esearch(
    db="sra",              # Base de datos Sequence Read Archive
    term=search_term,
    retmax=2000            # Máximo de resultados a procesar
)
search_results = Entrez.read(handle)
handle.close()

# Mostrar estadísticas iniciales
print(f" Encontrados {search_results['Count']} estudios totales")
print(f" Procesando los primeros {len(search_results['IdList'])} estudios...\n")

# ============================================
# PASO 2: RECORRER ESTUDIOS Y EVITAR DUPLICADOS
# ============================================

# Set para almacenar BioProjects ya procesados (evita duplicados)
bioprojects_vistos = set()
# Lista para almacenar resultados finales
resultados = []

# Iterar sobre cada ID de estudio
for idx, estudio_id in enumerate(search_results['IdList']):
    print("\n" + "=" * 60)
    print(f" Estudio {idx+1} | ID: {estudio_id}")
    print("=" * 60)

    # Obtener resumen del estudio en formato XML
    handle = Entrez.esummary(db="sra", id=estudio_id, rettype="docsum")
    estudio = Entrez.read(handle)[0]
    handle.close()

    # Verificar si tiene datos XML (algunos registros pueden estar incompletos)
    if "ExpXml" not in estudio:
        print(" Este registro no tiene ExpXml, se omite.")
        continue

    # Extraer y preparar XML para parsing
    exp_xml = estudio["ExpXml"]
    # Envolver en elemento raíz para que sea XML válido
    wrapped_xml = f"<Root>{exp_xml}</Root>"

    # Parsear XML
    try:
        root = ET.fromstring(wrapped_xml)
    except ET.ParseError as e:
        print(f" Error al parsear XML: {e}")
        continue

    # Extraer BioProject ID (identificador único del estudio biológico)
    bioproject = root.findtext(".//Bioproject", "")

    # ============================================
    # DEDUPLICACIÓN: EVITAR PROCESAR EL MISMO BIOPROJECT
    # ============================================
    if bioproject in bioprojects_vistos:
        print(f" BioProject {bioproject} ya procesado, se salta.")
        continue
    # Marcar este BioProject como procesado
    bioprojects_vistos.add(bioproject)

    # ============================================
    # EXTRACCIÓN DE METADATOS DEL XML
    # ============================================
    
    # Información general del estudio
    title = root.findtext(".//Title", "")  # Título del experimento
    study_elem = root.find(".//Study")
    study_name = study_elem.attrib.get("name", "") if study_elem is not None else ""
    
    # Información del organismo
    organism_elem = root.find(".//Organism")
    organism = organism_elem.attrib.get("ScientificName", "") if organism_elem is not None else ""
    
    # Información sobre la librería de secuenciación
    library_strategy = root.findtext(".//Library_descriptor/LIBRARY_STRATEGY", "")  # Ej: RNA-Seq, WXS
    library_source = root.findtext(".//Library_descriptor/LIBRARY_SOURCE", "")      # Ej: TRANSCRIPTOMIC
    library_selection = root.findtext(".//Library_descriptor/LIBRARY_SELECTION", "") # Ej: cDNA
    
    # Tipo de secuenciación (single-end vs paired-end)
    library_layout = (
        "PAIRED"
        if root.find(".//Library_descriptor/LIBRARY_LAYOUT/PAIRED") is not None
        else "SINGLE"
    )
    
    # Tejido o célula de donde se extrajo la muestra
    tissue = root.findtext(".//LIBRARY_NAME", "")
    # ID de la muestra biológica
    biosample = root.findtext(".//Biosample", "")

    # Mostrar información extraída
    print(f"Title: {title}")
    print(f"Study: {study_name}")
    print(f"Organism: {organism}")
    print(f"Library strategy: {library_strategy}")
    print(f"Library source: {library_source}")
    print(f"Library selection: {library_selection}")
    print(f"Library layout: {library_layout}")
    print(f"Tissue (from library name): {tissue}")
    print(f"BioProject: {bioproject}")
    print(f"BioSample: {biosample}")

    # ============================================
    # PASO 3: ANÁLISIS CON LLM LOCAL
    # ============================================
    
    # Preparar texto con los metadatos para enviar al LLM
    texto_para_llm = f"""
Title: {title}
Study: {study_name}
Organism: {organism}
Library strategy: {library_strategy}
Library source: {library_source}
Library selection: {library_selection}
Library layout: {library_layout}
Tissue (from library name): {tissue}
BioProject: {bioproject}
BioSample: {biosample}
"""

    # Prompt para que el LLM extraiga información estructurada
    prompt = f"""
Analyze this biological sequencing study and extract information.

{texto_para_llm}

Return a JSON with these fields:
- experimental_conditions: list all treatments or conditions studied
- is_time_series: true or false
- tissues_studied: list of tissues or cell types

Only return the JSON, no explanations.
"""

    # Llamar al LLM localmente (ollama corriendo con qwen2.5:14b)
    try:
        response = client.generate(
            model='qwen2.5:14b',
            prompt=prompt,
            options={'temperature': 0.1}  # Determinístico para resultados consistentes
        )
        result_json = json.loads(response['response'])
    except Exception as e:
        print(f" Error en LLM: {e}")
        result_json = {"error": str(e)}

    # Guardar este resultado con sus análisis
    resultados.append({
        "bioproject": bioproject,
        "study": study_name,
        "organism": organism,
        "tissue": tissue,
        "conditions": result_json.get("experimental_conditions", []),
        "is_time_series": result_json.get("is_time_series", None),
        "tissues_studied": result_json.get("tissues_studied", []),
    })

# ============================================
# PASO 4: EXPORTAR RESULTADOS A CSV
# ============================================

if resultados:
    # Crear DataFrame de pandas con los resultados
    df = pd.DataFrame(resultados)
    # Guardar como CSV para análisis posterior en Excel, Python, etc.
    df.to_csv("Pyner_SRA_arabidopsis_drought_unique.csv", index=False)
    print(f"\n Resultados guardados en Pyner_SRA_arabidopsis_drought_unique.csv")
else:
    print("\n No se generaron resultados.")
