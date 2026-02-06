#!/usr/bin/env python3
"""
Pyner v0.3 - Búsqueda automática en NCBI SRA (sin repetir BioProjects)
"""

import time
from Bio import Entrez
import ollama
import json
import pandas as pd
import xml.etree.ElementTree as ET

# ============================================
# CONFIGURACIÓN BÁSICA
# ============================================

Entrez.email = "lucianofranco.a@gmail.com"
Entrez.api_key = "4579ad4ab4c2144aa84a340b2e4374111308"

client = ollama.Client()

# ============================================
# PASO 1: BUSCAR EN NCBI SRA
# ============================================

print(" Buscando en NCBI SRA: arabidopsis + sequía...")

search_term = 'arabidopsis[Organism] AND (drought OR "water stress" OR "water deficit" OR dehydration OR "sequía")'

handle = Entrez.esearch(
    db="sra",
    term=search_term,
    retmax=2000  # procesaremos hasta 20 registros para ver variedad
)
search_results = Entrez.read(handle)
handle.close()

print(f" Encontrados {search_results['Count']} estudios totales")
print(f" Procesando los primeros {len(search_results['IdList'])} estudios...\n")

# ============================================
# PASO 2: RECORRER ESTUDIOS SIN REPETIR BIOPROJECT
# ============================================

bioprojects_vistos = set()
resultados = []

for idx, estudio_id in enumerate(search_results['IdList']):
    print("\n" + "=" * 60)
    print(f" Estudio {idx+1} | ID: {estudio_id}")
    print("=" * 60)

    handle = Entrez.esummary(db="sra", id=estudio_id, rettype="docsum")
    estudio = Entrez.read(handle)[0]
    handle.close()

    if "ExpXml" not in estudio:
        print(" Este registro no tiene ExpXml, se omite.")
        continue

    exp_xml = estudio["ExpXml"]
    wrapped_xml = f"<Root>{exp_xml}</Root>"

    try:
        root = ET.fromstring(wrapped_xml)
    except ET.ParseError as e:
        print(f" Error al parsear XML: {e}")
        continue

    bioproject = root.findtext(".//Bioproject", "")

    # 🚫 Evitar repetir BioProjects
    if bioproject in bioprojects_vistos:
        print(f" BioProject {bioproject} ya procesado, se salta.")
        continue
    bioprojects_vistos.add(bioproject)

    # Extraer campos
    title = root.findtext(".//Title", "")
    study_elem = root.find(".//Study")
    study_name = study_elem.attrib.get("name", "") if study_elem is not None else ""
    organism_elem = root.find(".//Organism")
    organism = organism_elem.attrib.get("ScientificName", "") if organism_elem is not None else ""
    library_strategy = root.findtext(".//Library_descriptor/LIBRARY_STRATEGY", "")
    library_source = root.findtext(".//Library_descriptor/LIBRARY_SOURCE", "")
    library_selection = root.findtext(".//Library_descriptor/LIBRARY_SELECTION", "")
    library_layout = (
        "PAIRED"
        if root.find(".//Library_descriptor/LIBRARY_LAYOUT/PAIRED") is not None
        else "SINGLE"
    )
    tissue = root.findtext(".//LIBRARY_NAME", "")
    biosample = root.findtext(".//Biosample", "")

    # Mostrar resumen
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
    # PASO 3: ANÁLISIS CON LLM
    # ============================================

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

    prompt = f"""
Analyze this biological sequencing study and extract information.

{texto_para_llm}

Return a JSON with these fields:
- experimental_conditions: list all treatments or conditions studied
- is_time_series: true or false
- tissues_studied: list of tissues or cell types

Only return the JSON, no explanations.
"""

    try:
        response = client.generate(
            model='qwen2.5:14b',
            prompt=prompt,
            options={'temperature': 0.1}
        )
        result_json = json.loads(response['response'])
    except Exception as e:
        print(f" Error en LLM: {e}")
        result_json = {"error": str(e)}

    # Guardar resultados
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
# GUARDAR RESULTADOS
# ============================================

if resultados:
    df = pd.DataFrame(resultados)
    df.to_csv("Pyner_SRA_arabidopsis_drought_unique.csv", index=False)
    print(f"\n Resultados guardados en Pyner_SRA_arabidopsis_drought_unique.csv")
else:
    print("\n No se generaron resultados.")
