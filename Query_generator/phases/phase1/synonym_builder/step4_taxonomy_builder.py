#!/usr/bin/env python3
"""
Paso 4: Taxonomy Builder (Opcional - Recomendado)
=================================================

Objetivo:
  Integrar el dump de taxonomía del NCBI (names.dmp) a nuestro `final_synonym_dictionary.json`.
  La taxonomía oficial de NCBI elimina el sesgo clínico que tiene MeSH para organismos
  (ej. evitar que PubMedBERT asocie "wheat" con "Wheat Hypersensitivity").

Entrada:
  - taxdump/names.dmp (separado por \t|\t)
  - output/final_synonym_dictionary.json (el diccionario que ya contiene MeSH + variantes)

Salida:
  - output/final_synonym_dictionary.json (actualizado in-place con organismos)

Lógica:
  1. Lee names.dmp y agrupa por tax_id.
  2. Identifica el "scientific name" como preferred_term.
  3. Convierte "genbank common name", "common name" y "synonym" en keywords/alias.
  4. Agrega o actualiza cada organismo en el diccionario principal.
"""

import os
import json
import time
from pathlib import Path

# =====================================================================
# CONFIGURACIÓN
# =====================================================================
BASE_DIR    = Path(__file__).resolve().parent
NAMES_DMP   = BASE_DIR / "taxdump" / "names.dmp"
DICT_PATH   = BASE_DIR / "output" / "final_synonym_dictionary.json"

# Nombres a excluir por ser muy inespecíficos o problemáticos
EXCLUDE_NAMES = {
    "all", "root", "environmental samples", "unidentified", "unclassified",
    "unknown", "other", "synthetic construct", "artificial sequences"
}

def clean_name(name: str) -> str:
    """Limpia un nombre eliminando tags raros y bajando a minúsculas para la key."""
    # Remover tags como <bacteria> al final
    if '<' in name and name.endswith('>'):
        name = name[:name.rfind('<')].strip()
    return name

def process_names_dmp(file_path: Path) -> dict:
    """
    Lee names.dmp y construye un diccionario:
    {
      tax_id: {
        "scientific_name": "...",
        "synonyms": [{"term": "...", "type": "..."}, ...]
      }
    }
    """
    print(f"Leyendo taxonomía desde {file_path.name}...")
    
    if not file_path.exists():
        raise FileNotFoundError(f"No se encontró {file_path}. Debes descargar taxdump.tar.gz primero.")

    taxa = {}
    line_count = 0
    
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_count += 1
            if line_count % 1_000_000 == 0:
                print(f"  ...procesadas {line_count:,} líneas")
                
            parts = [p.strip() for p in line.split('\t|\t')]
            if len(parts) >= 4:
                tax_id = parts[0]
                name = parts[1]
                # parts[2] = unique name (no siempre usado)
                name_class = parts[3].replace('\t|', '').strip()
                
                clean_n = clean_name(name)
                if not clean_n or clean_n.lower() in EXCLUDE_NAMES:
                    continue
                
                if tax_id not in taxa:
                    taxa[tax_id] = {"scientific_name": None, "synonyms": []}
                
                # Clasificar el nombre según su name_class
                if name_class == "scientific name":
                    taxa[tax_id]["scientific_name"] = clean_n
                elif name_class in ["synonym", "equivalent name"]:
                    taxa[tax_id]["synonyms"].append({"term": clean_n, "type": "scientific_synonym"})
                elif name_class in ["common name", "genbank common name"]:
                    taxa[tax_id]["synonyms"].append({"term": clean_n, "type": "common_name"})
                elif name_class == "blast name":
                    taxa[tax_id]["synonyms"].append({"term": clean_n, "type": "common_name"})
                # Ignoramos "authority", "type material", "in-part", etc. para evitar ruido

    print(f"✓ {len(taxa):,} taxones extraídos en total")
    return taxa

def merge_taxonomy_into_dict(taxa: dict, dict_path: Path):
    """
    Toma los taxones extraídos y los inyecta en el final_synonym_dictionary.json.
    """
    print(f"Cargando diccionario existente desde {dict_path.name}...")
    t0 = time.time()
    
    with open(dict_path, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    dictionary = data["dictionary"]
    original_size = len(dictionary)
    print(f"✓ Diccionario original tenía {original_size:,} entradas")
    
    # Procesar taxones
    added_count = 0
    updated_count = 0
    
    # Solo procesaremos organismos que tengan algún nombre común o alias útil
    # Si es solo el scientific name estricto sin variantes conocidas, NCBI de todos modos lo entiende.
    # El valor real está en los common names ("wheat", "mouse", "human").
    
    useful_taxa = {
        tid: info for tid, info in taxa.items() 
        if info["scientific_name"] and info["synonyms"]
    }
    print(f"Filtrados {len(useful_taxa):,} taxones que tienen sinónimos o nombres comunes útiles")

    for tax_id, info in useful_taxa.items():
        scientific = info["scientific_name"]
        
        # Las "keys" (entradas directas en el dict) serán el nombre científico 
        # y todos sus nombres comunes.
        
        # 1. Preparar el bloque de sinónimos taxonómicos
        tax_synonyms = []
        for syn in info["synonyms"]:
            tax_synonyms.append({
                "term": syn["term"],
                "source": "ncbi_taxonomy",
                "confidence": "official" if syn["type"] == "common_name" else 0.95
            })
            
        keys_to_populate = [scientific.lower().strip()] + [ s["term"].lower().strip() for s in info["synonyms"] ]
        keys_to_populate = list(set([k for k in keys_to_populate if len(k) > 2]))
        
        for key in keys_to_populate:
            if key in dictionary:
                # Actualizar entrada existente: Si ya existía por MeSH, le agregamos los sinónimos taxonómicos
                entry = dictionary[key]
                existing_terms = { s["term"].lower() for s in entry.get("synonyms", []) }
                
                # Agregar scientific name si no estaba
                if scientific.lower() != entry["preferred_term"].lower() and scientific.lower() not in existing_terms:
                    entry["synonyms"].append({
                        "term": scientific,
                        "source": "ncbi_taxonomy",
                        "confidence": "official"
                    })
                    existing_terms.add(scientific.lower())
                
                # Agregar las otras variantes
                for tsyn in tax_synonyms:
                    if tsyn["term"].lower() not in existing_terms and tsyn["term"].lower() != entry["preferred_term"].lower():
                        entry["synonyms"].append(tsyn)
                        existing_terms.add(tsyn["term"].lower())
                
                updated_count += 1
            else:
                # Crear nueva entrada puramente taxonómica
                dictionary[key] = {
                    "preferred_term": scientific,
                    "tax_id": tax_id,
                    "synonyms": [s for s in tax_synonyms if s["term"].strip().lower() != scientific.strip().lower()],
                    "synonym_count": {
                        "total": len(tax_synonyms),
                        "mesh": 0,
                        "text": 0,
                        "taxonomy": len(tax_synonyms)
                    }
                }
                added_count += 1

    print(f"✓ Integración completa:")
    print(f"  • Entradas originales actualizadas (MeSH + Taxon): {updated_count:,}")
    print(f"  • Nuevas entradas agregadas (Solo Taxon):          {added_count:,}")
    print(f"  • Total entradas en el diccionario ahora:          {len(dictionary):,}")
    
    # Actualizar metadata
    data["metadata"]["sources"]["taxonomy"] = "ncbi_taxdump (names.dmp)"
    stats = data["metadata"].setdefault("statistics", {})
    stats["total_terms"] = len(dictionary)
    stats["taxonomy_added"] = added_count
    data["metadata"]["version"] = "3.0 (MeSH + Variants + Taxonomy)"
    
    # Guardar
    print(f"\nGuardando diccionario actualizado en {dict_path.name} (esto puede tomar 1 min)...")
    with open(dict_path, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(',', ':')) # compact size
        
    print(f"✅ ¡Diccionario de sinónimos actualizado! (Tardó {time.time()-t0:.1f}s)")

if __name__ == "__main__":
    print("=" * 60)
    print(" STEP 4: NCBI TAXONOMY BUILDER ")
    print("=" * 60)
    
    try:
        t0 = time.time()
        taxa = process_names_dmp(NAMES_DMP)
        merge_taxonomy_into_dict(taxa, DICT_PATH)
        print("=" * 60)
        print(f"Pipeline completado en {time.time()-t0:.1f} segundos.")
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
