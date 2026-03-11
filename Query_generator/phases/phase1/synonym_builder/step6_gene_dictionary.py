#!/usr/bin/env python3
"""
Paso 6: Gene Dictionary Builder
===============================

Objetivo:
  Construir un diccionario unificado de genes combinando 3 fuentes oficiales:
  1. NCBI Gene (gene_info.gz)
  2. UniProtKB (uniprot_sprot.dat.gz o idmapping)
  3. Ensembl (aliases)

  El diccionario final servirá para expandir términos génicos (ej. "BRCA1")
  a todos sus alias oficiales, nombres de proteínas y símbolos sin depender de MeSH.

NOTA: Este script asume que los archivos crudos ya han sido descargados en 
el directorio 'gene_data'. Dado que pesan decenas de GB en total, este script
procesará un archivo a la vez construyendo un índice estructurado.
"""

import os
import gzip
import json
import time
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "gene_data"
OUTPUT_FILE = BASE_DIR / "output" / "final_gene_dictionary.json"

NCBI_FILE = DATA_DIR / "gene_info.gz"
UNIPROT_FILE = DATA_DIR / "uniprot_sprot.dat.gz"
# Dependiendo de qué archivo de Ensembl se use (ej. un mart export o el homo_sapiens.xyz)

def process_ncbi_gene(file_path):
    """
    Parsea gene_info.gz de NCBI (separado por tabuladores)
    Columnas principales:
    0: tax_id
    1: GeneID
    2: Symbol
    3: LocusTag
    4: Synonyms
    5: dbXrefs
    6: chromosome
    7: map_location
    8: description
    """
    print(f"Procesando NCBI Gene: {file_path.name}")
    genes = {}
    
    if not file_path.exists():
        print(f"⚠ Archivo no encontrado: {file_path}")
        return genes

    with gzip.open(file_path, 'rt', encoding='utf-8') as f:
        next(f, None) # saltar header
        
        count = 0
        for line in f:
            count += 1
            if count % 5_000_000 == 0:
                print(f"  ... {count:,} líneas procesadas")
                
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 9:
                tax_id = parts[0]
                gene_id = parts[1]
                symbol = parts[2]
                synonyms = parts[4]
                description = parts[8]
                
                if symbol == '-':
                    continue
                    
                syn_list = [s.strip() for s in synonyms.split('|')] if synonyms != '-' else []
                syn_list = [s for s in syn_list if s]

                # Nueva Clave: tax_id + gene_id para evitar colapsar organismos distintos
                key = f"{tax_id}_{gene_id}"
                
                genes[key] = {
                    "tax_id": tax_id,
                    "gene_id": gene_id,
                    "preferred_symbol": symbol,
                    "description": description if description != '-' else "",
                    "synonyms": set(syn_list),
                    "sources": {"ncbi_gene"}
                }

    print(f"✓ NCBI Gene procesado: {len(genes):,} genes únicos extraídos")
    return genes

def process_uniprot(file_path, genes_dict):
    """
    Parsea uniprot_sprot.dat.gz (Swiss-Prot).
    Extrae Gene Names (GN), Protein Names (DE) y dbXrefs de GeneID NCBI.
    Cruza por GeneID cuando es posible, de lo contrario guarda huérfanos.
    """
    print(f"Procesando UniProtKB (Swiss-Prot): {file_path.name}")
    if not file_path.exists():
        print(f"⚠ Archivo no encontrado: {file_path}")
        return genes_dict
        
    current_gene_names = []
    current_protein_names = []
    current_taxid = ""
    current_geneid = ""
    count = 0
    
    with gzip.open(file_path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith("ID "):
                count += 1
                current_gene_names = []
                current_protein_names = []
                current_taxid = ""
                current_geneid = ""
                if count % 100_000 == 0:
                    print(f"  ... {count:,} entradas procesadas")
                    
            elif line.startswith("OX   NCBI_TaxID="):
                # OX   NCBI_TaxID=9606;
                current_taxid = line.split("NCBI_TaxID=")[1].split(";")[0].strip()
                
            elif line.startswith("DR   GeneID;"):
                # DR   GeneID; 10458; -
                current_geneid = line.split(";")[1].strip()
                    
            elif line.startswith("DE "):
                if "RecName: Full=" in line or "AltName: Full=" in line:
                    pname = line.split("Full=")[1].split(";")[0].strip()
                    if "{" in pname: pname = pname.split("{")[0].strip()
                    current_protein_names.append(pname)
                    
            elif line.startswith("GN "):
                parts = line[5:].split(";")
                for p in parts:
                    p = p.strip()
                    if p.startswith("Name="):
                        gname = p[5:]
                        if "{" in gname: gname = gname.split("{")[0].strip()
                        current_gene_names.append(gname)
                    elif p.startswith("Synonyms="):
                        syns = p[9:].split(",")
                        for s in syns:
                            s = s.strip()
                            if "{" in s: s = s.split("{")[0].strip()
                            current_gene_names.append(s)
                            
            elif line.startswith("//"): 
                if not current_gene_names:
                    continue
                    
                all_aliases = set(current_gene_names[1:] + current_protein_names)
                
                # Cross-reference with NCBI if possible
                key = f"{current_taxid}_{current_geneid}"
                
                if current_geneid and key in genes_dict:
                    genes_dict[key]["synonyms"].update(all_aliases)
                    genes_dict[key]["sources"].add("uniprot")
                    if not genes_dict[key]["description"] and current_protein_names:
                        genes_dict[key]["description"] = current_protein_names[0]
                else: # Uniprot entry with no NCBI GeneID mapped
                    orph_key = f"uniprot_{count}"
                    genes_dict[orph_key] = {
                        "tax_id": current_taxid,
                        "gene_id": "",
                        "preferred_symbol": current_gene_names[0],
                        "description": current_protein_names[0] if current_protein_names else "",
                        "synonyms": all_aliases,
                        "sources": {"uniprot"}
                    }
                    
    print(f"✓ UniProt procesado. Proteínas totales procesadas aportadas a genes/huérfanos.")
    return genes_dict

def merge_and_save(genes_dict, output_path):
    print("Consolidando diccionario final indexado por Símbolo/Alias...")
    final_dict = {}
    
    # 1. Limpieza de set -> list
    genes_list = list(genes_dict.values())
    
    # 2. Generar un diccionario invertido para búsqueda ultra-rápida (Symbol -> List[Gene_Data])
    for gene in genes_list:
        p_symbol = gene["preferred_symbol"]
        p_lower = p_symbol.lower()
        
        # Limpieza de sinónimos redundantes
        clean_syns = []
        for s in gene["synonyms"]:
            if s.lower() != p_lower and len(s) > 2:
                clean_syns.append(s)
                
        # Estructura limpia (sin redundancias inútiles)
        gene_obj = {
            "preferred_symbol": p_symbol,
            "tax_id": gene.get("tax_id"),
            "description": gene["description"],
            "synonyms": list(set(clean_syns)),
            "sources": list(gene["sources"])
        }
        
        # Indexar bajo el símbolo principal
        if p_lower not in final_dict:
            final_dict[p_lower] = {
                "preferred_symbol": p_symbol,
                "synonyms": list(set(clean_syns)),
                "sources": list(gene["sources"])
            }
            # No guardamos las variantes taxonómicas múltiples si la query es ciega.
            # Simplemente guardamos un "mega-objeto" agnóstico de taxID, sumando info.
        else:
            final_dict[p_lower]["synonyms"].extend(clean_syns)
            final_dict[p_lower]["synonyms"] = list(set(final_dict[p_lower]["synonyms"]))
            final_dict[p_lower]["sources"].extend(gene["sources"])
            final_dict[p_lower]["sources"] = list(set(final_dict[p_lower]["sources"]))
            
    # Filtro opcional: Para reducir el tamaño del diccionario, podríamos omitir aquellos keys 
    # que sean solo combinaciones excesivamente genéricas, pero por ahora conservamos todos.

    final_dict_size = len(final_dict)
    print(f"Guardando {final_dict_size:,} símbolos primarios globales en {output_path.name}...")
    
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(final_dict, f, separators=(',', ':'))
    print("✅ Diccionario completo guardado de forma compacta.")

def main():
    print("=" * 60)
    print(" PASO 6: CONSTRUCCIÓN DEL DICCIONARIO UNIFICADO DE GENES ")
    print("=" * 60)
    
    DATA_DIR.mkdir(exist_ok=True)
    
    if not NCBI_FILE.exists() and not UNIPROT_FILE.exists():
        print(f"⚠ Faltan los archivos crudos en: {DATA_DIR}")
        print("Debes descargar:")
        print(" 1) curl -O https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz")
        print(" 2) curl -O https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz")
        return
        
    t0 = time.time()
    
    genes = {}
    genes = process_ncbi_gene(NCBI_FILE)
    genes = process_uniprot(UNIPROT_FILE, genes)
    
    # Ensembl processing can be added here once we define the specific export format
    # (Ensembl doesn't have a single flat file for all species like NCBI/UniProt, 
    # usually requires BioMart lists per species, e.g., using a text TSV dump).
    
    if genes:
        merge_and_save(genes, OUTPUT_FILE)
        
    print(f"Completado en {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()
