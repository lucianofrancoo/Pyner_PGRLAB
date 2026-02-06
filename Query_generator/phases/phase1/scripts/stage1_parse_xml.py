"""
Pyner Phase 1 - Stage 1: XML Parser
===================================
Parsear archivos XML de NCBI SRA (experiment, sample, run)

Stage 1 Overview:
1. Parsear experiment.xml (título, estrategia, librería)
2. Parsear sample.xml (organismo, atributos)
3. Parsear run.xml (información de ejecución)
4. Normalizar campos
5. Guardar en índices

Ejecución:
    python scripts/stage1_parse_xml.py --max-files 1000 --test-mode
    python scripts/stage1_parse_xml.py --max-files 50000 --checkpoint

"""

import json
import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime
import time
from collections import defaultdict
from typing import Dict, List, Any, Optional, Tuple
import sys

# Add parent directory to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import (
    NCBI_SRA_PATH, OUTPUT_DIR, MAX_FILES_PHASE1_STAGE1,
    CHECKPOINT_INTERVAL, DEBUG_PRINT_INTERVAL
)
from utils import (
    setup_logging, print_section_header, CheckpointManager,
    ProgressTracker, SystemMonitor, print_stage_summary
)


# ============================================
# CONSTANTS
# ============================================

EXPECTED_ORGANISMS = 4500  # Aprox. cuántos organismos esperamos
EXPECTED_STRATEGIES = 67   # Estrategias NCBI SRA estándar
BUFFER_SIZE = 100  # Números de archivos antes de guardar buffer


# ============================================
# DATA STRUCTURES
# ============================================

class XMLParser:
    """Parser robusto de XMLs de NCBI SRA"""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.errors = {
            "parse_error": 0,
            "missing_fields": 0,
            "size_error": 0,
            "encoding_error": 0
        }
    
    def parse_experiment_xml(self, xml_path: Path) -> Optional[Dict[str, Any]]:
        """
        Parsear experiment.xml
        
        Campos extraídos:
        - ACCESSION (ID del experimento)
        - TITLE (título del estudio)
        - STUDY_REF → BioProject
        - LIBRARY_STRATEGY (estrategia: RNA-Seq, WGS, etc.)
        - LIBRARY_SOURCE (fuente: TRANSCRIPTOMIC, GENOMIC, etc.)
        - LIBRARY_SELECTION (selección: cDNA, RANDOM, etc.)
        - LIBRARY_LAYOUT (PAIRED o SINGLE)
        """
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            experiments = root.findall(".//EXPERIMENT")
            results = []
            
            for exp in experiments:
                accession = exp.get("accession", "")
                
                # Extraer campos
                title = exp.findtext(".//TITLE", "").strip()
                
                bioproject = exp.findtext(".//EXTERNAL_ID[@label='BioProject']", "") or \
                            exp.findtext(".//EXTERNAL_ID[@namespace='BioProject']", "")
                
                library_strategy = exp.findtext(".//LIBRARY_STRATEGY", "").strip()
                library_source = exp.findtext(".//LIBRARY_SOURCE", "").strip()
                library_selection = exp.findtext(".//LIBRARY_SELECTION", "").strip()
                
                # Library layout
                layout = "SINGLE"
                if exp.find(".//LIBRARY_LAYOUT/PAIRED") is not None:
                    layout = "PAIRED"
                
                # Library name (a veces indica tejido)
                library_name = exp.findtext(".//LIBRARY_NAME", "").strip()
                
                # Library construction protocol (puede contener info importante)
                protocol = exp.findtext(".//LIBRARY_CONSTRUCTION_PROTOCOL", "").strip()
                
                result = {
                    "accession": accession,
                    "title": title,
                    "bioproject": bioproject,
                    "library_strategy": library_strategy,
                    "library_source": library_source,
                    "library_selection": library_selection,
                    "library_layout": layout,
                    "library_name": library_name,
                    "protocol_snippet": protocol[:200] if protocol else ""  # Primeros 200 chars
                }
                
                # Validación
                if not all([accession, title, library_strategy]):
                    self.errors["missing_fields"] += 1
                    continue
                
                results.append(result)
            
            return {"experiments": results, "count": len(results)}
        
        except ET.ParseError as e:
            self.errors["parse_error"] += 1
            self.logger.debug(f"Parse error: {xml_path.name} - {str(e)[:100]}")
            return None
        except Exception as e:
            self.errors["parse_error"] += 1
            self.logger.debug(f"Unexpected error parsing {xml_path.name}: {str(e)[:100]}")
            return None
    
    def parse_sample_xml(self, xml_path: Path) -> Optional[Dict[str, Any]]:
        """
        Parsear sample.xml
        
        Campos extraídos:
        - ACCESSION (ID de la muestra)
        - TAXON_ID
        - SCIENTIFIC_NAME (organismo)
        - SAMPLE_ATTRIBUTES (key-value pairs con metadata)
        """
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            samples = root.findall(".//SAMPLE")
            results = []
            
            for sample in samples:
                accession = sample.get("accession", "")
                
                # Información del organismo
                taxon_id = sample.findtext(".//TAXON_ID", "")
                scientific_name = sample.findtext(".//SCIENTIFIC_NAME", "").strip()
                
                # Atributos de la muestra
                attributes = {}
                for attr in sample.findall(".//SAMPLE_ATTRIBUTE"):
                    tag = attr.findtext("TAG", "").lower().strip()
                    value = attr.findtext("VALUE", "").strip()
                    if tag and value:
                        if tag not in attributes:
                            attributes[tag] = []
                        attributes[tag].append(value)
                
                result = {
                    "accession": accession,
                    "taxon_id": taxon_id,
                    "scientific_name": scientific_name,
                    "attributes": attributes,
                    "attribute_keys": list(attributes.keys())
                }
                
                # Validación
                if not scientific_name:
                    self.errors["missing_fields"] += 1
                    continue
                
                results.append(result)
            
            return {"samples": results, "count": len(results)}
        
        except ET.ParseError as e:
            self.errors["parse_error"] += 1
            self.logger.debug(f"Parse error: {xml_path.name} - {str(e)[:100]}")
            return None
        except Exception as e:
            self.errors["parse_error"] += 1
            self.logger.debug(f"Unexpected error parsing {xml_path.name}: {str(e)[:100]}")
            return None
    
    def parse_run_xml(self, xml_path: Path) -> Optional[Dict[str, Any]]:
        """
        Parsear run.xml
        
        Campos extraídos:
        - RUN accession
        - Fecha de ejecución
        - Center name
        """
        try:
            tree = ET.parse(xml_path)
            root = tree.getroot()
            
            runs = root.findall(".//RUN")
            results = []
            
            for run in runs:
                accession = run.get("accession", "")
                run_date = run.get("run_date", "")
                center_name = run.get("run_center", "") or run.get("center_name", "")
                
                result = {
                    "accession": accession,
                    "run_date": run_date,
                    "center_name": center_name
                }
                
                results.append(result)
            
            return {"runs": results, "count": len(results)}
        
        except ET.ParseError as e:
            self.errors["parse_error"] += 1
            self.logger.debug(f"Parse error: {xml_path.name} - {str(e)[:100]}")
            return None
        except Exception as e:
            self.errors["parse_error"] += 1
            return None


# ============================================
# INDEXING
# ============================================

class IndexBuilder:
    """Construir índices mientras se parsean archivos"""
    
    def __init__(self):
        self.organisms = defaultdict(lambda: {"count": 0, "studies": set()})
        self.strategies = defaultdict(lambda: {"count": 0, "organisms": defaultdict(int)})
        self.sources = defaultdict(int)
        self.selections = defaultdict(int)
        self.attributes_seen = defaultdict(int)
        
        # Estadísticas
        self.total_experiments = 0
        self.total_samples = 0
        self.total_runs = 0
    
    def add_experiment_data(self, exp_data: Dict[str, Any], study_id: str):
        """Agregar datos de experimento a índices"""
        strategy = exp_data.get("library_strategy", "UNKNOWN").strip()
        source = exp_data.get("library_source", "UNKNOWN").strip()
        selection = exp_data.get("library_selection", "UNKNOWN").strip()
        
        self.strategies[strategy]["count"] += 1
        self.sources[source] += 1
        self.selections[selection] += 1
        
        self.total_experiments += 1
    
    def add_sample_data(self, sample_data: Dict[str, Any], study_id: str):
        """Agregar datos de muestra a índices"""
        organism = sample_data.get("scientific_name", "UNKNOWN").strip()
        
        self.organisms[organism]["count"] += 1
        self.organisms[organism]["studies"].add(study_id)
        
        # Registrar atributos
        for attr_key in sample_data.get("attribute_keys", []):
            self.attributes_seen[attr_key] += 1
        
        self.total_samples += 1
    
    def get_indices(self) -> Dict[str, Dict]:
        """Obtener índices finales"""
        return {
            "organisms": dict(self.organisms),
            "strategies": dict(self.strategies),
            "sources": dict(self.sources),
            "selections": dict(self.selections),
            "attributes": dict(self.attributes_seen),
            "stats": {
                "total_experiments": self.total_experiments,
                "total_samples": self.total_samples,
                "total_runs": self.total_runs,
                "unique_organisms": len(self.organisms),
                "unique_strategies": len(self.strategies),
                "unique_sources": len(self.sources),
                "unique_selections": len(self.selections),
                "unique_attributes": len(self.attributes_seen)
            }
        }


# ============================================
# MAIN EXECUTION
# ============================================

def main():
    """Ejecutar Stage 1 - XML Parsing"""
    
    # Setup
    logger = setup_logging("stage1_parse_xml")
    parser = XMLParser(logger)
    index_builder = IndexBuilder()
    checkpoint_mgr = CheckpointManager("stage1", logger)
    
    start_time = time.time()
    
    print_section_header(logger, "STAGE 1: XML PARSING")
    
    # ============ CHECKPOINT: Recuperar estado anterior si existe ============
    logger.info("🔄 Verificando checkpoint anterior...")
    checkpoint = checkpoint_mgr.get_last_checkpoint()
    
    processed_from_checkpoint = 0
    if checkpoint:
        processed_from_checkpoint = checkpoint["data"].get("processed", 0)
        logger.info(f"✅ Recuperando desde archivo {processed_from_checkpoint:,}")
    
    # ============ DISCOVERY: Encontrar todos los BioProject IDs ============
    print_section_header(logger, "DISCOVERY: Encontrando archivos XML")
    
    logger.info(f"📁 Buscando XMLs en: {NCBI_SRA_PATH}")
    bioproject_dirs = sorted([d for d in NCBI_SRA_PATH.iterdir() if d.is_dir()])
    
    logger.info(f"📊 Total de BioProjects encontrados: {len(bioproject_dirs):,}")
    
    # Limitar para esta etapa
    max_files = MAX_FILES_PHASE1_STAGE1
    bioproject_dirs = bioproject_dirs[:max_files]
    
    logger.info(f"🎯 Procesando primeros: {len(bioproject_dirs):,} BioProjects")
    
    # ============ PARSING: Procesar archivos XML ============
    print_section_header(logger, "PARSING: Extrayendo metadatos")
    
    progress = ProgressTracker(len(bioproject_dirs), logger, "XML Parsing")
    monitor = SystemMonitor(logger)
    
    for i, bioproject_dir in enumerate(bioproject_dirs):
        
        # Debug print cada N archivos
        if i % DEBUG_PRINT_INTERVAL == 0:
            measure = monitor.measure()
            logger.debug(
                f"[{i:,}/{len(bioproject_dirs):,}] "
                f"RAM: {measure['ram_mb']/1024:.1f}GB | "
                f"CPU: {measure['cpu_percent']:.1f}%"
            )
        
        bioproject_id = bioproject_dir.name
        
        try:
            # Parsear los 3 tipos de XML
            exp_file = bioproject_dir / f"{bioproject_id}.experiment.xml"
            sample_file = bioproject_dir / f"{bioproject_id}.sample.xml"
            run_file = bioproject_dir / f"{bioproject_id}.run.xml"
            
            if exp_file.exists():
                exp_data = parser.parse_experiment_xml(exp_file)
                if exp_data:
                    for exp in exp_data.get("experiments", []):
                        index_builder.add_experiment_data(exp, bioproject_id)
            
            if sample_file.exists():
                sample_data = parser.parse_sample_xml(sample_file)
                if sample_data:
                    for sample in sample_data.get("samples", []):
                        index_builder.add_sample_data(sample, bioproject_id)
            
            if run_file.exists():
                run_data = parser.parse_run_xml(run_file)
                if run_data:
                    index_builder.total_runs += run_data.get("count", 0)
            
            progress.update(1)
        
        except Exception as e:
            logger.error(f"❌ Error procesando {bioproject_id}: {str(e)[:100]}")
            progress.update(1, error=True)
        
        # Checkpoint cada N archivos
        if (i + 1) % CHECKPOINT_INTERVAL == 0:
            logger.info(f"💾 Guardando checkpoint en archivo {i+1:,}...")
            checkpoint_mgr.save_checkpoint(
                {"processed": i + 1, "indices": index_builder.get_indices()},
                i + 1
            )
    
    progress.print_progress(force=True)
    
    # ============ SAVING: Guardar índices ============
    print_section_header(logger, "SAVING: Escribiendo índices")
    
    logger.info("💾 Guardando índices en JSON...")
    final_indices = index_builder.get_indices()
    
    # Convertir sets a listas para JSON serialization
    for org_name, org_data in final_indices["organisms"].items():
        org_data["studies"] = list(org_data["studies"])
    
    output_file = OUTPUT_DIR / "stage1_indices.json"
    with open(output_file, "w") as f:
        json.dump(final_indices, f, indent=2)
    
    logger.info(f"✅ Índices guardados: {output_file}")
    
    # ============ FINAL STATS ============
    elapsed = time.time() - start_time
    
    stats = {
        "Total Experiments": f"{final_indices['stats']['total_experiments']:,}",
        "Total Samples": f"{final_indices['stats']['total_samples']:,}",
        "Unique Organisms": f"{final_indices['stats']['unique_organisms']:,}",
        "Unique Strategies": f"{final_indices['stats']['unique_strategies']:,}",
        "Parse Errors": parser.errors["parse_error"],
        "Missing Fields": parser.errors["missing_fields"],
        "Files Processed": f"{len(bioproject_dirs):,}",
    }
    
    print_stage_summary(logger, "STAGE 1: XML PARSING", stats, elapsed)
    
    logger.info("✅ Stage 1 completado exitosamente")
    logger.info(f"📊 Archivo de índices: {output_file}")
    logger.info(f"📝 Logs guardados en: {LOGS_DIR}")


if __name__ == "__main__":
    main()
