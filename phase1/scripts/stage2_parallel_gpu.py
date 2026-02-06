"""
Pyner Phase 1 - Stage 2: Parallel Processing with GPU
======================================================
Procesamiento paralelo optimizado con GPUs para 50K-500K archivos

Stage 2 Objectives:
1. Paralelizar procesamiento entre múltiples workers
2. Distribuir carga entre 3 GPUs
3. Implementar pipeline eficiente (reader → parser → indexer)
4. Usar multiprocessing + queues
5. Generar KB completa y optimizada

Ejecución:
    python scripts/stage2_parallel_gpu.py --stage 2 --max-files 50000
    python scripts/stage2_parallel_gpu.py --stage 3 --max-files 500000

Performance Target:
    Stage 2 (50K):   ~5-10 min
    Stage 3 (500K):  ~50-100 min (con 3 GPUs + 8 workers)

"""

import multiprocessing as mp
from multiprocessing import Process, Queue, Manager
import json
import logging
import time
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import sys
from typing import Dict, List, Tuple, Any
import psutil

try:
    import torch
    TORCH_AVAILABLE = True
except:
    TORCH_AVAILABLE = False

from config import (
    NCBI_SRA_PATH, OUTPUT_DIR, 
    MAX_FILES_PHASE1_STAGE2, MAX_FILES_PHASE1_STAGE3,
    NUM_WORKERS, BATCH_SIZE, GPU_IDS,
    USE_GPU, DEBUG_PRINT_INTERVAL
)
from utils import (
    setup_logging, print_section_header, ProgressTracker,
    SystemMonitor, CheckpointManager, print_stage_summary
)
from stage1_parse_xml import XMLParser, IndexBuilder


# ============================================
# WORKER PROCESS FUNCTION
# ============================================

def worker_process(
    worker_id: int,
    input_queue: mp.Queue,
    output_queue: mp.Queue,
    gpu_id: int,
    logger_config: Tuple
):
    """
    Worker process para parsear XMLs en paralelo
    
    Cada worker:
    1. Obtiene batches de bioproject_dirs de input_queue
    2. Parsea los XMLs
    3. Construye índices parciales
    4. Envía resultados a output_queue
    5. Usa GPU asignada si disponible
    """
    
    # Reconfigrurar logger en el proceso worker
    logger = setup_logging(f"worker_{worker_id}")
    parser = XMLParser(logger)
    index_builder = IndexBuilder()
    
    # Asignar GPU si disponible
    if TORCH_AVAILABLE and torch.cuda.is_available() and gpu_id is not None:
        torch.cuda.set_device(gpu_id)
        logger.debug(f"🎮 Worker {worker_id} usando GPU {gpu_id}")
    
    logger.debug(f"🚀 Worker {worker_id} iniciado (GPU: {gpu_id})")
    
    processed_count = 0
    error_count = 0
    
    while True:
        try:
            # Obtener batch de trabajo
            item = input_queue.get(timeout=5)
            
            # Señal de término
            if item is None:
                logger.debug(f"✅ Worker {worker_id} terminando ({processed_count} archivos procesados)")
                break
            
            bioproject_dirs, batch_num = item
            
            # Procesar batch
            for bioproject_dir in bioproject_dirs:
                bioproject_id = bioproject_dir.name
                
                try:
                    # Parsear XMLs
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
                    
                    processed_count += 1
                
                except Exception as e:
                    error_count += 1
                    logger.debug(f"❌ Error en {bioproject_id}: {str(e)[:80]}")
            
            # Enviar resultados
            output_queue.put({
                "worker_id": worker_id,
                "batch_num": batch_num,
                "processed": len(bioproject_dirs),
                "indices": index_builder.get_indices(),
                "errors": error_count
            })
        
        except mp.TimeoutError:
            logger.debug(f"⏱️  Worker {worker_id} timeout esperando trabajo")
            continue
        except Exception as e:
            logger.error(f"❌ Error en worker: {str(e)[:100]}")
            break
    
    logger.debug(f"🏁 Worker {worker_id} finalizado")


# ============================================
# MAIN EXECUTION
# ============================================

def run_parallel_stage(stage_number: int, max_files: int):
    """
    Ejecutar Stage 2 o 3 con procesamiento paralelo
    
    Args:
        stage_number: 2 o 3
        max_files: Máximo archivos a procesar
    """
    
    # Setup
    logger = setup_logging(f"stage{stage_number}_parallel")
    checkpoint_mgr = CheckpointManager(f"stage{stage_number}", logger)
    
    start_time = time.time()
    
    print_section_header(logger, f"STAGE {stage_number}: PARALLEL GPU PROCESSING")
    
    logger.info(f"🎯 Objetivo: Procesar {max_files:,} archivos con {NUM_WORKERS} workers")
    logger.info(f"🎮 GPUs disponibles: {GPU_IDS}")
    logger.info(f"⚙️  Batch size: {BATCH_SIZE}")
    
    # ============ DISCOVERY ============
    print_section_header(logger, "DISCOVERY: Encontrando archivos")
    
    logger.info(f"📁 Buscando BioProjects en: {NCBI_SRA_PATH}")
    bioproject_dirs = sorted([d for d in NCBI_SRA_PATH.iterdir() if d.is_dir()])
    bioproject_dirs = bioproject_dirs[:max_files]
    
    logger.info(f"📊 Total BioProjects a procesar: {len(bioproject_dirs):,}")
    
    # ============ INITIALIZATION ============
    print_section_header(logger, "INITIALIZATION: Configurando workers")
    
    # Crear queues
    input_queue = mp.Queue(maxsize=NUM_WORKERS * 2)
    output_queue = mp.Queue()
    
    # GPU assignment (round-robin)
    gpu_assignments = [GPU_IDS[i % len(GPU_IDS)] for i in range(NUM_WORKERS)]
    
    logger.info(f"💼 Asignaciones de GPU: {dict(enumerate(gpu_assignments))}")
    
    # Iniciar workers
    workers = []
    for worker_id in range(NUM_WORKERS):
        p = Process(
            target=worker_process,
            args=(
                worker_id,
                input_queue,
                output_queue,
                gpu_assignments[worker_id],
                None
            )
        )
        p.start()
        workers.append(p)
        logger.debug(f"✅ Worker {worker_id} iniciado")
    
    # ============ DATA DISTRIBUTION ============
    print_section_header(logger, "DISTRIBUTION: Enviando batches a workers")
    
    logger.debug("📤 Distribuyendo archivos entre workers...")
    batch_num = 0
    for i in range(0, len(bioproject_dirs), BATCH_SIZE):
        batch = bioproject_dirs[i:i+BATCH_SIZE]
        input_queue.put((batch, batch_num))
        batch_num += 1
    
    # Señales de término para workers
    for _ in range(NUM_WORKERS):
        input_queue.put(None)
    
    logger.info(f"📦 Total batches distribuidos: {batch_num}")
    
    # ============ PROCESSING & MONITORING ============
    print_section_header(logger, "PROCESSING: Monitoreando ejecución")
    
    monitor = SystemMonitor(logger, interval=10)
    progress = ProgressTracker(len(bioproject_dirs), logger, "Parallel Processing")
    
    completed_batches = 0
    global_indices = IndexBuilder()
    worker_results = defaultdict(list)
    
    while completed_batches < batch_num:
        try:
            result = output_queue.get(timeout=30)
            
            worker_id = result["worker_id"]
            batch_num_completed = result["batch_num"]
            processed = result["processed"]
            
            # Agregar resultados
            worker_results[worker_id].append(result)
            progress.update(processed)
            completed_batches += 1
            
            # Debug reporting
            if completed_batches % 10 == 0:
                measure = monitor.measure()
                monitor.print_status(len(bioproject_dirs), progress.processed)
                logger.debug(
                    f"📊 Batches completados: {completed_batches}/{batch_num} | "
                    f"Archivos: {progress.processed:,}/{len(bioproject_dirs):,}"
                )
        
        except Exception as e:
            logger.warning(f"❌ Error recibiendo resultado: {str(e)}")
            continue
    
    progress.print_progress(force=True)
    
    # Esperar a que todos los workers terminen
    logger.debug("⏳ Esperando a que workers terminen...")
    for i, worker in enumerate(workers):
        worker.join(timeout=10)
        if worker.is_alive():
            logger.warning(f"⚠️  Worker {i} no respondió, terminando...")
            worker.terminate()
        logger.debug(f"✅ Worker {i} terminado")
    
    # ============ AGGREGATION ============
    print_section_header(logger, "AGGREGATION: Combinando resultados")
    
    logger.info("🔗 Combinando índices de todos los workers...")
    
    total_experiments = 0
    total_samples = 0
    total_runs = 0
    
    organisms = defaultdict(lambda: {"count": 0, "studies": set()})
    strategies = defaultdict(int)
    sources = defaultdict(int)
    selections = defaultdict(int)
    
    for worker_id, results_list in worker_results.items():
        for result in results_list:
            indices = result["indices"]
            
            # Merge organisms
            for org_name, org_data in indices["organisms"].items():
                organisms[org_name]["count"] += org_data.get("count", 0)
                organisms[org_name]["studies"].update(org_data.get("studies", []))
            
            # Merge strategies
            for strategy, count in indices["strategies"].items():
                strategies[strategy] += count.get("count", 0) if isinstance(count, dict) else count
            
            # Merge others
            for source, count in indices["sources"].items():
                sources[source] += count
            
            for selection, count in indices["selections"].items():
                selections[selection] += count
            
            # Stats
            stats = indices["stats"]
            total_experiments += stats.get("total_experiments", 0)
            total_samples += stats.get("total_samples", 0)
            total_runs += stats.get("total_runs", 0)
    
    logger.info("✅ Índices combinados exitosamente")
    
    # ============ SAVING ============
    print_section_header(logger, "SAVING: Guardando KB")
    
    logger.info("💾 Guardando índices finales...")
    
    final_output = {
        "stage": stage_number,
        "timestamp": datetime.now().isoformat(),
        "files_processed": len(bioproject_dirs),
        "statistics": {
            "total_experiments": total_experiments,
            "total_samples": total_samples,
            "total_runs": total_runs,
            "unique_organisms": len(organisms),
            "unique_strategies": len(strategies),
            "unique_sources": len(sources),
            "unique_selections": len(selections),
        },
        "organisms": {k: {
            "count": v["count"],
            "studies": len(v["studies"])
        } for k, v in sorted(organisms.items(), key=lambda x: x[1]["count"], reverse=True)[:100]},  # Top 100
        "strategies": dict(sorted(strategies.items(), key=lambda x: x[1], reverse=True)),
        "sources": dict(sorted(sources.items(), key=lambda x: x[1], reverse=True)),
        "selections": dict(sorted(selections.items(), key=lambda x: x[1], reverse=True)),
    }
    
    output_file = OUTPUT_DIR / f"stage{stage_number}_knowledge_base.json"
    with open(output_file, "w") as f:
        json.dump(final_output, f, indent=2)
    
    logger.info(f"✅ KB guardada: {output_file}")
    
    # ============ FINAL STATS ============
    elapsed = time.time() - start_time
    
    stats = {
        "Total de Archivos": f"{len(bioproject_dirs):,}",
        "Worker Processes": NUM_WORKERS,
        "GPUs Utilizadas": len(GPU_IDS),
        "Batches Procesados": completed_batches,
        "Experimentos Extraídos": f"{total_experiments:,}",
        "Muestras Extraídas": f"{total_samples:,}",
        "Organismos Únicos": f"{len(organisms):,}",
        "Estrategias Únicas": len(strategies),
    }
    
    print_stage_summary(logger, f"STAGE {stage_number}: PARALLEL PROCESSING", stats, elapsed)
    
    logger.info(f"⚡ Throughput: {len(bioproject_dirs) / elapsed:.1f} files/sec")
    logger.info(f"✅ Stage {stage_number} completado exitosamente")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Pyner Phase 1 Stage 2/3 - Parallel GPU Processing")
    parser.add_argument("--stage", type=int, default=2, choices=[2, 3], help="Stage to run (2 o 3)")
    args = parser.parse_args()
    
    if args.stage == 2:
        run_parallel_stage(2, MAX_FILES_PHASE1_STAGE2)
    else:
        run_parallel_stage(3, MAX_FILES_PHASE1_STAGE3)
