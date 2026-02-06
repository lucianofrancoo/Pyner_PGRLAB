"""
Pyner Phase 1 - Configuration Module
====================================
Configuración centralizada para toda la Fase 1

Recursos:
- GPUs: 3x NVIDIA RTX 4000 Ada (~24GB VRAM cada una)
- RAM: 251GB
- CPUs: Multicore

"""

import os
from pathlib import Path

# ============================================
# PATHS
# ============================================

PROJECT_ROOT = Path(__file__).parent.parent
PHASE1_ROOT = PROJECT_ROOT / "phase1"
SCRIPTS_DIR = PHASE1_ROOT / "scripts"
CHECKPOINTS_DIR = PHASE1_ROOT / "checkpoints"
LOGS_DIR = PHASE1_ROOT / "logs"
OUTPUT_DIR = PHASE1_ROOT / "output"
TESTS_DIR = PHASE1_ROOT / "tests"

# Crear directorios si no existen
for dir_path in [CHECKPOINTS_DIR, LOGS_DIR, OUTPUT_DIR, TESTS_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

# Data sources
NCBI_SRA_PATH = Path("/home/lahumada/disco1/NCBI_Metadata/SRA")

# ============================================
# PROCESSING PARAMETERS
# ============================================

# Número de archivos a procesar en esta fase
MAX_FILES_PHASE1_STAGE1 = 1000  # Stage 1: Prueba con 1K archivos
MAX_FILES_PHASE1_STAGE2 = 50000  # Stage 2: Escalar a 50K
MAX_FILES_PHASE1_STAGE3 = 500000  # Stage 3: Producción 500K

# Paralelización
NUM_WORKERS = 8  # Procesos paralelos (ajusta según CPU cores)
BATCH_SIZE = 100  # Archivos por batch
QUEUE_SIZE = 50  # Tamaño de cola entre procesos

# GPU
USE_GPU = True
GPU_IDS = [0, 1, 2]  # IDs de 3 GPUs disponibles
GPU_MEMORY_FRACTION = 0.8  # Usar 80% de VRAM de GPU

# ============================================
# CHECKPOINT & RECOVERY
# ============================================

# Sistema de checkpoints
CHECKPOINT_INTERVAL = 10000  # Guardar checkpoint cada 10K archivos procesados
ENABLE_CHECKPOINTS = True
CHECKPOINT_FORMAT = "parquet"  # Formato eficiente para serialización

# ============================================
# LOGGING
# ============================================

LOG_LEVEL = "DEBUG"  # DEBUG, INFO, WARNING, ERROR
LOG_FORMAT = "[%(asctime)s] %(levelname)-8s | %(name)s | %(message)s"
LOG_FILE = LOGS_DIR / f"phase1_execution.log"

# ============================================
# OUTPUT & KB
# ============================================

# Formatos de salida
OUTPUT_FORMAT = "json"  # JSON o parquet
OUTPUT_FILES = {
    "organisms": OUTPUT_DIR / "organisms_index.json",
    "strategies": OUTPUT_DIR / "strategies_index.json",
    "sources": OUTPUT_DIR / "sources_index.json",
    "selections": OUTPUT_DIR / "selections_index.json",
    "tissues": OUTPUT_DIR / "tissues_index.json",
    "treatments": OUTPUT_DIR / "treatments_index.json",
    "global_index": OUTPUT_DIR / "global_index.json",
    "metadata": OUTPUT_DIR / "metadata.json",
}

# ============================================
# DEBUG & MONITORING
# ============================================

ENABLE_PROFILING = True  # Habilitar profiling de CPU/GPU
PROFILING_INTERVAL = 5  # Reporte cada 5 segundos
DEBUG_PRINT_INTERVAL = 100  # Print debug cada 100 archivos procesados
MEMORY_MONITORING = True  # Monitorear RAM y VRAM

# ============================================
# VALIDATION
# ============================================

VALIDATE_XML = True  # Validar XMLs durante parsing
SKIP_BROKEN_XML = True  # Contar pero no procesar XMLs rotos
MAX_XML_SIZE_MB = 50  # Límite máximo de tamaño XML (para detectar anomalías)

# ============================================
# PRINT CONFIGURATION (para verificar setup)
# ============================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("PYNER PHASE 1 - CONFIGURATION")
    print("="*70 + "\n")
    
    print("📁 PATHS:")
    print(f"  Project Root:     {PROJECT_ROOT}")
    print(f"  NCBI SRA Data:    {NCBI_SRA_PATH}")
    print(f"  Checkpoints:      {CHECKPOINTS_DIR}")
    print(f"  Logs:             {LOGS_DIR}")
    print(f"  Output:           {OUTPUT_DIR}\n")
    
    print("⚙️  PROCESSING:")
    print(f"  Max Files (Stage 1): {MAX_FILES_PHASE1_STAGE1:,}")
    print(f"  Max Files (Stage 2): {MAX_FILES_PHASE1_STAGE2:,}")
    print(f"  Max Files (Stage 3): {MAX_FILES_PHASE1_STAGE3:,}")
    print(f"  Workers (Parallel):  {NUM_WORKERS}")
    print(f"  Batch Size:          {BATCH_SIZE}\n")
    
    print("🎮 GPU:")
    print(f"  Use GPU:             {USE_GPU}")
    print(f"  GPU IDs:             {GPU_IDS}")
    print(f"  GPU Memory Frac:     {GPU_MEMORY_FRACTION * 100:.0f}%\n")
    
    print("🔄 CHECKPOINTS:")
    print(f"  Enabled:             {ENABLE_CHECKPOINTS}")
    print(f"  Interval:            {CHECKPOINT_INTERVAL:,} files")
    print(f"  Format:              {CHECKPOINT_FORMAT}\n")
    
    print("📊 DEBUG:")
    print(f"  Log Level:           {LOG_LEVEL}")
    print(f"  Profiling:           {ENABLE_PROFILING}")
    print(f"  Memory Monitor:      {MEMORY_MONITORING}")
    print(f"  Debug Print Interval: {DEBUG_PRINT_INTERVAL} files\n")
    
    print("="*70 + "\n")
