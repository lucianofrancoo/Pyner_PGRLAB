"""
Pyner Phase 1 - Utilities Module
================================
Utilidades para logging, GPU management, checkpoints y monitoreo

"""

import logging
import json
import psutil
import time
import pickle
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional
import sys

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False

from config import (
    LOGS_DIR, CHECKPOINTS_DIR, LOG_LEVEL, LOG_FORMAT,
    GPU_IDS, USE_GPU, ENABLE_PROFILING, PROFILING_INTERVAL,
    MEMORY_MONITORING, DEBUG_PRINT_INTERVAL
)


# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(stage_name: str) -> logging.Logger:
    """
    Configurar logging para una fase específica
    
    Args:
        stage_name: Nombre de la fase (ej: "parse_xml", "build_kb")
    
    Returns:
        Logger configurado
    """
    logger = logging.getLogger(f"Pyner.Phase1.{stage_name}")
    logger.setLevel(getattr(logging, LOG_LEVEL))
    
    # Handler para archivo
    log_file = LOGS_DIR / f"{stage_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    fh = logging.FileHandler(log_file)
    fh.setLevel(getattr(logging, LOG_LEVEL))
    
    # Handler para consola
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(getattr(logging, LOG_LEVEL))
    
    # Formatter
    formatter = logging.Formatter(LOG_FORMAT)
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger


# ============================================
# GPU MANAGEMENT
# ============================================

def check_gpu_availability() -> Dict[int, Dict[str, Any]]:
    """
    Verificar disponibilidad y características de GPUs
    
    Returns:
        Diccionario con info de cada GPU
    """
    gpu_info = {}
    
    if not TORCH_AVAILABLE:
        return {"error": "PyTorch no instalado"}
    
    if not torch.cuda.is_available():
        return {"available": False, "count": 0}
    
    for gpu_id in GPU_IDS:
        try:
            props = torch.cuda.get_device_properties(gpu_id)
            memory_total = props.total_memory / 1e9  # GB
            memory_allocated = torch.cuda.memory_allocated(gpu_id) / 1e9
            memory_cached = torch.cuda.memory_cached(gpu_id) / 1e9
            memory_free = memory_total - memory_allocated - memory_cached
            
            gpu_info[gpu_id] = {
                "name": props.name,
                "capability": f"{props.major}.{props.minor}",
                "total_memory_gb": round(memory_total, 2),
                "allocated_gb": round(memory_allocated, 2),
                "cached_gb": round(memory_cached, 2),
                "free_gb": round(memory_free, 2),
                "utilization_percent": round((memory_allocated / memory_total) * 100, 1)
            }
        except Exception as e:
            gpu_info[gpu_id] = {"error": str(e)}
    
    return gpu_info


def print_gpu_status(logger: logging.Logger):
    """Imprimir estado de GPUs en log"""
    gpu_info = check_gpu_availability()
    
    if "error" in gpu_info:
        logger.warning(f"GPU Error: {gpu_info['error']}")
        return
    
    logger.info("GPU STATUS:")
    for gpu_id, info in gpu_info.items():
        if "error" in info:
            logger.warning(f"  GPU {gpu_id}: {info['error']}")
        else:
            logger.info(
                f"  GPU {gpu_id} ({info['name']}): "
                f"{info['free_gb']:.1f}GB free / {info['total_memory_gb']:.1f}GB total "
                f"({info['utilization_percent']:.0f}% util)"
            )


# ============================================
# MEMORY & CPU MONITORING
# ============================================

class SystemMonitor:
    """Monitor de recursos del sistema (CPU, RAM, GPU)"""
    
    def __init__(self, logger: logging.Logger, interval: float = PROFILING_INTERVAL):
        self.logger = logger
        self.interval = interval
        self.start_time = time.time()
        self.process = psutil.Process()
        self.measurements = []
    
    def measure(self) -> Dict[str, Any]:
        """Tomar una medición"""
        elapsed = time.time() - self.start_time
        
        # CPU/RAM
        cpu_percent = self.process.cpu_percent(interval=0.1)
        ram_info = self.process.memory_info()
        ram_mb = ram_info.rss / 1e6
        
        # GPU (si disponible)
        gpu_memory = {}
        if TORCH_AVAILABLE and torch.cuda.is_available():
            for gpu_id in GPU_IDS:
                try:
                    gpu_memory[gpu_id] = round(
                        torch.cuda.memory_allocated(gpu_id) / 1e9, 2
                    )
                except:
                    pass
        
        measurement = {
            "elapsed_sec": round(elapsed, 1),
            "cpu_percent": cpu_percent,
            "ram_mb": round(ram_mb, 1),
            "gpu_memory_gb": gpu_memory
        }
        
        self.measurements.append(measurement)
        return measurement
    
    def print_status(self, file_count: int, processed_count: int):
        """Imprimir estado de recursos"""
        if not self.measurements:
            return
        
        latest = self.measurements[-1]
        
        ram_gb = latest["ram_mb"] / 1024
        self.logger.debug(
            f"📊 Resources | "
            f"CPU: {latest['cpu_percent']:.1f}% | "
            f"RAM: {ram_gb:.1f}GB | "
            f"Processed: {processed_count:,}/{file_count:,} files"
        )
        
        if latest["gpu_memory_gb"]:
            gpu_str = ", ".join(
                f"GPU{gid}: {mem:.1f}GB"
                for gid, mem in latest["gpu_memory_gb"].items()
            )
            self.logger.debug(f"   GPU Memory: {gpu_str}")


# ============================================
# CHECKPOINT SYSTEM
# ============================================

class CheckpointManager:
    """Sistema de checkpoints para recuperación ante fallos"""
    
    def __init__(self, stage_name: str, logger: logging.Logger):
        self.stage_name = stage_name
        self.logger = logger
        self.checkpoint_dir = CHECKPOINTS_DIR / stage_name
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
    
    def get_last_checkpoint(self) -> Optional[Dict[str, Any]]:
        """Obtener último checkpoint si existe"""
        checkpoint_files = sorted(self.checkpoint_dir.glob("checkpoint_*.pkl"))
        
        if not checkpoint_files:
            return None
        
        latest = checkpoint_files[-1]
        try:
            with open(latest, "rb") as f:
                checkpoint = pickle.load(f)
            self.logger.info(f"✅ Checkpoint recuperado: {latest.name}")
            return checkpoint
        except Exception as e:
            self.logger.error(f"❌ Error cargando checkpoint: {e}")
            return None
    
    def save_checkpoint(self, data: Dict[str, Any], file_index: int):
        """Guardar checkpoint"""
        checkpoint_file = self.checkpoint_dir / f"checkpoint_{file_index:06d}.pkl"
        
        try:
            # Agregar metadata
            checkpoint = {
                "stage": self.stage_name,
                "file_index": file_index,
                "timestamp": datetime.now().isoformat(),
                "data": data
            }
            
            with open(checkpoint_file, "wb") as f:
                pickle.dump(checkpoint, f)
            
            self.logger.debug(f"💾 Checkpoint guardado: {checkpoint_file.name}")
        except Exception as e:
            self.logger.error(f"❌ Error guardando checkpoint: {e}")
    
    def save_metadata(self, metadata: Dict[str, Any]):
        """Guardar metadata de la ejecución"""
        metadata_file = LOGS_DIR / f"metadata_{self.stage_name}.json"
        
        try:
            with open(metadata_file, "w") as f:
                json.dump(metadata, f, indent=2)
            self.logger.info(f"📋 Metadata guardada: {metadata_file.name}")
        except Exception as e:
            self.logger.error(f"❌ Error guardando metadata: {e}")


# ============================================
# PROGRESS TRACKING
# ============================================

class ProgressTracker:
    """Rastreador de progreso con ETA"""
    
    def __init__(self, total: int, logger: logging.Logger, name: str = "Progress"):
        self.total = total
        self.logger = logger
        self.name = name
        self.processed = 0
        self.start_time = time.time()
        self.errors = 0
    
    def update(self, count: int = 1, error: bool = False):
        """Actualizar progreso"""
        self.processed += count
        if error:
            self.errors += 1
    
    def print_progress(self, force: bool = False):
        """Imprimir progreso con ETA"""
        if self.processed % DEBUG_PRINT_INTERVAL != 0 and not force:
            return
        
        elapsed = time.time() - self.start_time
        rate = self.processed / elapsed if elapsed > 0 else 0
        eta_sec = (self.total - self.processed) / rate if rate > 0 else 0
        
        percent = (self.processed / self.total) * 100 if self.total > 0 else 0
        
        bar_len = 40
        filled = int(bar_len * self.processed / self.total)
        bar = "█" * filled + "░" * (bar_len - filled)
        
        self.logger.info(
            f"⏳ {self.name} | "
            f"[{bar}] {percent:.1f}% "
            f"({self.processed:,}/{self.total:,}) | "
            f"Rate: {rate:.1f} files/sec | "
            f"ETA: {eta_sec/60:.1f} min | "
            f"Errors: {self.errors}"
        )


# ============================================
# DEBUG HELPERS
# ============================================

def print_section_header(logger: logging.Logger, section: str):
    """Imprimir headers de sección"""
    logger.info("\n" + "="*70)
    logger.info(f"📍 {section}")
    logger.info("="*70)


def print_stage_summary(
    logger: logging.Logger,
    stage_name: str,
    stats: Dict[str, Any],
    elapsed_time: float
):
    """Imprimir resumen de una etapa"""
    logger.info("\n" + "="*70)
    logger.info(f"✅ {stage_name} COMPLETADO")
    logger.info("="*70)
    
    for key, value in stats.items():
        if isinstance(value, float):
            logger.info(f"  {key:.<40} {value:.2f}")
        else:
            logger.info(f"  {key:.<40} {value}")
    
    logger.info(f"  {'Tiempo total':.<40} {elapsed_time:.2f} sec ({elapsed_time/60:.2f} min)")
    logger.info("="*70 + "\n")


if __name__ == "__main__":
    # Test utilities
    logger = setup_logging("test")
    logger.info("Testing utilities...")
    
    print_gpu_status(logger)
    
    monitor = SystemMonitor(logger)
    for _ in range(3):
        m = monitor.measure()
        logger.info(f"Measurement: {m}")
        time.sleep(1)
    
    logger.info("✅ Utilities test completed")
