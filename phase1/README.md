# 🔬 Pyner Phase 1: Knowledge Base Extraction

**Objetivo:** Extraer 1.3M de archivos XML de NCBI SRA y construir un Knowledge Base estructurado.

**Recursos Disponibles:**
- 3x NVIDIA RTX 4000 Ada (~24GB VRAM cada una)
- 251GB RAM
- Multicore CPU

---

## 📋 Estructura

```
phase1/
├── config.py                    ← Configuración centralizada
├── utils.py                     ← Utilidades (logging, GPU, checkpoints)
├── scripts/
│   ├── stage1_parse_xml.py     ← Parsear 1K archivos iniciales (prueba)
│   ├── stage2_parallel_gpu.py   ← Paralelizar 50K archivos
│   └── stage3_parallel_gpu.py   ← Escalar a 500K archivos (máx)
├── checkpoints/                 ← Recuperación ante fallos
├── logs/                        ← Logs detallados de ejecución
├── output/                      ← Índices JSON generados
└── tests/                       ← Fixtures para testing
```

---

## 🚀 Ejecución Rápida

### Prerequisitos

```bash
# Instalar dependencias
pip install -r requirements-dev.txt

# Verificar configuración
cd phase1
python config.py
```

### Stage 1: Prueba Inicial (1K archivos)

```bash
cd phase1

# Ejecutar parser básico
python scripts/stage1_parse_xml.py

# Output esperado:
# - phase1/output/stage1_indices.json
# - phase1/logs/stage1_parse_xml_YYYYMMDD_HHMMSS.log
# 
# Tiempo esperado: ~2-5 minutos
# Resultado: Fase 1 stage 1 completada ✅
```

**Monitor mientras se ejecuta:**
```bash
# En otra terminal
tail -f phase1/logs/stage1_parse_xml_*.log
```

### Stage 2: Paralelo GPU - 50K Archivos

```bash
cd phase1

# Ejecutar con 8 workers y 3 GPUs
python scripts/stage2_parallel_gpu.py --stage 2

# Output:
# - phase1/output/stage2_knowledge_base.json
# - phase1/checkpoints/stage2/ (recuperación)
#
# Tiempo esperado: ~5-10 minutos
# Throughput: ~100-150 archivos/seg
```

### Stage 3: Escalar a 500K Archivos

```bash
cd phase1

# Procesar 500K archivos (máximo para esta fase)
python scripts/stage2_parallel_gpu.py --stage 3

# Output:
# - phase1/output/stage3_knowledge_base.json
# - phase1/logs/stage3_parallel_*.log
#
# Tiempo esperado: ~50-100 minutos
# Sistema completo con checkpoint recovery
```

---

## 🎯 Fases y Checkpoints

### Estructura de Fases

```
PHASE 1: KNOWLEDGE BASE EXTRACTION
│
├── STAGE 1: Parse XML (1K test)
│   ├── ✅ Discover BioProjects
│   ├── ✅ Parse experiment.xml
│   ├── ✅ Parse sample.xml
│   ├── ✅ Parse run.xml
│   ├── ✅ Build initial indices
│   └── ✅ Save output JSON
│
├── STAGE 2: Parallel GPU (50K production)
│   ├── ✅ Discovery (50K BioProjects)
│   ├── ✅ Init 8 worker processes
│   ├── ✅ Assign 3 GPUs (round-robin)
│   ├── ✅ Batch distribution
│   ├── ✅ Parallel processing
│   ├── ✅ Aggregation
│   ├── ✅ Generate complete KB
│   └── ✅ Checkpoint every 10K files
│
└── STAGE 3: Full Scale (500K production)
     ├── ✅ Same as Stage 2
     ├── ✅ Longer runtime (~1-2 hours)
     ├── ✅ Recovery from checkpoints
     └── ✅ Final KB generation
```

### Checkpoint System

Si hay fallo, recuperación automática:

```bash
# Si Stage 2 falla en archivo 25000:
❌ Timeout/Error durante procesamiento

# Ejecutar de nuevo:
python scripts/stage2_parallel_gpu.py --stage 2

# Sistema automáticamente:
✅ Detecta checkpoint en archive 25000
✅ Resume desde ahí
✅ Continúa procesamiento
✅ Sin perder trabajo anterior
```

---

## 📊 Debug & Monitoring

Cada script imprime:
1. **Print de sección** - Qué está haciendo
2. **Print de progreso** - Cada 100 archivos
3. **Print de recursos** - RAM, CPU, GPU cada 5 segundos
4. **Print de errores** - Con contexto completo

### Ejemplo de Output:

```
======================================================================
📍 STAGE 1: XML PARSING
======================================================================
🔄 Verificando checkpoint anterior...
📁 Buscando XMLs en: /home/lahumada/disco1/NCBI_Metadata/SRA
📊 Total de BioProjects encontrados: 445,489

🎯 Procesando primeros: 1,000 BioProjects

======================================================================
📍 PARSING: Extrayendo metadatos
======================================================================
[DEBUG] @elapsed=0.1s | XML Parsing | [████░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 10.0% 
        | Generated: 1,234 organisms, 45 strategies

[DEBUG] @elapsed=5.2s | XML Parsing | RAM: 4.5GB | CPU: 45.2%
        | Processed: 100/1,000 files

[DEBUG] @elapsed=10.3s | XML Parsing | [████████░░░░░░░░░░░░░░░░░░░░░░░░░░░░] 20.0%
        | Rate: 50.3 files/sec | ETA: 3.2 min

... (cada 100 archivos)

======================================================================
✅ STAGE 1 COMPLETADO
======================================================================
  Total Experiments:...................... 5,432
  Unique Organisms:....................... 1,234
  Unique Strategies:...................... 45
  Parse Errors:........................... 23
  Tiempo total:........................... 125.45 sec (2.09 min)
======================================================================
```

---

## 🔍 Verificar Resultados

### Output JSON Structure:

```bash
# Ver tamaño de output
ls -lh phase1/output/

# Ejemplo:
# stage1_indices.json      ~2MB
# stage2_knowledge_base.json ~50MB
# stage3_knowledge_base.json ~500MB (final KB)

# Verificar contenido
head -100 phase1/output/stage1_indices.json | python -m json.tool

# Contar organismos únicos
python -c "
import json
data = json.load(open('phase1/output/stage1_indices.json'))
print(f'Organismos: {data[\"stats\"][\"unique_organisms\"]}')
print(f'Estrategias: {data[\"stats\"][\"unique_strategies\"]}')
"
```

---

## ⚙️ Configuración Avanzada

### Ajustar parámetros

Editar `phase1/config.py`:

```python
# Línea ~30: Cambiar número máximo de archivos
MAX_FILES_PHASE1_STAGE1 = 5000   # Aumentar a 5K para prueba más grande

# Línea ~35: Paralelización
NUM_WORKERS = 16  # Más workers = más rápido (si hay CPU)
BATCH_SIZE = 50   # Batches más pequeños = mejor distribución

# Línea ~45: GPU
GPU_IDS = [0, 1, 2]  # Usar todas 3 GPUs
GPU_MEMORY_FRACTION = 0.9  # Usar 90% de VRAM

# Línea ~57: Checkpoint
CHECKPOINT_INTERVAL = 5000  # Checkpoint cada 5K files

# Línea ~70: Debug
DEBUG_PRINT_INTERVAL = 50  # Print cada 50 files (más verbose)
```

Luego ejecutar:

```bash
python scripts/stage2_parallel_gpu.py --stage 2
```

---

## 🐛 Troubleshooting

### Problema: "GPU out of memory"

```
torch.cuda.OutOfMemoryError: CUDA out of memory

Solución:
1. Reducir GPU_MEMORY_FRACTION en config.py
2. Reducir BATCH_SIZE
3. Reducir NUM_WORKERS
```

### Problema: "Stuck on checkpoint"

```
⏱️ Worker timeout esperando trabajo

Solución:
1. Aumentar timeout: config.py línea ~50
2. Reducir NUM_WORKERS
3. Killer manualmente: pkill -f stage2_parallel_gpu.py
```

### Problema: "Output JSON demasiado grande"

```
El archivo JSON es > 1GB

Solución:
1. Usar stage2 en lugar de stage3 (50K vs 500K)
2. Cambiar output format a parquet (más compactado)
3. Dividir en múltiples JSONs
```

---

## 📈 Performance Expectations

| Stage | Files | Workers | GPUs | Time | Throughput |
|-------|-------|---------|------|------|------------|
| 1 | 1K | 1 | 0 | 2-5 min | ~200 files/sec |
| 2 | 50K | 8 | 3 | 5-10 min | ~100-150 files/sec |
| 3 | 500K | 8 | 3 |50-100 min | ~80-150 files/sec |

Con 3x RTX 4000 Ada y 8 workers, esperamos **~100-150 files/sec**.

---

## 🎓 Output Knowledge Base (KB)

Después de cualquier stage, habrá JSON con estructura:

```json
{
  "stage": 1,
  "timestamp": "2026-02-06T14:30:00",
  "files_processed": 1000,
  "statistics": {
    "total_experiments": 5432,
    "unique_organisms": 1234,
    "unique_strategies": 45
  },
  "organisms": {
    "arabidopsis_thaliana": {"count": 234, "studies": 145},
    "solanum_lycopersicum": {"count": 89, "studies": 54}
  },
  "strategies": {
    "RNA-Seq": 2345,
    "WGS": 1234
  }
}
```

---

## ✅ Checklist de Validación

- [ ] `python config.py` ejecuta sin errores
- [ ] Logs se generan en `phase1/logs/`
- [ ] Checkpoints se crean en `phase1/checkpoints/` cada 10K archivos
- [ ] Output JSON generado en `phase1/output/`
- [ ] JSON es válido (no corrupto)
- [ ] GPU memory se libera después de terminar
- [ ] Recuperación automática funciona si se interrumpe

---

## 📞 Debug

Si algo va mal:

```bash
# 1. Ver último log
tail -50 phase1/logs/*.log

# 2. Ver checkpoints disponibles
ls -la phase1/checkpoints/stage2/

# 3. Ver output generado
ls -lh phase1/output/

# 4. Validar JSON
python -m json.tool phase1/output/stage1_indices.json > /dev/null

# 5. Ver estadísticas en tiempo real
watch -n 1 'ls -la phase1/logs/stage2_parallel_*.log | tail -1'
```

---

**¡Listo para ejecutar! 🚀** 

Comienza con Stage 1 para validar setup:
```bash
cd phase1
python scripts/stage1_parse_xml.py
```
