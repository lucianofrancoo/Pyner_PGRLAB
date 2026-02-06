# 🚀 Pyner - Project Structure & Quick Start

Después de reorganizar completamente el proyecto, la estructura es ahora:

```
Pyner_PGRLAB/
│
├── 📖 PLANNING/ (Documentación de planificación)
│   ├── DEVELOPMENT_PLAN.md      ← Plan técnico detallado
│   ├── ROADMAP.md               ← Timeline visual
│   ├── CONTRIBUTING.md          ← Guía para contribuyentes
│   └── KB_SCHEMA.md             ← Esquema del Knowledge Base
│
├── 🔬 PHASE1/ (ACTIVO: Extracción KB)
│   ├── config.py                ← Configuración centralizada
│   ├── utils.py                 ← Logging, GPU, checkpoints
│   ├── README.md                ← 👈 LEER ESTO PRIMERO
│   │
│   ├── scripts/
│   │   ├── stage1_parse_xml.py             ← Prueba (1K archivos)
│   │   └── stage2_parallel_gpu.py          ← Paralelo (50K-500K)
│   │
│   ├── checkpoints/             ← Recovery si falla
│   ├── logs/                    ← Debug detallado
│   ├── output/                  ← Índices JSON generados
│   └── tests/                   ← Test fixtures
│
├── 📁 PHASE2/ (Próxima: Query Optimizer)
├── 📁 PHASE3/ (Última: Optimización)
│
├── 🔬 SCRIPTS ACTUALES (v0.1, v0.2)
│   ├── Pyner_search_v0.1.py
│   ├── Pyner_v0.1.py
│   └── Pyner_v0.2.py
│
├── 📚 DOCUMENTACIÓN
│   ├── README.md                ← Intro general del proyecto
│   ├── INDEX.md                 ← Guía de 5 min
│   │
│   └── docs/
│       └── KB_SCHEMA.md         ← Estructura KB (ahora en planning/)
│
├── requirements.txt             ← Dependencias base
├── requirements-dev.txt         ← Dev/testing
│
└── .github/                     ← Templates GitHub
    └── ISSUE_TEMPLATE/
        ├── bug_report.md
        └── feature_request.md
```

---

## ⚡ Empezar YA (5 minutos)

### 1. Leer Documentación

```bash
# Entendimiento rápido
cat phase1/README.md

# Leer PRIMERO: Descripción general
cat README.md
```

### 2. Instalar Dependencias

```bash
pip install -r requirements.txt
pip install -r requirements-dev.txt

# Si quieres GPU
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
```

### 3. Ejecutar Phase 1 - Etapa 1 (Prueba)

```bash
cd phase1

# Ver configuración
python config.py

# Ejecutar parser (1K archivos = 2-5 minutos)
python scripts/stage1_parse_xml.py

# Monitor en otra terminal
tail -f logs/*.log
```

### 4. Ver Resultados

```bash
# Output generado
ls -lh output/

# Ver índices
python -c "
import json
data = json.load(open('output/stage1_indices.json'))
print('Organismos:', data['stats']['unique_organisms'])
print('Estrategias:', data['stats']['unique_strategies'])
"
```

---

## 🎯 3 Stages de Phase 1

| Stage | Archivos | Tiempo | GPU | Propósito |
|-------|----------|--------|-----|-----------|
| **1** | 1,000 | 2-5 min | No | Validar setup |
| **2** | 50,000 | 5-10 min | Sí | Validar paralelismo |
| **3** | 500,000 | 50-100 min | Sí | KB final |

### Ejecutar cada stage:

```bash
cd phase1

# Stage 1 (validación)
python scripts/stage1_parse_xml.py

# Stage 2 (paralelismo)
python scripts/stage2_parallel_gpu.py --stage 2

# Stage 3 (producción)
python scripts/stage2_parallel_gpu.py --stage 3
```

---

## 🔄 Qué Ver en Logs

Cada ejecución imprime:

```
📊 RECURSOS cada 100 archivos:
  - CPU usage %
  - RAM en GB
  - GPU memory en GB
  - Progreso (%)
  - ETA (minutos)

💾 CHECKPOINTS cada 10K archivos:
  - Estado guardado para recuperación
  
✅ RESUMEN final:
  - Organismos encontrados
  - Estrategias encontradas
  - Errores/warnings
  - Tiempo total
```

---

## 💾 Outputs Esperados

Después de cada stage, en `phase1/output/`:

```
stage1_indices.json          ~2-5 MB
stage2_knowledge_base.json   ~50-100 MB
stage3_knowledge_base.json   ~500 MB (final KB completa)
```

JSON contiene:
- Organismos y sus frecuencias
- Estrategias de secuenciación
- Fuentes de librería
- Selecciones
- Atributos de muestras

---

## 🚨 Si Algo Falla

### Fallo: GPU out of memory

```python
# Editar phase1/config.py
GPU_MEMORY_FRACTION = 0.5  # Reducir a 50%
NUM_WORKERS = 4            # Reducir workers
```

### Fallo: Timeout en workers

```bash
# Ver logs
tail phase1/logs/*stage2*

# Aumentar timeout en config.py
# Luego rerun (se recupera automáticamente desde checkpoint)
python scripts/stage2_parallel_gpu.py --stage 2
```

### Fallo: JSON corrupto

```bash
# Validar
python -m json.tool phase1/output/stage1_indices.json > /dev/null

# Si no valida, ver logs para errores
tail -100 phase1/logs/*.log | grep -i error
```

---

## 📊 Monitoreo en Vivo

En otra terminal, mientras se ejecuta:

```bash
# Ver logs actualizándose
watch -n 1 'tail -20 phase1/logs/*.log'

# Ver CPU/RAM en vivo
watch -n 1 'ps aux | grep stage'

# Ver GPUs
watch -n 1 'nvidia-smi'

# Ver checkpoints siendo creados
watch -n 5 'ls -la phase1/checkpoints/*/checkpoint_*.pkl | tail -5'
```

---

## ✅ Después de completar Phase 1

Una vez que `stage3_knowledge_base.json` está generado:

1. **Validar KB**
   ```bash
   python -m json.tool phase1/output/stage3_knowledge_base.json | head -50
   ```

2. **Contar estadísticas**
   ```bash
   python << 'EOF'
   import json
   data = json.load(open('phase1/output/stage3_knowledge_base.json'))
   print(f"Organismos: {len(data['organisms']):,}")
   print(f"Estrategias: {len(data['strategies']):,}")
   print(f"Archivos: {data['files_processed']:,}")
   EOF
   ```

3. **Ir a Phase 2** (Query Optimizer)
   - Espera carpeta `phase2/` pronto
   - Usará KB de Phase 1 para generar queries mejores

---

## 🎓 Estructura de Fases

```
PHASE 1: KB Extraction ✅ (ACTIVO)
├─ Stage 1: Parse XML (1K test)     ← START HERE
├─ Stage 2: Parallel (50K prod)     ← Next (GPU)
└─ Stage 3: Full Scale (500K final) ← Last (GPU)

PHASE 2: Query Optimizer (Próximo, 4-5 semanas)
├─ Build Query Builder
├─ Test con 50+ cases
└─ Integrate con Pyner_search_v0.1

PHASE 3: Production (Después, 2-4 semanas)
├─ Scale a 1.3M archivos
├─ Evaluate precision
└─ Release v0.3
```

---

## 📖 Documentación Completa

| Archivo | Para |
|---------|------|
| [README.md](README.md) | Intro del proyecto |
| [INDEX.md](INDEX.md) | Quick start 5 min |
| [phase1/README.md](phase1/README.md) | ← **LEER PRIMERO** |
| [planning/DEVELOPMENT_PLAN.md](planning/DEVELOPMENT_PLAN.md) | Detalles técnicos |
| [planning/ROADMAP.md](planning/ROADMAP.md) | Timeline |
| [planning/CONTRIBUTING.md](planning/CONTRIBUTING.md) | Cómo contribuir |

---

## 🎯 Próximo Paso

### Ahora:

```bash
cd phase1
python scripts/stage1_parse_xml.py
```

Esto te dará:
- ✅ Files XML parseados
- ✅ Índices generados
- ✅ Confianza en los scripts
- ⏱️ ~2-5 minutos

### Si todo funciona:

```bash
python scripts/stage2_parallel_gpu.py --stage 2
```

Esto escala a:
- 50K archivos
- 8 workers en paralelo
- 3 GPUs distribuidas
- ⏱️ ~5-10 minutos

### Si stage 2 funciona:

```bash
python scripts/stage2_parallel_gpu.py --stage 3
```

Final:
- 500K archivos
- KB completa
- ⏱️ ~50-100 minutos
- 💾 ~500MB JSON con sinónimos

---

## 🚀 ¡Listoooo!

**El proyecto está totalmente estructurado y listo para ejecutar.**

Comienza:
```bash
cd phase1
python scripts/stage1_parse_xml.py
```

Documentación adicional en: [phase1/README.md](phase1/README.md)

¿Preguntas? Ver logs en `phase1/logs/` o revisar `phase1/config.py` para ajustes avanzados.

---

**Actualizado:** Febrero 6, 2026  
**Status:** 🟢 Listo para ejecutar Phase 1 Stage 1
