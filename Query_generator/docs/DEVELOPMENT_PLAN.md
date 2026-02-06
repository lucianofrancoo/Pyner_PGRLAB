# 📋 Plan de Desarrollo - Pyner PGRLAB v0.3

**Última actualización:** Febrero 6, 2026  
**Estado:** 🔄 En Planificación  
**Objetivo Principal:** Mejorar generación de queries de búsqueda usando Knowledge Base local de NCBI SRA

---

## 🎯 Visión General

El proyecto Pyner actualmente usa un LLM genérico (Qwen2.5) para generar queries NCBI con sinónimos básicos. **El desafío:** generar queries *especializadas* asimilando el conocimiento de **1.3 millones de archivos XML** del NCBI SRA local.

**Solución propuesta:** Extraer, estructurar y usar esa metainformación para crear un modelo contextualizado que entienda:
- Organismos reales estudiados en SRA
- Estrategias de secuenciación existentes
- Atributos de muestras (tejidos, tratamientos, etc.)
- Palabras clave del dominio bioinformático

**Resultado esperado:** Queries 3-5x más precisas → Menos resultados irrelevantes → Análisis más rápido

---

## 📊 Magnitud de los Datos

```
Total de archivos:    1,340,469 archivos XML
Estructura:          
  ├── 445,489 × experiment.xml  (títulos, estrategias, librerías)
  ├── 445,489 × sample.xml      (organismos, atributos)
  └── 445,489 × run.xml         (metadatos de ejecución)

Primeros 5 directorios analizados:
  DRA000001 - DRA000008 ✅ (archivo binario de estructura confirmada)
```

---

## 🚀 Plan de Desarrollo - 3 Fases

### ⭐ FASE 1: Extracción de Knowledge Base (KB)
**Duración:** 1-2 semanas  
**Responsable:** Principal developer  
**Status:** 🟡 Por comenzar

#### Objetivo
Procesar primeros ~100k-500k archivos para extraer:
- Organismos únicos y frecuencia
- Estrategias de secuenciación (LIBRARY_STRATEGY)
- Fuentes de biblioteca (LIBRARY_SOURCE)
- Selecciones de biblioteca (LIBRARY_SELECTION)
- Tipos de tejido/muestra (SAMPLE_ATTRIBUTES)
- Palabras clave de títulos y descripciones

#### Salida Esperada
```json
{
  "organisms": {
    "arabidopsis_thaliana": {
      "freq": 12453,
      "aliases": ["arabidopsis", "Arabidopsis thaliana", "ath"],
      "synonyms": ["small cress", "mouse-ear cress"]
    },
    "solanum_lycopersicum": {...}
  },
  
  "strategies": {
    "RNA-Seq": {
      "freq": 89234,
      "aliases": ["transcriptomics", "transcript-seq", "whole transcriptome"],
      "related_keywords": ["expression", "transcription", "rna"]
    },
    "WGS": {...},
    "AMPLICON": {...}
  },
  
  "treatments_and_conditions": {
    "drought": {
      "freq": 5234,
      "variants": ["water stress", "water deficit", "dehydration", "desiccation"],
      "organism_specific": {
        "arabidopsis": ["drought", "water stress"],
        "rice": ["drought", "water deficit"]
      }
    },
    "nitrogen": {...},
    "temperature": {...}
  },
  
  "tissue_types": {
    "leaf": {"freq": 34521, "aliases": ["leaves", "foliage"]},
    "root": {"freq": 28901, "aliases": ["roots", "rhizome"]},
    ...
  }
}
```

#### Tareas Específicas

**1.1 Crear parser XML robusto**
```python
# Script: extract_metadata.py
- Manejar 100k archivos en paralelo
- Error handling para XMLs corrompidos
- Extraer con regex/xpath campos específicos
- Normalizar nombres (minúsculas, espacios)
- Contar frecuencias
```

**1.2 Procesar sample.xml**
```python
# Información crítica:
- TAXON_ID → SCIENTIFIC_NAME (mapping)
- SAMPLE_ATTRIBUTES → extraer valores de tags
- Buscar patrones: tejido, tratamiento, tiempo, genotipo
- Crear índice invertido (tag → lista de valores únicos)
```

**1.3 Procesar experiment.xml**
```python
# Información crítica:
- TITLE → tokenizar y extraer conceptos
- LIBRARY_STRATEGY, LIBRARY_SOURCE, LIBRARY_SELECTION
- LIBRARY_CONSTRUCTION_PROTOCOL → palabras clave importantes
```

**1.4 Crear índice invertido global**
```python
# Para cada concepto (ej: "drought"):
- Variantes encontradas en SRA
- Frecuencia en cada organismo
- Campos donde aparece (título, atributo, etc.)
- Contexto (¿junto con qué se menciona?)
```

#### Checkpoint 1
- [ ] Parser XML funcional probado en 1000 archivos
- [ ] Índice invertido inicial (100k archivos)
- [ ] JSON de KB exportado
- [ ] Estadísticas principales reportadas

---

### 🤖 FASE 2: Entrenamiento de Modelo Contextualizado
**Duración:** 2-3 semanas  
**Responsable:** ML engineer / Bioinformatician  
**Status:** 🔴 Bloqueado por Fase 1

#### Objetivo
Desarrollar modelo que entienda contextos biológicos y traduzca queries naturales a queries NCBI óptimas.

#### Enfoque 1: Basado en Reglas + KB (Determinístico) ⭐⭐⭐
```python
# RECOMENDADO: Rápido, reproducible, no requiere GPU

class QueryOptimizer:
    def __init__(self, kb):
        self.kb = kb  # Knowledge base de Fase 1
    
    def generate_query(self, user_input):
        """
        "arabidopsis drought" → 
        arabidopsis[Organism] AND 
        (drought OR "water stress" OR "water deficit"[All Fields])
        """
        # 1. Parsear input
        tokens = tokenize(user_input)
        
        # 2. Normalizar contra KB
        normalized = []
        for token in tokens:
            if token in kb["organisms"]:
                normalized.append(("organism", token, kb["organisms"][token]))
            elif token in kb["treatments"]:
                normalized.append(("treatment", token, kb["treatments"][token]))
            # ...
        
        # 3. Construir query booleana
        query_parts = []
        for type_, term, data in normalized:
            if type_ == "organism":
                query_parts.append(f"{data['scientific_name']}[Organism]")
            elif type_ == "treatment":
                variants = " OR ".join(f'"{v}"' for v in data["variants"])
                query_parts.append(f"({variants}[All Fields])")
        
        return " AND ".join(query_parts)
```

#### Enfoque 2: Fine-tune del LLM (Opcional, después)
```python
# Para futuro: si queremos usar LLM local + KB

# Training data:
# Input: "arabidopsis drought"
# Output: "arabidopsis[Organism] AND (drought OR \"water stress\"...)"

# Con ejemplos reales de SRA, el LLM aprenderá mejor
```

#### Tareas Específicas

**2.1 Implementar Query Builder determinístico**
```python
# Script: query_optimizer.py
- Mapeo: token user input → campos NCBI específicos
- Manejo de ambigüedades (ej: "rice" = múltiples especies)
- Expansion de sinónimos automática basada en KB
- Testing contra ejemplos conocidos
```

**2.2 Crear test suite**
```python
# Test cases:
- "arabidopsis drought" → debe expandirse correctamente
- "Solanum lycopersicum nitrogen" → mapear a campos existentes
- "human cancer transcriptomics" → detectar que es incompatible con SRA
- Casos edge: abreviaturas, typos, términos raros
```

**2.3 Integración con Pyner_search_v0.1.py**
```python
# Reemplazar LLM genérico con QueryOptimizer

# Antes:
- LLM genérico genera query genérica

# Después:
- QueryOptimizer usa KB local
- Más preciso, más rápido, sin dependencias de LLM
```

#### Checkpoint 2
- [ ] QueryOptimizer implementado
- [ ] Test suite con 50+ casos
- [ ] Resultados 3-5x más precisos que baseline
- [ ] Documentación de campos NCBI soportados

---

### 📈 FASE 3: Fine-tuning y Optimización
**Duración:** 2-4 semanas  
**Responsable:** Equipo completo (feedback loops)  
**Status:** 🔴 Bloqueado por Fase 2

#### Objetivo
Validar, refinar, publicar y preparar para producción.

#### Tareas Específicas

**3.1 Evaluación exhaustiva**
```python
# Comparar:
- Queries antigas (LLM genérico) vs nuevas (QB + KB)
- Precisión: ¿Cuántos resultados relevantes?
- Recall: ¿Se pierden estudios importantes?
- Velocidad de búsqueda

# Métricas:
- Mean Rank de estudios relevantes
- % de "ruido" en resultados
- Tiempo de búsqueda
```

**3.2 Feedback de usuarios**
```
- Ejecutar Pyner_v0.2.py con nueva query
- Validar: ¿Son los resultados útiles?
- Ajustar KB basado en feedback
- Iterar
```

**3.3 Escalabilidad**
```python
# Procesar todos los 1.3M archivos (no solo 100k)
- Parallelización más agresiva
- Optimización de índices
- Caché distribuida
- Deploy en servidor si es necesario
```

**3.4 Documentación y Release**
```
- Wiki con ejemplos de queries
- API documentation
- Casos de uso de usuario
- Release v0.3 en GitHub
```

#### Checkpoint 3
- [ ] Evaluación completa completada
- [ ] KB procesada al 100% (1.3M archivos)
- [ ] Release v0.3 publicado
- [ ] Documentación de usuario lista

---

## 📂 Estructura de Carpetas (Propuesta)

```
Pyner_PGRLAB/
├── README.md                          # Explicación general ✅
├── DEVELOPMENT_PLAN.md                # Este archivo ← TÚ ESTÁS AQUÍ
│
├── scripts/
│   ├── Pyner_search_v0.1.py           # Query builder (existente)
│   ├── Pyner_v0.1.py                  # Búsqueda simple (existente)
│   ├── Pyner_v0.2.py                  # Búsqueda avanzada (existente)
│   │
│   ├── extract_metadata.py            # 🆕 FASE 1: Parser XML
│   ├── build_kb.py                    # 🆕 FASE 1: Construir Knowledge Base
│   ├── query_optimizer.py             # 🆕 FASE 2: Query Builder determinístico
│   └── evaluate_queries.py            # 🆕 FASE 3: Evaluación de queries
│
├── data/
│   ├── kb/
│   │   ├── organisms_index.json       # 🆕 Output de Fase 1
│   │   ├── strategies_index.json      # 🆕
│   │   ├── treatments_index.json      # 🆕
│   │   ├── tissues_index.json         # 🆕
│   │   └── global_index.json          # 🆕 Índice invertido completo
│   │
│   └── results/
│       ├── baseline_queries.csv       # Queries antiguas (para comparación)
│       └── optimized_queries.csv      # Queries nuevas
│
├── tests/
│   ├── test_extract_metadata.py       # Tests parser XML
│   ├── test_query_optimizer.py        # Tests query builder
│   └── test_cases.json                # Casos de prueba
│
└── docs/
    ├── API.md                         # 🆕 Documentación de QueryOptimizer
    ├── KB_SCHEMA.md                   # 🆕 Esquema del Knowledge Base
    └── CONTRIBUTING.md                # 🆕 Guía para contribuyentes
```

---

## 👥 Roles y Responsabilidades

| Rol | Tareas | Requisitos |
|-----|--------|-----------|
| **Principal Dev** | Fase 1 (Parser XML) | Python, XML/XPath, paralelización |
| **ML Engineer** | Fase 2 (Query Optimizer) | Lógica de búsqueda booleana, algo bioinformatica |
| **QA** | Fase 3 (Testing) | Casos de prueba, evaluación |
| **DevOps** (opcional) | Infraestructura | Docker, servidor (si escalamos) |

---

## 📅 Timeline Estimado

```
FASE 1 (Extracción KB)
├─ Semana 1:  Parser XML + tests básicos
├─ Semana 2:  Procesar 100k archivos + validar índice
└─ Semana 2:  Optimización de rendimiento
                ✅ KB inicial completado

FASE 2 (Query Optimizer)
├─ Semana 3:  Query Builder v1
├─ Semana 4:  Integration + test suite (50+ casos)
└─ Semana 4:  Validación contra baseline LLM
                ✅ QueryOptimizer productivo

FASE 3 (Refinamiento)
├─ Semana 5-6: Scalability (todos los 1.3M archivos)
├─ Semana 6-7: Evaluación exhaustiva + feedback
└─ Semana 8:   Release v0.3 + Documentación
                ✅ Release publico

TOTAL: ~8 semanas (2 meses)
Con parallelización: Podría reducirse a 4-6 semanas
```

---

## 🎓 Cómo Contribuir

### Para Contribuyentes que quieren participar:

**FASE 1 - Extracción KB**
```bash
# Si quieres ayudar:
1. Fork el proyecto
2. Crea rama: git checkout -b feature/extract-metadata
3. Implementa parser para un tipo de XML (sample.xml, experiment.xml)
4. Tests unitarios incluidos
5. Pull request a main con descripción

# Tareas disponibles:
- [ ] Parser de SAMPLE_ATTRIBUTES
- [ ] Parser de LIBRARY_DESCRIPTOR
- [ ] Normalización de nombres de organismos
- [ ] Deduplicación de sinónimos
```

**FASE 2 - Query Optimizer**
```bash
# Si tienes experiencia en bioinformatica:
1. Entender estructura de queries NCBI
2. Implementar mapeos: conceptos usuario → campos NCBI
3. Testing exhaustivo
4. Documentación de lógica de mapeo
```

**FASE 3 - Testing y Docs**
```bash
# Si prefieres QA/documentación:
1. Crear casos de prueba en test_cases.json
2. Evaluar precisión de queries
3. Escribir documentación de usuario
4. Feedback y mejoras
```

### Checklist de Contribución
- [ ] Código comentado y formateado (PEP8)
- [ ] Tests incluidos (mínimo 70% coverage)
- [ ] Documentación clara
- [ ] Commit message descriptivo
- [ ] Referencia issues/PRs relacionados

---

## 🔧 Requisitos Técnicos

### Hardware (Recomendado)
```
Para Fase 1 (procesar 1.3M archivos):
- Mínimo: 16GB RAM, procesador multicore
- Ideal: 32GB RAM, GPU (opcional)

Usuario tiene ya disponible:
✅ GPU 3x NVIDIA RTX 4000 Ada (~24GB VRAM cada una)
✅ 251GB RAM
✅ Procesador multicore
→ PERFECTO para escalar a todos los 1.3M archivos
```

### Software
```
- Python 3.8+
- BioPython (Entrez)
- Ollama (para LLM si lo usamos)
- Pandas (análisis)
- NumPy (operaciones)
```

### Dependencias
```bash
pip install -r requirements-dev.txt
# Includes: pytest, black, flake8 (testing/linting)
```

---

## 📊 Métricas de Éxito (KPIs)

Al final de Fase 3, esperamos:

| Métrica | Baseline | Target |
|---------|----------|--------|
| Precisión de queries | 65% | 85-90% |
| # de resultados "ruido" | -40% | -10-15% |
| Tiempo generación query | 3-5s (LLM) | <500ms (KB) |
| Cobertura de organismos | 10 especies | 500+ especies |
| Documentación | 2 archivos | 5+ + API docs |

---

## 🚨 Riesgos y Mitigación

| Riesgo | Probabilidad | Impacto | Mitigación |
|--------|--------------|--------|-----------|
| Memory overflow (1.3M files) | Media | Alto | Procesar en chunks, usar generadores |
| XML corrompidos | Alta | Bajo | Error handling robusto + logging |
| Lógica booleana incorrecta | Media | Alto | Test suite completo + validación |
| Performance (búsqueda lenta) | Baja | Medio | Indexación, caché, optimización |

---

## 📞 Contacto y Preguntas

- **Issues:** [GitHub Issues](https://github.com/lucianofrancoo/Pyner_PGRLAB/issues)
- **Discussions:** [GitHub Discussions](https://github.com/lucianofrancoo/Pyner_PGRLAB/discussions)
- **Email:** lucianofranco.a@gmail.com

---

## 📝 Cambios y Actualizaciones

| Fecha | Quién | Cambio |
|-------|-------|--------|
| 2026-02-06 | GitHub Copilot | Plan inicial creado |
| - | - | - |

---

**Estado General:** 🟡 Listo para Fase 1  
**Última review:** 2026-02-06  
**Próxima milestone:** Fase 1, Checkpoint 1 ✅
