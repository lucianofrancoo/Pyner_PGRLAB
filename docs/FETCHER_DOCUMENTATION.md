# PYNER - Guía Completa del Sistema Integrado

## 📋 Índice
1. [Visión General](#visión-general)
2. [Arquitectura del Sistema](#arquitectura-del-sistema)
3. [Inicio Rápido](#inicio-rápido)
4. [Uso Detallado](#uso-detallado)
5. [Estructura de Datos](#estructura-de-datos)
6. [Módulos del Sistema](#módulos-del-sistema)
7. [Resolución de Problemas](#resolución-de-problemas)

---

## 🎯 Visión General

PYNER es un sistema integrado que permite buscar información científica utilizando lenguaje natural, obteniendo:
- **Proyectos genómicos** (BioProject)
- **Datos experimentales** (SRA)
- **Publicaciones científicas** (PubMed)

Todo en un solo flujo automatizado.

### Workflow Completo

```
Lenguaje Natural → Boolean Query → Base de Datos → SRA → PubMed → CSV
     ↓                    ↓              ↓            ↓       ↓        ↓
"Arabidopsis    →  "Arabidopsis  →  BioProject  →  42 exp  → 2 papers → results.csv
 phosphate          AND              PRJNA...       90 exp    0 papers
 stress"            phosphate"       PRJDB...       
```

---

## 🏗️ Arquitectura del Sistema

### Componentes Principales

```
PYNER/
├── Query_generator/phases/phase3/     # IA: Lenguaje natural → Boolean
│   └── api/main.py
│
├── Fetcher_NCBI/                      # Motor de búsqueda y linking
│   ├── boolean_fetcher_integrated.py  # Orquestador principal
│   ├── ncbi_fetcher_sra_fixed.py     # Fetch de SRA
│   ├── ncbi_linkout.py               # Linking a PubMed
│   └── config.py                      # API keys y configuración
│
└── test_fetcher_integrator.sh        # Script principal (interfaz)
```

### Flujo de Datos

#### Opción 1: BioProject (RECOMENDADO)
```
┌─────────────────────────────────────────────────────────────┐
│ 1. LENGUAJE NATURAL                                         │
│    "Arabidopsis under phosphate deficiency"                 │
└────────────┬────────────────────────────────────────────────┘
             ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. IA GENERA BOOLEAN QUERY                                  │
│    "Arabidopsis AND phosphate AND deficiency"               │
└────────────┬────────────────────────────────────────────────┘
             ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. BÚSQUEDA EN BIOPROJECT                                   │
│    Encuentra: PRJDB35626, PRJNA1179470, etc.              │
└────────────┬────────────────────────────────────────────────┘
             ↓
┌─────────────────────────────────────────────────────────────┐
│ 4. PARA CADA BIOPROJECT:                                    │
│                                                             │
│    A) BUSCAR EN SRA                                        │
│       ✓ Extrae todos los SRX (experiments)                │
│       ✓ Extrae todos los SAMN (biosamples)                │
│       ✓ Metadata completa de cada experiment              │
│                                                             │
│    B) CASCADE SEARCH EN PUBMED                             │
│       Level 1: Buscar BIOPROJECT directo                  │
│                Query: "PRJNA1179470"[All Fields]          │
│                                                             │
│       Level 2: Buscar cada BIOSAMPLE                      │
│                Query: "SAMN44494209"[All Fields]          │
│                       "SAMN44494208"[All Fields]          │
│                                                             │
│       Level 3: Buscar cada SRA ACCESSION                  │
│                Query: "SRX26886995"[All Fields]           │
│                       "SRX26886994"[All Fields]           │
│                                                             │
│    C) SI NO ENCUENTRA PAPERS → Marcar "NA"                │
└────────────┬────────────────────────────────────────────────┘
             ↓
┌─────────────────────────────────────────────────────────────┐
│ 5. OUTPUT CSV + JSON                                        │
│    - bioproject, title, organism, description              │
│    - sra_experiments_count, biosamples_count               │
│    - publications_found, dois, pmids                       │
│    - "NA" cuando no hay papers                             │
└─────────────────────────────────────────────────────────────┘
```

#### Opción 2: PubMed (Búsqueda Directa)
```
Lenguaje Natural → Boolean Query → PubMed directo → Papers → CSV
```

**Workflow PubMed:**
```
┌─────────────────────────────────────────────────────────────┐
│ 1. LENGUAJE NATURAL                                         │
│    "Arabidopsis phosphate deficiency"                       │
└────────────┬────────────────────────────────────────────────┘
             ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. IA GENERA BOOLEAN QUERY                                  │
│    "Arabidopsis AND phosphate AND deficiency"               │
└────────────┬────────────────────────────────────────────────┘
             ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. BÚSQUEDA DIRECTA EN PUBMED                              │
│    Query: ("Arabidopsis thaliana"[Organism] OR             │
│           "arabidopsis"[All Fields]) AND                    │
│           "phosphate deficiency"[All Fields]                │
│    Encontrados: 165 papers                                  │
└────────────┬────────────────────────────────────────────────┘
             ↓
┌─────────────────────────────────────────────────────────────┐
│ 4. EXTRACCIÓN DE METADATA                                  │
│    Para cada paper (max según usuario):                    │
│    ✓ PMID, Title, Year                                     │
│    ✓ Authors, DOI, URL                                     │
│    ✓ Abstract                                              │
└────────────┬────────────────────────────────────────────────┘
             ↓
┌─────────────────────────────────────────────────────────────┐
│ 5. OUTPUT CSV + JSON                                        │
│    - pmid, title, year, authors                            │
│    - doi, url, abstract                                    │
│    - Número total encontrado vs extraído                   │
└─────────────────────────────────────────────────────────────┘
```

**Ventajas de PubMed directo:**
- ✅ Rápido (solo búsqueda en una DB)
- ✅ Obtiene papers directamente
- ✅ Ideal para revisiones bibliográficas
- ⚠️ No incluye datos SRA/BioProject

---

## 🚀 Inicio Rápido

### Requisitos Previos

```bash
# 1. Python 3.8+
python --version

# 2. Dependencias instaladas
cd /home/lahumada/disco1/Pyner_PGRLAB/Fetcher_NCBI
pip install -r requirements.txt

# 3. API Key configurada (verificar)
grep NCBI_API_KEY config.py
```

### Ejecución Simple

```bash
cd /home/lahumada/disco1/Pyner_PGRLAB
bash test_fetcher_integrator.sh
```

**Ejemplo de sesión:**
```
================================================================
         PYNER - INTEGRATED BOOLEAN FETCHER
================================================================

[1/5] Ingresa tu consulta en lenguaje natural:
> Arabidopsis under phosphate stress

[2/5] Selecciona la base de datos de búsqueda:
  1) BioProject (recomendado para proyectos genómicos)
  2) PubMed (búsqueda directa en publicaciones)
Selección [1-2]: 1
  ✓ Seleccionado: BioProject

[3/5] Generando query booleano con IA...

NCBI Query:
 Arabidopsis AND phosphate AND stress

Query generado: Arabidopsis AND phosphate AND stress
¿Continuar con este query? [S/n]: S

[4/5] ¿Cuántos BioProjects/resultados quieres procesar?
Número máximo (default 20): 10

[5/5] Ejecutando workflow integrado BioProject → SRA → PubMed...

[Processing logs...]

✓ WORKFLOW COMPLETADO
Resultados guardados en:
  📄 CSV:  results_20260212_143000.csv
  📄 JSON: results_20260212_143000.json

Estadísticas:
  Total BioProjects: 10
  Con publicaciones: 2
  Sin publicaciones: 8
```

---

## 📘 Uso Detallado

### Modo Interactivo (Recomendado)

```bash
bash test_fetcher_integrator.sh
```

**Ventajas:**
- Interfaz guiada paso a paso
- Selección de base de datos
- Confirmación del query generado
- Control sobre número de resultados

### Modo Directo (Avanzado)

Si quieres saltar la IA y usar directamente un boolean query:

```bash
cd Fetcher_NCBI

# Búsqueda simple
python boolean_fetcher_integrated.py "Arabidopsis phosphate" --max 20 --output-csv results.csv

# Búsqueda avanzada con operadores
python boolean_fetcher_integrated.py \
    "Arabidopsis AND (phosphate OR phosphorus) AND stress" \
    --max 50 \
    --output-csv results.csv \
    --output-json results.json
```

### Parámetros del Fetcher

| Parámetro | Descripción | Ejemplo | Default |
|-----------|-------------|---------|---------|
| `query` | Boolean query | `"Arabidopsis phosphate"` | - |
| `--max` | Máximo de BioProjects | `50` | 50 |
| `--output-csv` | Archivo CSV | `results.csv` | Auto-generado |
| `--output-json` | Archivo JSON | `results.json` | No generado |

---

## 📊 Estructura de Datos

### CSV Output - BioProject Workflow

**Columnas del CSV:**

```csv
bioproject,title,submission_date,organism,project_type,description,
sra_experiments_count,biosamples_count,sra_accessions_count,
sra_experiments,biosamples,sra_accessions,
publications_found,search_method,papers_summary,dois,pmids
```

#### Descripción de Campos (BioProject)

| Campo | Tipo | Descripción | Ejemplo |
|-------|------|-------------|---------|
| `bioproject` | String | ID del BioProject | PRJNA1179470 |
| `title` | String | Título del proyecto | "Transcriptomic Profiling..." |
| `submission_date` | Date | Fecha de registro | 2024-12-01 |
| `organism` | String | Organismo(s) | Arabidopsis thaliana |
| `project_type` | String | Tipo de proyecto | Transcriptome |
| `description` | Text | Descripción completa | "In Arabidopsis..." |
| `sra_experiments_count` | Integer | Número de experiments SRA | 42 |
| `biosamples_count` | Integer | Número de biosamples únicos | 42 |
| `sra_accessions_count` | Integer | Número de accessions SRA | 42 |
| `sra_experiments` | String | IDs de experimentos (separados por `;`) | "SRX31557536; SRX31557537; ..." |
| `biosamples` | String | IDs de biosamples (separados por `;`) | "SAMN54118006; SAMN54118007; ..." |
| `sra_accessions` | String | IDs de accessions (separados por `;`) | "SRX31557536; SRX31557537; ..." |
| `publications_found` | Integer | Papers encontrados | 2 o 0 |
| `search_method` | String | Método que encontró papers | direct/biosamples/sra_accessions/NA |
| `papers_summary` | String | Resumen de resultados | "2 paper(s) - direct" o NA |
| `dois` | String | DOIs separados por `;` | "10.1234/abc; 10.5678/def" o NA |
| `pmids` | String | PMIDs separados por `;` | "35123456; 35123457" o NA |

### Ejemplos de Registros

#### Caso 1: BioProject CON publicaciones
```csv
PRJNA1329195,"Differentially expressed microRNAs between WT and AtC3H3 OX",2025-09-16,Arabidopsis thaliana,Transcriptome,"CCCH-type zinc finger proteins...",12,12,12,2,direct,"2 paper(s) - direct","10.1234/abc; 10.5678/def","38123456; 38123457"
```

#### Caso 2: BioProject SIN publicaciones (más común)
```csv
PRJNA1179470,"Transcriptomic response to drought",2024-12-01,Arabidopsis thaliana,Transcriptome,"Study of drought stress...",42,42,42,0,NA,NA,NA,NA
```

---

### CSV Output - PubMed Direct Search

**Columnas del CSV:**
```csv
pmid,title,year,journal,publication_type,authors,doi,pmcid,url,abstract,fetched_at
```

**Descripción de campos:**
- `pmid`: PubMed ID único
- `title`: Título de la publicación
- `year`: Año de publicación
- `journal`: Nombre de la revista
- `publication_type`: Tipo de artículo (Journal Article, Review, Meta-Analysis, etc.)
- `authors`: Primeros 3 autores (separados por ";")
- `doi`: Digital Object Identifier
- `pmcid`: PubMed Central ID (si está disponible en texto completo, sino "NA")
- `url`: Enlace a PubMed
- `abstract`: **Resumen completo** (sin truncar)
- `fetched_at`: Timestamp de descarga

**Ejemplo:**
```csv
35099557,"The tomato OST1-VOZ1 module...",2022,"The Plant cell","Journal Article; Research Support, Non-U.S. Gov't","Chong L; Xu R; Huang P",10.1093/plcell/koac026,PMC9048945,https://www.ncbi.nlm.nih.gov/pubmed/35099557,"Flowering is a critical...",2026-02-12T14:23:38
37373193,"Molecular Mechanisms...",2023,"International journal of molecular sciences","Journal Article; Review","Liang W; Ma X; Wan P",10.3390/ijms241310637,PMC10298849,https://www.ncbi.nlm.nih.gov/pubmed/37373193,"Water is one of the most critical factors...",2026-02-12T14:24:36
```

**Tipos de publicación comunes:**
- `Journal Article` - Artículo de investigación
- `Review` - Revisión bibliográfica
- `Meta-Analysis` - Meta-análisis
- `Research Support, Non-U.S. Gov't` - Con financiamiento no-US

**Nota sobre PMCID:** El PMCID permite acceder al texto completo gratuito en PubMed Central cuando está disponible:
- Con PMCID: `https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9048945/`
- Artículos en PubMed Central son de acceso abierto

#### Comparación BioProject vs PubMed

| Feature | BioProject | PubMed Direct |
|---------|------------|---------------|
| Tiempo | 10-60 min | 1-5 min |
| Datos SRA | ✅ Sí | ❌ No |
| Papers | Via cascade | Directos |
| Metadata | Journal/Type ❌ | Journal/Type ✅ |
| Uso ideal | Proyectos genómicos | Revisión bibliográfica |

---

### JSON Output

```json
{
  "metadata": {
    "total_results": 10,
    "date": "2026-02-12T14:30:00",
    "with_publications": 2,
    "without_publications": 8
  },
  "results": [
    {
      "bioproject": "PRJNA1179470",
      "title": "...",
      "sra_experiments_count": 42,
      "publications_found": 0,
      "search_method": "NA",
      ...
    }
  ]
}
```

---

## 🔧 Módulos del Sistema

### 1. Query Generator (Phase 3)

**Ubicación:** `Query_generator/phases/phase3/api/main.py`

**Función:** Convierte lenguaje natural en boolean query usando IA (Ollama)

**Ejemplo:**
```python
# Input
"plants under phosphate stress"

# Output
"plants AND phosphate AND stress"
```

**Uso directo:**
```bash
cd Query_generator/phases/phase3
python api/main.py --quick "Arabidopsis phosphate stress"
```

### 2. Boolean Fetcher Integrated

**Ubicación:** `Fetcher_NCBI/boolean_fetcher_integrated.py`

**Clase principal:** `BooleanFetcherIntegrated`

**Métodos clave:**

```python
class BooleanFetcherIntegrated:
    def search_bioproject_boolean(query, retmax=100)
        # Busca BioProjects con query booleano
        
    def fetch_sra_for_bioproject(bioproject_id)
        # Extrae SRA experiments y biosamples
        
    def search_pubmed_publications(bioproject, biosamples, sra_accs, cascade=True)
        # Cascade search en PubMed (3 niveles)
        
    def run_workflow(query, max_bioproject=50)
        # Orquesta el flujo completo
```

**Uso programático:**
```python
from boolean_fetcher_integrated import BooleanFetcherIntegrated

fetcher = BooleanFetcherIntegrated()
results = fetcher.run_workflow("Arabidopsis phosphate", max_bioproject=20)
fetcher.save_results_csv(Path("results.csv"))
```

### 3. SRA Fetcher

**Ubicación:** `Fetcher_NCBI/ncbi_fetcher_sra_fixed.py`

**Clase principal:** `SRAFetcher`

**Función:** Extrae metadata completa de SRA experiments

**Métodos:**
```python
class SRAFetcher:
    def search_sra_by_bioproject(bioproject_id)
        # Encuentra IDs de SRA experiments
        
    def fetch_sra_metadata(sra_id)
        # Extrae 12 campos de metadata
        
    def fetch_all_by_bioproject(bioproject_id, max_per_bioproject=None)
        # Fetch completo de todos los experiments
```

**Metadata extraída (12 campos):**
- exp_accession, study_accession, biosample, run_accession
- organism, title, library_strategy, library_source
- library_layout, platform, instrument, total_runs

### 4. PubMed Linkout

**Ubicación:** `Fetcher_NCBI/ncbi_linkout.py`

**Clase principal:** `LinkoutFetcher`

**Función:** Vincula datos SRA con publicaciones en PubMed

**Estrategia Cascade (3 niveles):**

```python
class LinkoutFetcher:
    # Level 1: Direct BioProject
    def search_publications_for_bioproject(bioproject_id)
    
    # Level 2: Via BioSamples
    def search_publications_for_biosamples(biosamples_list)
    
    # Level 3: Via SRA Accessions
    def search_publications_for_sra_accessions(sra_accessions_list)
```

**Metadata de publicaciones extraída:**
- pmid, title, abstract, authors
- year, doi, url

---

## 🔍 Cascade Search: Explicación Detallada

### ¿Por qué Cascade?

Muchos BioProjects no aparecen directamente en PubMed porque:
1. Los datos se suben **antes** de publicar el paper
2. El paper cita los **BioSamples** o **SRA accessions**, no el BioProject
3. Múltiples papers pueden usar el mismo dataset

### Estrategia de Búsqueda

```
Level 1: BIOPROJECT DIRECTO
────────────────────────────
Query: "PRJNA1179470"[All Fields]
├─ ✓ Encontrado → STOP (retorna papers)
└─ ✗ No encontrado → CASCADE a Level 2

Level 2: BIOSAMPLES
────────────────────
Para cada biosample (máx 10):
  Query: "SAMN44494209"[All Fields]
  Query: "SAMN44494208"[All Fields]
  ...
├─ ✓ Alguno encontrado → STOP (retorna papers)
└─ ✗ Ninguno encontrado → CASCADE a Level 3

Level 3: SRA ACCESSIONS
────────────────────────
Para cada SRA acc (máx 10):
  Query: "SRX26886995"[All Fields]
  Query: "SRX26886994"[All Fields]
  ...
├─ ✓ Alguno encontrado → STOP (retorna papers)
└─ ✗ Ninguno encontrado → Marcar "NA"
```

### Ejemplo Real

**BioProject:** PRJNA1179470

**Resultado:**
```
2026-02-12 12:43:02 - INFO - Level 1: Searching BIOPROJECT directly...
2026-02-12 12:43:02 - INFO - Query: "PRJNA1179470"[All Fields]
2026-02-12 12:43:03 - INFO - Found 0 PubMed records
2026-02-12 12:43:03 - INFO - ✗ No results for direct bioproject search

2026-02-12 12:43:03 - INFO - Level 2: Searching 42 BIOSAMPLES...
2026-02-12 12:43:04 - INFO - Query: "SAMN44494209"[All Fields]
2026-02-12 12:43:04 - INFO - × No PubMed records
[... búsqueda de otros biosamples ...]
2026-02-12 12:43:10 - INFO - ✗ No biosamples found in PubMed

2026-02-12 12:43:10 - INFO - Level 3: Searching 42 SRA ACCESSIONS...
2026-02-12 12:43:11 - INFO - Query: "SRX26886995"[All Fields]
2026-02-12 12:43:11 - INFO - × No PubMed records
[... búsqueda de otros SRA ...]
2026-02-12 12:43:20 - INFO - ✗ No SRA accessions found in PubMed

2026-02-12 12:43:20 - INFO - ⊘ No publications found
```

**Resultado final en CSV:**
```csv
PRJNA1179470,...,...,42,42,42,0,NA,NA,NA,NA
```

---

## ⚙️ Configuración

### API Keys

**Archivo:** `Fetcher_NCBI/config.py`

```python
NCBI_EMAIL = "your.email@example.com"
NCBI_API_KEY = "your_api_key_here"  # Obtener en NCBI
RATE_LIMIT = 0.1  # 10 req/sec con API key
```

**Obtener API Key:**
1. Ir a https://www.ncbi.nlm.nih.gov/account/
2. Settings → API Key Management
3. Create new key
4. Copiar y pegar en `config.py`

**Ventajas del API Key:**
- Sin key: 3 requests/segundo
- Con key: **10 requests/segundo** ✅

### Rate Limiting

El sistema respeta automáticamente los límites de NCBI:
```python
time.sleep(RATE_LIMIT)  # 0.1 seg = 10 req/sec
```

---

## 📈 Performance y Optimización

### Tiempos Estimados

| Operación | Tiempo (1 BioProject) | Tiempo (50 BioProjects) |
|-----------|----------------------|-------------------------|
| Boolean search | 1-2 seg | 1-2 seg |
| Fetch SRA (40 exp) | 40 seg | 33 min |
| PubMed cascade | 10-30 seg | 8-25 min |
| **TOTAL** | **~1 min** | **~40-60 min** |

### Recomendaciones

1. **Para testing:** Usar `--max 5` o `--max 10`
2. **Para producción:** Usar `--max 50` (límite recomendado)
3. **Para exploración:** Empezar con 3-5, revisar resultados, luego escalar

### Límites del Sistema

| Parámetro | Límite | Razón |
|-----------|--------|-------|
| BioProjects por búsqueda | 100 (recomendado 50) | Rate limiting |
| BioSamples en cascade | 10 primeros | Evitar sobrecarga |
| SRA accessions en cascade | 10 primeros | Evitar sobrecarga |
| Papers por hit | 3-5 | Suficiente para linking |

---

## 🐛 Resolución de Problemas

### Error: "No BioProjects found"

**Síntoma:**
```
✓ Found 0 BioProjects
WARNING - No BioProjects found
```

**Causa:** Query demasiado específico o términos incorrectos

**Solución:**
```bash
# En vez de:
"Arabidopsis thaliana Col-0 ecotype phosphate starvation 7 days"

# Usar:
"Arabidopsis phosphate"
```

### Error: "No SRA experiments found"

**Síntoma:**
```
⚠️ Found 0 SRA experiments
```

**Causa:** BioProject sin datos SRA asociados (puede tener solo BioSamples o Assembly)

**Impacto:** El sistema continúa con búsqueda directa de BioProject en PubMed

**No requiere acción:** Es normal, algunos BioProjects no tienen SRA

### Error: "Rate limit exceeded"

**Síntoma:**
```
HTTP Error 429: Too Many Requests
```

**Causa:** Demasiadas requests sin API key o problema temporal

**Solución:**
1. Verificar API key en `config.py`
2. Esperar 1-2 minutos y reintentar
3. Reducir `--max` a un número menor

### Error: "Connection timeout"

**Síntoma:**
```
URLError: <urlopen error [Errno 110] Connection timed out>
```

**Causa:** Problema de red o NCBI temporalmente no disponible

**Solución:**
1. Verificar conexión a internet
2. Reintentar en 5-10 minutos
3. Revisar status de NCBI: https://www.ncbi.nlm.nih.gov/

### Logs Detallados

**Ubicación:** `Fetcher_NCBI/logs/`

Para debugging:
```bash
cd Fetcher_NCBI
tail -f logs/fetcher_*.log
```

---

## 📚 Ejemplos de Uso

### Ejemplo 1: Búsqueda Simple

```bash
bash test_fetcher_integrator.sh
```

**Input:**
```
Arabidopsis phosphate deficiency
```

**Steps:**
1. IA genera: `"Arabidopsis AND phosphate AND deficiency"`
2. Busca en BioProject: 15 encontrados
3. Para cada uno: extrae SRA + busca papers
4. Genera CSV con resultados

### Ejemplo 2: Búsqueda Avanzada

```bash
cd Fetcher_NCBI
python boolean_fetcher_integrated.py \
    "Arabidopsis AND (drought OR water) AND stress AND response" \
    --max 30 \
    --output-csv drought_results.csv \
    --output-json drought_results.json
```

### Ejemplo 3: Análisis Post-Procesamiento

**En Python:**
```python
import pandas as pd

# Cargar resultados
df = pd.read_csv('results_20260212_143000.csv')

# Filtrar solo con publicaciones
with_papers = df[df['publications_found'] > 0]
print(f"BioProjects con papers: {len(with_papers)}")

# Estadísticas de SRA
print(f"Promedio de experiments: {df['sra_experiments_count'].mean():.1f}")
print(f"Total experiments: {df['sra_experiments_count'].sum()}")

# Exportar solo con papers
with_papers.to_csv('only_published.csv', index=False)
```

### Ejemplo 4: Búsqueda por Organismo Específico

```bash
python boolean_fetcher_integrated.py \
    "Medicago truncatula AND nitrogen" \
    --max 20 \
    --output-csv medicago_nitrogen.csv
```

---

## 🎓 Mejores Prácticas

### 1. Diseño de Queries

✅ **Buenos queries:**
- `"Arabidopsis phosphate"`
- `"drought AND stress"`
- `"nitrogen AND deficiency"`

❌ **Queries problemáticos:**
- `"plants"` (demasiado amplio)
- `"AT5G12345"` (demasiado específico)
- `"Arabidopsis thaliana Columbia ecotype grown in..."` (demasiado detallado)

### 2. Número de Resultados

| Objetivo | --max recomendado |
|----------|-------------------|
| Testing | 3-5 |
| Exploración | 10-20 |
| Análisis completo | 50 |
| Dataset grande | 100 (ejecutar overnight) |

### 3. Interpretación de Resultados

**publications_found = 0 (NA):**
- ✅ **Normal:** Muchos datasets se publican antes del paper
- ✅ **Esperado:** Estudios en progreso
- ⚠️ **Revisar:** Si TODOS son 0, verificar query

**publications_found > 0:**
- ✅ **Excelente:** Linking exitoso
- 📄 **Acción:** Revisar DOIs y PMIDs

### 4. Workflow Recomendado

```bash
# 1. Test inicial (5 min)
bash test_fetcher_integrator.sh
# Input: "Arabidopsis phosphate"
# Max: 5

# 2. Revisar resultados
head -20 results_*.csv

# 3. Si resultados OK, escalar
cd Fetcher_NCBI
python boolean_fetcher_integrated.py \
    "Arabidopsis phosphate" \
    --max 50 \
    --output-csv full_results.csv

# 4. Análisis
python analyze_results.py full_results.csv
```

---

## 📞 Soporte y Contacto

### Recursos

- **Documentación NCBI:** https://www.ncbi.nlm.nih.gov/home/documentation/
- **BioPython:** https://biopython.org/
- **E-utilities:** https://www.ncbi.nlm.nih.gov/books/NBK25501/

### Mantenimiento

**Logs del sistema:**
```bash
cd Fetcher_NCBI/logs
ls -lth  # Ver logs más recientes
```

**Limpiar cache:**
```bash
cd Fetcher_NCBI
rm -rf __pycache__
rm -f *.pyc
```

---

## 🔄 Actualizaciones y Versionado

**Versión actual:** 1.0.0  
**Fecha:** 2026-02-12  
**Estado:** Producción

### Changelog

**v1.0.0 (2026-02-12)**
- ✅ Workflow completo integrado
- ✅ Cascade PubMed search (3 niveles)
- ✅ Interfaz interactiva
- ✅ Selección de base de datos
- ✅ Export CSV + JSON

---

## 📄 Licencia

PGRLAB Internal Use

---

**Última actualización:** 2026-02-12  
**Mantenedor:** PGRLAB Team  
**Versión:** 1.0.0
# 🎯 Quick Reference: SRA Hierarchy Refactoring

## What Changed?

### ✅ Before (Wrong):
- Output field: `sra_accessions` contained **SRX*** codes (experiments)
- Didn't preserve hierarchy properly
- Misleading about what data was actually stored

### ✅ After (Correct):  
- Output field: `sra_runs` contains **SRR*** codes (actual sequence runs)
- Properly preserves: BioProject → BioSample → Experiment → Run
- Clear and accurate naming

---

## The Hierarchy Explained

```
┌─ BioProject (PRJNA*)
│  └─ BioSample (SAMN*)
│     └─ Experiment (SRX*)
│        └─ Run (SRR*)  ← This is the actual sequence data!
```

**In CSV Output:**
- `sra_experiments` → SRX* codes (metadata)
- `biosamples` → SAMN* codes (samples)
- `sra_runs` → SRR* codes (actual data) ✅

---

## Code Changes Summary

| Component | Old | New | Why |
|-----------|-----|-----|-----|
| Return Type | `(exp, bs_set, srx_set)` | `(exp, bs_dict, srr_list)` | Proper hierarchy |
| Accessions | SRX codes | SRR codes | Real data not metadata |
| BioSamples | Set | Dict | Can store titles/experiments |
| CSV Field | `sra_accessions_count` | `sra_runs_count` | Accurate naming |
| CSV Field | `sra_accessions` | `sra_runs` | Correct codes |

---

## Files Modified

1. **Fetcher_NCBI/boolean_fetcher_integrated.py**
   - `fetch_sra_for_bioproject()` - Returns new structure
   - `search_pubmed_publications()` - Accepts new types
   - `process_bioproject()` - Unpacks and stores new structure
   - `save_results_csv()` - Uses correct field names

2. **README.md**
   - Updated CSV format documentation
   - Added hierarchy explanation
   - Updated examples with SRR codes

3. **New Documentation**
   - `REFACTORING_SUMMARY.md` - Full technical details
   - `HIERARCHICAL_STRUCTURE_CHANGES.txt` - Final report
   - `test_hierarchical_refactoring.py` - Test suite

---

## Testing

✅ All tests passed:
```python
✓ Experiments returned as list
✓ BioSamples returned as dict with experiment titles
✓ SRA Runs returned as list of SRR* codes
✓ Hierarchy preserved end-to-end
✓ CSV output correct
✓ process_bioproject() works with new structure
```

**Run tests:**
```bash
python3 test_hierarchical_refactoring.py
```

---

## Impact on CSV Output

### Example Before (Incorrect):
```
PRJNA1381306,12,12,12,"SRX31557547;...","SAMN54118015;...","SRX31557547;..."
                                                              ↑ Wrong! These are experiments
```

### Example After (Correct):
```
PRJNA1381306,12,12,12,"SRX31557547;...","SAMN54118015;...","SRR36541090;..."
                                                              ↑ Correct! These are runs with actual sequence data
```

---

## Key Takeaways

| Concept | Code | What It Is |
|---------|------|-----------|
| Project ID | PRJNA* | Umbrella for everything |
| Sample ID | SAMN* | Biological sample |
| Experiment ID | SRX* | Sequencing experiment metadata |
| Run ID | SRR* | **Actual sequence data** ← This is what matters! |

---

## For Next Use

When you run Option 2 (BioProject search) again:
- ✅ It will fetch the correct SRR run codes
- ✅ It will store SRX experiment metadata
- ✅ It will reference SAMN biosample codes
- ✅ CSV will have accurate, meaningful data

---

## Questions?

See detailed documentation in:
- `REFACTORING_SUMMARY.md` - Full technical reference
- `README.md` - User-facing documentation
- `test_hierarchical_refactoring.py` - Working code examples

✅ **Status: Ready for Production**
# Hierarchical SRA Data Structure Refactoring - Summary

## Overview
Successfully refactored the SRA data structure in `boolean_fetcher_integrated.py` to properly model the biological hierarchy: **BioProject → BioSample → Experiment (SRX) → Run (SRR)**

**Status**: ✅ COMPLETED AND TESTED

---

## What Changed

### 1. **Data Structure Changes**

#### Before (Incorrect):
```python
def fetch_sra_for_bioproject(self, bioproject_id: str):
    # Returned separate sets
    return experiments, biosamples_set, sra_accessions_set
    # Where sra_accessions was actually SRX (experiment) codes, not SRR (run) codes
```

#### After (Correct):
```python
def fetch_sra_for_bioproject(self, bioproject_id: str):
    # Returns properly structured data
    return experiments, biosamples_dict, sra_runs
    
    # biosamples_dict = {
    #     'SAMN54118006': {
    #         'samples': set(),
    #         'experiment_titles': ['Experiment Title 1', 'Experiment Title 2']
    #     },
    #     ...
    # }
    
    # sra_runs = ['SRR36541090', 'SRR36541091', ...]  (actual run codes)
```

---

## Biological Hierarchy Clarified

```
BioProject (e.g., PRJNA1381306)
  ├── BioSample 1 (e.g., SAMN54118006)
  │   ├── Experiment 1 (e.g., SRX31557536)
  │   │   ├── Run 1 (e.g., SRR36541090)
  │   │   ├── Run 2 (e.g., SRR36541091)
  │   │   └── ...
  │   └── Experiment 2 (SRX...)
  │       └── ...
  ├── BioSample 2 (SAMN...)
  │   └── ...
  └── ...
```

**Key Relationships:**
- BioSample (SAMN*) ← links to → Experiment (SRX*) via `exp['biosample']`
- Experiment (SRX*) ← contains → Runs (SRR*) via `exp['runs']` list
- Each BioSample can have multiple experiments, each experiment has multiple runs

---

## Files Modified

### 1. `boolean_fetcher_integrated.py`

#### Method `fetch_sra_for_bioproject()` (Lines 154-198)
**Changes:**
- ✅ Returns `(experiments, biosamples_dict, sra_runs)` instead of `(experiments, biosamples_set, sra_accessions_set)`
- ✅ `biosamples_dict` is now a dict with structure: `{biosample_id: {'samples': set(), 'experiment_titles': []}}`
- ✅ `sra_runs` now contains **SRR codes** (actual sequencing runs) not SRX codes (experiments)
- ✅ Properly associates experiment titles with their biosamples

#### Method `search_pubmed_publications()` (Lines 201-248)
**Changes:**
- ✅ Updated parameters: `biosamples_dict` (dict) instead of `biosamples` (set)
- ✅ Updated parameters: `sra_runs` (list) instead of `sra_accessions` (set)
- ✅ Uses `list(biosamples_dict.keys())` to extract biosample IDs for searching
- ✅ Uses `sra_runs` list directly for searching
- ✅ Updated log message to reference "SRA runs" instead of "SRA accessions"

#### Method `process_bioproject()` (Lines 280-308)
**Changes:**
- ✅ Updated unpacking: `experiments, biosamples_dict, sra_runs = self.fetch_sra_for_bioproject(bioproject_id)`
- ✅ Updated counts: `biosamples_count` now reflects unique SAMN codes (not sets)
- ✅ Added `sra_runs_count` with unique SRR codes
- ✅ Updated result storage: `result['biosamples'] = sorted(list(biosamples_dict.keys()))`
- ✅ Added `result['sra_runs'] = sorted(list(set(sra_runs)))`
- ✅ Passes correct data types to `search_pubmed_publications()`

#### Method `save_results_csv()` (Lines 405-420)
**Changes:**
- ✅ Updated fieldnames: Changed `sra_accessions_count` → `sra_runs_count`
- ✅ Updated fieldnames: Changed `sra_accessions` → `sra_runs`
- ✅ CSV now outputs actual run codes (SRR*) not experiment codes (SRX*)

### 2. `README.md` (Lines 59-77)
**Changes:**
- ✅ Updated CSV header to show `sra_runs_count` and `sra_runs` fields
- ✅ Updated example to show SRR codes instead of SRX codes
- ✅ Added explanation section "Notas sobre la jerarquía" clarifying:
  - sra_experiments = SRX* codes
  - biosamples = SAMN* codes  
  - sra_runs = SRR* codes (the actual sequencing data)
- ✅ Documented hierarchy: BioProject → BioSample → Experiment (SRX) → Run (SRR)

---

## Testing & Validation

### Test Created: `test_hierarchical_refactoring.py`
**Tests:**
1. ✅ `test_sra_data_structure()` 
   - Validates return types (list, dict, list)
   - Confirms experiments have proper structure
   - Confirms SRR codes (not SRX) in runs
   - Confirms biosamples have SAMN prefix
   - Confirms dict structure has `experiment_titles` and `samples` keys

2. ✅ `test_process_bioproject()`
   - Tests that process_bioproject handles new data structures
   - Validates all required fields present
   - Confirms counts are integers
   - Tests end-to-end workflow

**Test Results:**
```
✓ experiments is list: 12 items
✓ biosamples_dict is dict: 12 items
✓ sra_runs is list: 12 items (12 unique)
✓ Runs are SRR codes (not SRX)
✓ BioSample IDs are SAMN codes
✓ All runs are SRR codes
✓ process_bioproject returns correct fields
✓ All counts are integers

✓✓✓ ALL TESTS PASSED ✓✓✓
Hierarchical refactoring is working correctly!
```

---

## Before & After Example

### Input: BioProject PRJNA1381306

#### Before (Incorrect Output):
```csv
bioproject,sra_experiments_count,biosamples_count,sra_accessions_count,sra_experiments,biosamples,sra_accessions
PRJNA1381306,12,12,12,"SRX31557547; SRX31557546; ...","SAMN54118015; SAMN54118016; ...","SRX31557547; SRX31557546; ..."
```
❌ `sra_accessions` contains SRX codes (experiments), not actual run data

#### After (Correct Output):
```csv
bioproject,sra_experiments_count,biosamples_count,sra_runs_count,sra_experiments,biosamples,sra_runs
PRJNA1381306,12,12,12,"SRX31557547; SRX31557546; ...","SAMN54118015; SAMN54118016; ...","SRR36541090; SRR36541091; ..."
```
✅ `sra_runs` contains SRR codes (actual sequencing runs with real sequence data)

---

## Data Flow

```
BioProject ID (PRJNA*)
    ↓
fetch_sra_for_bioproject()
    ├─→ Calls SRAFetcher.fetch_all_by_bioproject()
    ├─→ extract_sra_experiment_metadata() returns:
    │  {
    │    'exp_accession': 'SRX31557547',
    │    'title': 'Experiment title',
    │    'biosample': 'SAMN54118015',
    │    'runs': ['SRR36541090', 'SRR36541091']  ← SRR codes!
    │    ...
    │  }
    ├─→ Aggregates into:
    │  - experiments: [exp1, exp2, ...]
    │  - biosamples_dict: {SAMN: {titles: [...]}}
    │  - sra_runs: [SRR, SRR, SRR, ...]
    ↓
process_bioproject()
    ├─→ Counts: experiments, biosamples, runs (all unique)
    ├─→ Stores lists: sra_experiments (SRX), biosamples (SAMN), sra_runs (SRR)
    ├─→ Calls search_pubmed_publications()
    │   ├─→ Extracts SAMN IDs from biosamples_dict
    │   ├─→ Uses sra_runs directly
    │   └─→ Searches PubMed cascading from SAMN → SRR
    ├─→ Enriches with publication metadata
    ↓
save_results_csv()
    └─→ Outputs CSV with rows containing:
        - sra_experiments_count, biosamples_count, sra_runs_count
        - sra_experiments, biosamples, sra_runs (all semicolon-separated)
```

---

## Key Improvements

1. **Correct Hierarchy**: Now properly represents BioProject → BioSample → Experiment → Run structure
2. **Accurate Run Codes**: Uses SRR codes (actual sequencing runs) not SRX codes (experiments)
3. **Better Data Organization**: Biosamples and experiments properly associated
4. **Consistent Naming**: Field names now accurately describe their content
5. **Experiment Titles**: Associated with their biosamples for traceability
6. **Tested & Validated**: Comprehensive test suite confirms correct functionality

---

## Backward Compatibility

⚠️ **Breaking Changes:**
- CSV output field names changed: `sra_accessions` → `sra_runs`
- CSV output field names changed: `sra_accessions_count` → `sra_runs_count`
- Data types changed: sets → dict/list
- Any existing CSV files won't have the new field names

✅ **Action Needed:**
- Update any dependent tools/scripts that read old CSV format
- Use new field names for analysis

---

## Next Steps (Optional Enhancements)

Future improvements that could be added:
1. Extract BioSample titles/names from NCBI (currently only storing IDs)
2. Add run statistics (e.g., read count, base pairs) to output
3. Include library strategy/source in output for better filtering
4. Add hierarchical display option in CSV output
5. Create nested JSON output option for better structure preservation

---

## Summary

✅ **Hierarchical structure refactored**
✅ **Correct biological codes (SRR, SAMN, SRX)**
✅ **All dependent methods updated**
✅ **Comprehensive tests passed**
✅ **Documentation updated**
✅ **Ready for production use**

The system now correctly models and stores the SRA hierarchy, enabling better analysis and linking to publications based on accurate biological sample and run information.
================================================================================
   HIERARCHICAL SRA DATA STRUCTURE REFACTORING - FINAL REPORT
================================================================================

COMPLETION STATUS: ✅ 100% COMPLETE - ALL TESTS PASSED

================================================================================
WHAT WAS FIXED
================================================================================

1. DATA STRUCTURE HIERARCHY
   OLD (Incorrect):
   - sra_accessions was actually SRX* codes (experiments), not run data
   - Didn't properly model: BioProject → BioSample → Experiment → Run
   
   NEW (Correct):
   - sra_runs now contains SRR* codes (actual sequencing data)
   - biosamples_dict properly structures SAMN (sample) codes
   - Properly models: BioProject → BioSample → Experiment (SRX) → Run (SRR)

2. CODE IDENTIFIERS CLARIFIED
   SRX* = Experiment/Dataset ID (e.g., SRX31557547)
   SAMN* = BioSample ID (e.g., SAMN54118015)
   SRR* = Run ID (actual sequence data) (e.g., SRR36541090)

================================================================================
FILES MODIFIED
================================================================================

1. Fetcher_NCBI/boolean_fetcher_integrated.py
   - fetch_sra_for_bioproject()      ← Changed return structure
   - search_pubmed_publications()    ← Updated to handle new types
   - process_bioproject()             ← Updated unpacking & field handling
   - save_results_csv()               ← Updated field names

2. README.md (Lines 59-77)
   - Updated CSV header to show sra_runs instead of sra_accessions
   - Added hierarchy explanation section
   - Updated examples with SRR codes instead of SRX codes

3. Created: REFACTORING_SUMMARY.md
   - Comprehensive documentation of all changes
   - Before/after comparison
   - Testing results

================================================================================
TEST RESULTS
================================================================================

✓ fetch_sra_for_bioproject() returns correct types:
  - experiments: list (12 items)
  - biosamples_dict: dict (12 items with SAMN* keys)
  - sra_runs: list (12 SRR* items, all unique)

✓ SRA codes are correct:
  - Experiments start with SRX: SRX31557547 ✓
  - BioSamples start with SAMN: SAMN54118015 ✓
  - Runs start with SRR: SRR36541090 ✓

✓ process_bioproject() works with new structure:
  - Returns all required fields
  - Correctly counts experiments, biosamples, runs
  - Passes correct data types to search_pubmed_publications()

✓ CSV output uses correct field names:
  - sra_runs_count (not sra_accessions_count)
  - sra_runs with SRR codes (not SRX codes)

================================================================================
BEFORE & AFTER DATA COMPARISON
================================================================================

BEFORE (Incorrect):
  CSV contains: sra_accessions = "SRX31557547; SRX31557546; ..."
  Problem: These are EXPERIMENT codes, not actual sequence data

AFTER (Correct):
  CSV contains: sra_runs = "SRR36541090; SRR36541091; ..."
  Correct: These are RUN codes with actual sequence data

HIERARCHY NOW PRESERVED:
  BioProject (PRJNA1381306)
    └─ BioSample (SAMN54118015)
        └─ Experiment (SRX31557547)
            └─ Run (SRR36541090)  ← Actual sequence data

================================================================================
VALIDATION CHECKLIST
================================================================================

✅ Data structure types correct (list, dict, list)
✅ All SRA codes have correct prefixes (SRX, SAMN, SRR)
✅ Hierarchy properly preserved and retrievable
✅ CSV output uses correct fields
✅ Documentation updated
✅ All dependent methods updated
✅ Comprehensive test suite passed
✅ End-to-end workflow tested
✅ No Python syntax errors
✅ All counts are correct integers

================================================================================
NEXT TIME THE USER RUNS THE SYSTEM
================================================================================

Option 2 (BioProject search) will now:
  1. Fetch correct SRA runs (SRR*) not experiments
  2. Store experiments (SRX*) for reference
  3. Store biosamples (SAMN*) for linking
  4. Output CSV with accurate hierarchical structure
  5. Use SRR codes for PubMed cascade search

CSV output will now correctly show:
  - sra_experiments: List of SRX* codes
  - biosamples: List of SAMN* codes
  - sra_runs: List of SRR* codes (biological sequence data)

================================================================================
SUMMARY
================================================================================

The hierarchical structure is now correctly implemented and tested.
The system properly models: BioProject → BioSample → Experiment → Run

The distinction between SRX (experiments) and SRR (runs) is now clear,
and the CSV output will contain the actual sequence run data (SRR codes)
instead of experiment metadata codes (SRX codes).

System is ready for production use.

Generated: 2026-02-12
Status: ✅ COMPLETE
# 📋 Titles & Metadata Integration - Implementation Summary

## ✅ What's Now Stored

Your system now captures and stores **complete hierarchical SRA structure** with:

### 1. **Experiment Titles** (Previously Missing)
```
✓ "Solanum lycopersicum leaf FsK-treated Drought-stress replicate 1"
✓ "Solanum lycopersicum leaf FsK-treated no Drought-stress replicate 3"
```

### 2. **Complete Experiment Metadata**
```json
{
  "library_name": "lib_fsk_drought_rp1",
  "library_strategy": "RNA-Seq",
  "library_source": "TRANSCRIPTOMIC",
  "library_selection": "PolyA",
  "library_layout": "SINGLE",
  "instrument": "Illumina NovaSeq 6000"
}
```

### 3. **Sample Attributes** (From NCBI XML)
```json
{
  "isolate": "esculentum",
  "cultivar": "Moneymaker",
  "age": "29 days",
  "dev_stage": "Vegetative stage [PO:0009021]",
  "collection_date": "2021-06",
  "geographic_location": "Greece:Thessaly,Larissa",
  "tissue": "leaf",
  "treatment": "FsK-treated",
  "stress": "yes",
  "replicate": "1"
}
```

### 4. **Clear Hierarchical Relationships**
```
Which experiment is associated with which sample?
Which runs belong to which experiment?
→ All clearly documented in sra_hierarchy
```

---

## 📊 JSON Output Structure

Your JSON now contains a new field: **`sra_hierarchy`**

### Structure:
```json
{
  "bioproject": "PRJNA1381306",
  "sra_experiments_count": 12,
  "biosamples_count": 12,
  "sra_runs_count": 12,
  "sra_hierarchy": {
    "SAMN54118015": {
      "sample_id": "SAMN54118015",
      "experiments": [
        {
          "experiment_id": "SRX31557547",
          "title": "Solanum lycopersicum leaf FsK-treated...",
          "metadata": {
            "library_name": "lib_fsk_drought_rp1",
            "library_strategy": "RNA-Seq",
            "library_source": "TRANSCRIPTOMIC",
            "library_selection": "PolyA",
            "library_layout": "SINGLE",
            "instrument": "Illumina NovaSeq 6000"
          },
          "sample_attributes": {
            "isolate": "esculentum",
            "cultivar": "Moneymaker",
            "age": "29 days",
            "dev_stage": "Vegetative stage [PO:0009021]",
            "collection_date": "2021-06",
            "tissue": "leaf",
            "treatment": "FsK-treated"
          },
          "runs": ["SRR36541090"]
        }
      ]
    },
    "SAMN54118014": { ... },
    ...
  }
}
```

---

## 🔄 Data Flow Changes

### Before:
```
Experiment XML → Extract basic info → Store in list
                                   ✗ No titles
                                   ✗ No sample details
                                   ✗ No hierarchy
```

### After:
```
Experiment XML → extract_sra_experiment_metadata()
                 ├─ Extract exp_title ✓
                 ├─ Extract instrument ✓
                 ├─ Extract library_name ✓
                 ├─ Extract sample_attributes ✓
                 └─ Now includes all metadata
                 
fetch_sra_for_bioproject() → Build structure:
            {
              "experiments": [...],
              "biosamples_dict": {
                "SAMN*": {
                  "experiments": [
                    {
                      "exp_accession": "SRX*",
                      "title": "...",  ✓ NEW
                      "metadata": {...},  ✓ NEW
                      "sample_attributes": {...},  ✓ NEW
                      "runs": ["SRR*"]
                    }
                  ]
                }
              }
            }

build_hierarchical_sra_structure() → Organize into hierarchy
                                    → Store in sra_hierarchy field
                                    → Include in JSON output
```

---

## 🎯 Files Modified

### 1. `Fetcher_NCBI/ncbi_fetcher_sra_fixed.py`
**Changes in `extract_sra_experiment_metadata()`:**
- ✅ Added `exp_title` extraction (was called `title`)
- ✅ Added `instrument` extraction from XML
- ✅ Added `sample_attributes` extraction
- Returns hierarchical dict with all metadata

**New fields extracted:**
- instrument
- sample_attributes (dict of all SAMPLE_ATTRIBUTE tags)

### 2. `Fetcher_NCBI/boolean_fetcher_integrated.py`

**Updated `fetch_sra_for_bioproject()`:**
- ✅ Now builds complete experiment info dicts
- ✅ Includes titles, metadata, attributes
- ✅ Returns structured biosamples_dict with full metadata

**New method `build_hierarchical_sra_structure()`:**
- ✅ Creates clean hierarchy: BioSample → Experiment → metadata/runs
- ✅ Organizes all information for easy access
- ✅ Returns well-structured dict

**Updated `process_bioproject()`:**
- ✅ Calls `build_hierarchical_sra_structure()`
- ✅ Stores result in `sra_hierarchy` field
- ✅ Logs hierarchy summary

**Updated `save_results_json()`:**
- ✅ Now includes sra_hierarchy in JSON output
- ✅ Preserves all metadata when saving

### 3. Created: `test_titles_metadata.py`
Comprehensive test that validates:
- ✅ Experiment titles are extracted
- ✅ All metadata fields are present
- ✅ Sample attributes are stored
- ✅ Hierarchy structure is correct
- ✅ JSON output includes everything

---

## 📁 What You'll See When You Run It

When you search for a BioProject, the JSON output will now show:

```
PRJNA1381306
├─ 12 Experiments (with titles!)
├─ 12 BioSamples (with attributes!)
└─ 12 Runs (with clear associations!)

In sra_hierarchy:
├─ SAMN54118015
│  ├─ SRX31557547: "Solanum lycopersicum leaf FsK-treated..."
│  │  ├─ Library: lib_fsk_drought_rp1
│  │  ├─ Instrument: Illumina NovaSeq 6000
│  │  ├─ Isolate: esculentum
│  │  ├─ Treatment: FsK-treated
│  │  └─ Run: SRR36541090
│  └─ (more experiments for same sample...)
└─ SAMN54118014
   ├─ SRX31557546: "Solanum lycopersicum leaf FsK-treated no Drought..."
   │  ├─ Library: lib_fsk_nodrought_rp3
   │  └─ ...
```

---

## 🧪 Validation

✅ **All tests passed:**
- Experiment titles correctly extracted
- Metadata fields populated
- Sample attributes stored
- Hierarchical structure correct
- JSON output includes everything
- No missing or truncated data

---

## 📝 Next Steps (Optional)

The data is now complete! If you want to further improve it:

1. **Export to CSV with hierarchy** - Create a better CSV format that shows the hierarchy
2. **Add experiment descriptions** - Extract longer description fields
3. **Add run statistics** - Include read count, base pairs, etc.
4. **Create a viewing tool** - Display the hierarchy in a readable format
5. **Filter by metadata** - Search/filter by treatment, tissue, etc.

---

## Summary

Your system now stores **complete hierarchical SRA data with all titles and metadata**. The JSON output (`sra_hierarchy` field) contains:

- ✅ Experiment titles
- ✅ Library names and technical details
- ✅ Sequencing instrument information
- ✅ Sample attributes (isolate, cultivar, age, tissue, treatment, etc.)
- ✅ Clear associations between samples, experiments, and runs
- ✅ All organized in a clean hierarchical structure

Everything is ready for analysis, visualization, or export to other formats!

---

**Test it:**
```bash
# Quick test:
python3 test_titles_metadata.py

# Then check the JSON:
cat /tmp/test_hierarchy.json | python3 -m json.tool
```

**Status**: ✅ **COMPLETE** - Ready to use!
