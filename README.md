# PYNER - Sistema Integrado de Búsqueda Científica

Sistema completo para buscar y vincular información de proyectos genómicos (BioProject), datos experimentales (SRA) y publicaciones científicas (PubMed) usando lenguaje natural.

---

## � Instalación

### Opción 1: Script Automático (Recomendado)

```bash
cd Pyner_PGRLAB
bash install_dependencies.sh
```

El script instalará automáticamente BioPython y verificará que todo esté listo.

### Opción 2: Manual

```bash
pip install biopython
```

**Nota:** El script principal valida automáticamente todas las dependencias antes de ejecutarse.

📖 **Más detalles:** Ver [INSTALLATION.md](INSTALLATION.md) para guía completa de instalación y troubleshooting.

---

## 🚀 Inicio Rápido

### Script Principal: PYNER Miner

```bash
bash pyner_miner.sh
```

**Flujo completo integrado:**
```
[1] Select Mode:
    → Lite: Fast fetch only (basic CSV/JSON)
    → Pro:  Full AI analysis (53 metadata columns + relevance scoring)

[2] Enter natural language query:
    > "tomato drought stress RNA-Seq"

[3] Choose database:
    → [1] PubMed - Fast literature search
    → [2] BioProject - Full omics data with cascade linking

[4] AI generates boolean query with synonyms

[5] Configure max results (unlimited by default)

[6] Execution:
    → Fetch: Basic publication metadata
    → Pro mode: AI analysis with Ollama (organisms, tissues, conditions, etc.)
    → Export: CSV + JSON in output/ directory
```

**Output:**
- **Lite mode:** Basic CSV/JSON with PMID, title, abstract, DOI
- **Pro mode:** Comprehensive CSV with 53 experimental metadata columns + AI-extracted details

---

## 📚 Documentación

### 📖 **[docs/FETCHER_DOCUMENTATION.md](docs/FETCHER_DOCUMENTATION.md)** ← Documentación completa

Guía completa del sistema con ejemplos, arquitectura y troubleshooting.

### 🔧 **[INSTALLATION.md](INSTALLATION.md)** ← Guía de instalación

Instrucciones detalladas de instalación y resolución de problemas.

---

## 🛡️ Robustez y Validaciones

El sistema incluye validaciones automáticas para garantizar portabilidad entre entornos:

✅ **Validación de Python 3** - Verifica que `python3` esté disponible  
✅ **Validación de BioPython** - Verifica instalación antes de ejecutar  
✅ **Validación de archivos requeridos** - Verifica que todos los módulos existan  
✅ **Manejo de dependencias opcionales** - FastAPI/Pydantic solo para servidor API  
✅ **Mensajes de error claros** - Instrucciones específicas para resolver problemas  
✅ **Uso de `python3` explícito** - Evita problemas con alias de Python 2  

**Ejemplo de validación:**
```bash
$ bash test_fetcher_integrator.sh
ERROR: python3 not found
  Please install Python 3: sudo apt install python3

ERROR: BioPython not installed
  Install with: pip install biopython
```

---

## ✨ Características

### Lite Mode (Fast)
✅ **Lenguaje Natural → Boolean Query** (IA con Ollama)  
✅ **Búsqueda en PubMed** para revisión bibliográfica rápida  
✅ **Búsqueda en BioProject** con query booleano + cascade linking (SRA)  
✅ **Export básico** CSV/JSON con PMID, título, abstract, DOI  

### Pro Mode (Comprehensive)
✅ Todo lo de Lite mode +  
✅ **Análisis IA de papers** con Ollama (qwen2.5:14b)  
✅ **53 columnas de metadata** experimental extraídas automáticamente:
  - Organismos, especies, genotipos, tejidos, células
  - Condiciones ambientales, temperatura, luz, medios de cultivo
  - Moléculas extraídas (RNA/DNA/Protein con tipos específicos)
  - Diseño temporal, replicación, grupos de tratamiento
  - Métricas de calidad, umbrales estadísticos, normalización
✅ **Scoring de relevancia** (0-10) basado en tu consulta  
✅ **Token usage monitoring** en tiempo real (prompts, respuestas, velocidad)  
✅ **Full-text PMC integration** cuando disponible  
✅ **Tabla clasificada completa** lista para análisis downstream  

---

## 🏗️ Estructura

```
Pyner_PGRLAB/
├── pyner_miner.sh                     # 🚀 Script principal (Lite/Pro mode)
├── output/                            # 📁 All results stored here
├── Query_generator/phases/phase3/    # 🤖 IA: Natural → Boolean query
├── Fetcher_NCBI/                      # 🔍 Search and linking
│   ├── boolean_fetcher_integrated.py # BioProject cascade workflow
│   ├── pubmed_boolean_search.py      # PubMed direct search
│   └── ncbi_fetcher_sra_fixed.py     # SRA fetcher
├── Data_Analyzer/                     # 📊 AI analysis with Ollama (Pro mode)
│   ├── paper_analyzer.py              # 53-column extraction + scoring
│   ├── ollama_client.py               # LLM interface (qwen2.5:14b)
│   └── pmc_fetcher.py                 # Full-text fetcher from PMC
└── archive_old/                       # 🗄️ Deprecated scripts & data
```

**Integrated pipeline:**
```
1. Natural Language Input  →  Query_generator (AI with Ollama)
                               ↓
2. Boolean NCBI Query      →  Fetcher_NCBI (PubMed/BioProject)
                               ↓
3. Publications Retrieved  →  Data_Analyzer (Pro mode only)
                               ↓
4. Final Output            →  output/ directory (CSV + JSON)
```

**Modes:**
- **Lite:** Steps 1-2 only (fast fetch)
- **Pro:** Steps 1-3 (full AI analysis with 53 metadata columns)

---

## 📊 Output CSV

### Opción 1: BioProject

```csv
bioproject,title,organism,sra_experiments_count,biosamples_count,sra_runs_count,
sra_experiments,biosamples,sra_runs,publications_found,search_method,dois,pmids
```

**Ejemplo:**
```csv
PRJNA123456,"Study Title",Arabidopsis,42,12,156,"SRX123; SRX124; ...","SAMN123; SAMN124; ...","SRR1234567; SRR1234568; ...",2,direct,"10.1234/abc","35123456"
```

**Notas sobre la jerarquía:**
- **sra_experiments**: Lista de experimentos (SRX*) del BioProject
- **biosamples**: Lista de muestras biológicas (SAMN*) asociadas a los experimentos
- **sra_runs**: Lista de corridas de secuenciación (SRR*) que contienen los datos reales
- Estructura: `BioProject → BioSample → Experiment (SRX) → Run (SRR)`
- Todas las listas de IDs están separadas por punto y coma (`;`) para análisis posterior

### Opción 2: PubMed Direct

```csv
pmid,title,year,journal,publication_type,authors,doi,pmcid,url,abstract,fetched_at
```

**Ejemplo:**
```csv
35099557,"The tomato OST1-VOZ1...",2022,"The Plant cell","Journal Article; Research Support, Non-U.S. Gov't","Chong L; Xu R; Huang P",10.1093/plcell/koac026,PMC9048945,https://...,Abstract...,2026-02-12
```

---

## 🔧 Configuración

1. **Instalar dependencias:**
   ```bash
   cd Fetcher_NCBI
   pip install -r requirements.txt
   ```

2. **API Key NCBI** (opcional pero recomendado):
   - Obtener en: https://www.ncbi.nlm.nih.gov/account/
   - Editar `Fetcher_NCBI/config.py`:
     ```python
     NCBI_API_KEY = "your_api_key_here"
     ```

---

## 📖 Uso

### Modo Interactivo (Recomendado)
```bash
bash pyner_miner.sh
```

Interactive menu guides you through:
1. Mode selection (Lite/Pro)
2. Natural language query input
3. Database selection (PubMed/BioProject)
4. Query generation and confirmation
5. Max results configuration
6. Automated execution

### Modo Directo (Advanced)

**Fetch only (Lite equivalent):**

PubMed:
```bash
cd Fetcher_NCBI
python pubmed_boolean_search.py "Arabidopsis phosphate stress" --max 50 \
  --output-csv ../output/results.csv --output-json ../output/results.json
```

BioProject:
```bash
cd Fetcher_NCBI
python boolean_fetcher_integrated.py "Tomato drought RNA-Seq" --max 20 \
  --output-csv ../output/results.csv --output-json ../output/results.json
```

**Full analysis (Pro equivalent):**
```bash
# Step 1: Fetch
cd Fetcher_NCBI
python pubmed_boolean_search.py "query" --max 50 \
  --output-json ../output/fetch.json

# Step 2: Analyze
cd ../Data_Analyzer
python paper_analyzer.py ../output/fetch.json ../output/classified.csv
```

---

## ⚡ Performance

| Operación | 10 BioProjects | 50 BioProjects |
|-----------|----------------|----------------|
| Total | ~10 min | ~45 min |

**Recomendaciones:** Testing `--max 5`, Producción `--max 50`

---

## 📚 Más Información

Lee **[GUIA_COMPLETA.md](GUIA_COMPLETA.md)** para:
- Arquitectura detallada
- Cascade PubMed search
- Troubleshooting
- Ejemplos avanzados

---

**Versión:** 2.0.0 | **Fecha:** 2026-02-13 | **Estado:** ✅ Producción  
**Script principal:** `pyner_miner.sh` | **Modos:** Lite (fast) / Pro (comprehensive)
