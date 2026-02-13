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

```bash
bash test_fetcher_integrator.sh
```

**Ejemplo de uso:**
```
> Arabidopsis under phosphate stress
→ Choose database: [1] PubMed  [2] BioProject
→ Generates boolean query with synonyms
→ Searches selected database
→ Extracts data (SRA if BioProject, metadata if PubMed)
→ Exports CSV + JSON
```

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

✅ **Lenguaje Natural → Boolean Query** (IA con Ollama)  
✅ **Búsqueda en BioProject** con query booleano + extracción SRA  
✅ **Búsqueda directa en PubMed** para revisión bibliográfica rápida  
✅ **Cascade PubMed Linking** (3 niveles de búsqueda para BioProject)  
✅ **Marca "NA"** cuando no hay publicaciones  
✅ **Export CSV + JSON** con metadata completa  

---

## 🏗️ Estructura

```
Pyner_PGRLAB/
├── GUIA_COMPLETA.md                    # 📖 Documentación principal
├── test_fetcher_integrator.sh         # 🚀 Script principal (selección de DB)
├── Query_generator/phases/phase3/     # 🤖 IA: Natural → Boolean
└── Fetcher_NCBI/                       # 🔍 Búsqueda y linking
    ├── boolean_fetcher_integrated.py  # BioProject workflow
    ├── pubmed_boolean_search.py       # PubMed direct search
    ├── ncbi_fetcher_sra_fixed.py      # SRA fetcher
    └── ncbi_linkout.py                # PubMed linking
```

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

### Modo Interactivo
```bash
bash test_fetcher_integrator.sh
```

### Modo Directo

**BioProject:**
```bash
cd Fetcher_NCBI
python boolean_fetcher_integrated.py "Arabidopsis phosphate" --max 20 --output-csv results.csv
```

**PubMed:**
```bash
cd Fetcher_NCBI
python pubmed_boolean_search.py "Arabidopsis phosphate stress" --max 50 --output-csv pubmed.csv
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

**Versión:** 1.0.0 | **Fecha:** 2026-02-12 | **Estado:** ✅ Producción
