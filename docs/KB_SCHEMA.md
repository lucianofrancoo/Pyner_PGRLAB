# 📊 Knowledge Base Schema - Pyner PGRLAB

Este documento describe la estructura y esquema del Knowledge Base generado en **FASE 1** (Extracción de Metadatos).

---

## 📁 Archivos de Salida

```
data/kb/
├── organisms_index.json           # Organismos únicos + sinónimos
├── strategies_index.json          # Estrategias de secuenciación
├── sources_index.json             # Fuentes de librería
├── selections_index.json          # Selecciones de librería
├── tissues_index.json             # Tipos de tejido
├── treatments_index.json          # Tratamientos/condiciones
└── global_index.json              # Índice invertido completo
```

---

## 🧬 Schema Detallado

### 1. organisms_index.json

```json
{
  "arabidopsis_thaliana": {
    "scientific_name": "Arabidopsis thaliana",
    "common_names": ["arabidopsis", "thale cress", "mouse-ear cress"],
    "taxon_id": 3702,
    "frequency": 12453,
    "related_organisms": ["arabidopsis_lyrata", "capsella_rubella"],
    "studies": {
      "total": 12453,
      "rna_seq": 8934,
      "wgs": 2103,
      "amplicon": 1416
    },
    "attributes": {
      "tissues": ["leaf", "root", "flower", "seed"],
      "treatments": ["drought", "cold", "nitrogen", "light"],
      "genotypes": ["col-0", "ws-0", "ler-0"]
    },
    "last_update": "2026-02-06",
    "data_sources": ["sample.xml"]
  },
  
  "solanum_lycopersicum": {
    "scientific_name": "Solanum lycopersicum",
    "common_names": ["tomato", "garden tomato"],
    "taxon_id": 4081,
    "frequency": 3421,
    ...
  },
  
  "homo_sapiens": {
    "scientific_name": "Homo sapiens",
    "common_names": ["human"],
    "taxon_id": 9606,
    ...
  }
}
```

**Notas:**
- `frequency`: Número total de estudios con este organismo
- `studies`: Desglose por estrategia de secuenciación
- `attributes`: Metadatos comunes encontrados en SRA para este organismo
- Ordenado por frecuencia (DESC)

---

### 2. strategies_index.json

```json
{
  "RNA-Seq": {
    "library_strategy": "RNA-Seq",
    "aliases": ["transcriptomics", "transcript-seq", "whole transcriptome"],
    "description": "High-throughput RNA sequencing for gene expression profiling",
    "frequency": 89234,
    "organisms": {
      "arabidopsis_thaliana": 8934,
      "homo_sapiens": 34521,
      "mus_musculus": 21456,
      ...
    },
    "typical_library_sources": ["TRANSCRIPTOMIC"],
    "typical_library_selections": ["cDNA", "PolyA"],
    "keywords": ["expression", "transcription", "transcript", "rna", "gene"],
    "related_protocols": [
      "Illumina TruSeq RNA-Seq",
      "10x Genomics Single Cell",
      "Pacific Biosciences Iso-Seq"
    ]
  },
  
  "WGS": {
    "library_strategy": "WGS",
    "aliases": ["whole genome sequencing", "genomic", "genome-seq"],
    "description": "Whole genome sequencing for structural variation and SNP calling",
    "frequency": 45123,
    ...
  },
  
  "AMPLICON": {
    "library_strategy": "AMPLICON",
    "aliases": ["16S", "18S", "amplicon sequencing", "targeted amplicon"],
    "description": "Targeted PCR amplicon sequencing",
    "frequency": 23456,
    ...
  }
}
```

**Notas:**
- `aliases`: Variantes de nombres del usuario para el mismo concepto
- `keywords`: Palabras que típicamente aparecen en descripciones de este tipo
- `typical_library_*`: Lo que usualmente acompaña a esta estrategia

---

### 3. tissues_index.json

```json
{
  "leaf": {
    "common_names": ["leaves", "foliage", "frond"],
    "frequency": 34521,
    "organisms": {
      "arabidopsis_thaliana": 3456,
      "oryza_sativa": 2876,
      "zea_mays": 2145,
      ...
    },
    "ncbi_synonyms": ["plant leaf", "leaf tissue", "foliar"],
    "ncbi_biosamples": [
      "SAMD00016353",  // Ejemplos reales
      "SAMD00016354",
      ...
    ],
    "ontology": {
      "plant_part": "PO:0025034",  // Plant Ontology
      "description": "The above-ground lateral appendage of a plant stem"
    }
  },
  
  "root": {
    "common_names": ["roots", "root system", "rhizosphere"],
    "frequency": 28901,
    ...
  },
  
  "blood": {
    "common_names": ["whole blood", "blood sample", "serum"],
    "frequency": 15643,
    "organisms": {
      "homo_sapiens": 14521,
      "mus_musculus": 987,
      ...
    }
  }
}
```

---

### 4. treatments_index.json

```json
{
  "drought": {
    "treatment_name": "drought",
    "synonyms": [
      "water stress",
      "water deficit",
      "dehydration",
      "desiccation",
      "drought stress",
      "water-limited"
    ],
    "frequency": 5234,
    "organisms": {
      "arabidopsis_thaliana": 2145,
      "oryza_sativa": 876,
      "solanum_lycopersicum": 654,
      ...
    },
    "duration_range": {
      "min_hours": 1,
      "max_hours": 8736,
      "common_durations": [24, 48, 72, 168]
    },
    "measurement_types": [
      "stomatal_conductance",
      "leaf_water_potential",
      "relative_water_content"
    ],
    "related_treatments": [
      "heat_stress",
      "osmotic_stress",
      "saline_stress"
    ],
    "ncbi_attributes": {
      "stress_type": "environmental",
      "tag_variations": [
        "drought",
        "treatment: drought",
        "condition: drought stress"
      ]
    }
  },
  
  "nitrogen": {
    "treatment_name": "nitrogen",
    "synonyms": [
      "nitrogen fertilizer",
      "N fertilization",
      "nitrogen supply",
      "nitrogen availability"
    ],
    "frequency": 3876,
    ...
  }
}
```

**Notas:**
- `duration_range`: Estadísticas sobre duración típica del tratamiento
- `measurement_types`: Qué se típicamente mide con este tratamiento
- `ncbi_attributes`: Cómo aparece en los tags de NCBI

---

### 5. global_index.json

Índice invertido para búsqueda rápida (término → dónde aparece):

```json
{
  "arabidopsis": {
    "type": "organism_alias",
    "resolves_to": "arabidopsis_thaliana",
    "frequency": 12453,
    "confidence": 0.99
  },
  
  "drought": {
    "type": "treatment",
    "resolves_to": "drought",
    "frequency": 5234,
    "context": {
      "organisms": ["arabidopsis_thaliana", "oryza_sativa"],
      "strategies": ["RNA-Seq", "WGS"]
    }
  },
  
  "rna-seq": {
    "type": "strategy_alias",
    "resolves_to": "RNA-Seq",
    "frequency": 89234
  },
  
  "tomato": {
    "type": "organism_common_name",
    "resolves_to": "solanum_lycopersicum",
    "frequency": 3421
  }
}
```

---

## 📐 Estructura de Directorios de Datos

```
data/
├── kb/
│   ├── organisms_index.json          # 🔴 REQUIERE: Fase 1.1
│   ├── strategies_index.json         # 🔴 REQUIERE: Fase 1.1
│   ├── sources_index.json            # 🔴 REQUIERE: Fase 1.1
│   ├── selections_index.json         # 🔴 REQUIERE: Fase 1.1
│   ├── tissues_index.json            # 🔴 REQUIERE: Fase 1.2
│   ├── treatments_index.json         # 🔴 REQUIERE: Fase 1.2
│   ├── global_index.json             # 🟢 GENERADO: Fase 1.3
│   └── metadata/
│       ├── kb_stats.json             # Estadísticas generales
│       ├── processing_log.txt        # Log de procesamiento
│       └── timestamp.txt             # Cuándo se procesó
│
├── results/
│   ├── baseline_queries.csv          # Queries antiguas (LLM)
│   └── optimized_queries.csv         # Queries nuevas (QB)
│
└── raw/
    └── sample_xmls/                  # Copias de XMLs para testing
        ├── DRA000001.experiment.xml
        ├── DRA000001.sample.xml
        └── ...
```

---

## 🔄 Generación del KB

### Fase 1.1: Parsing básico
```bash
python scripts/extract_metadata.py \
    --input /home/lahumada/disco1/NCBI_Metadata/SRA \
    --max-files 100000 \
    --output data/kb/
```

**Genera:**
- organisms_index.json
- strategies_index.json
- sources_index.json
- selections_index.json

### Fase 1.2: Extracción de atributos
```bash
python scripts/extract_attributes.py \
    --kb-path data/kb/ \
    --input /home/lahumada/disco1/NCBI_Metadata/SRA
```

**Genera:**
- tissues_index.json (actualiza organisms_index.json)
- treatments_index.json
- kb_stats.json

### Fase 1.3: Índice invertido global
```bash
python scripts/build_global_index.py \
    --kb-path data/kb/ \
    --output data/kb/global_index.json
```

**Genera:**
- global_index.json

---

## 📊 Estadísticas Típicas (Esperadas)

```json
{
  "total_studies_processed": 100000,
  "total_files_processed": 300000,
  "organisms_unique": 4532,
  "strategies_unique": 67,
  "tissues_unique": 1245,
  "treatments_unique": 3456,
  "parsing_errors": 234,
  "parsing_success_rate": 99.92,
  "disk_space_kb": 487632,
  "processing_time_hours": 2.5,
  "processing_rate_files_per_sec": 33.3,
  "timestamp": "2026-02-06T14:30:00Z",
  "git_commit": "a041f81"
}
```

---

## 🔍 Cómo Usar el KB en Pyner_search

### Antes (sin KB, con LLM genérico)
```python
# User input
"arabidopsis drought"

# LLM:
# "Generate NCBI query for arabidopsis and drought"
# Output (genérico, sin contexto):
# arabidopsis[Organism] AND (drought OR water stress)

# Resultado: MUCHOS falsos positivos, searches imprecisas
```

### Después (con KB)
```python
# User input
"arabidopsis drought"

# QueryOptimizer:
# 1. Normalizador: "arabidopsis" → "arabidopsis_thaliana"
# 2. Busca en KB:
#    - organisms["arabidopsis_thaliana"]["scientific_name"]
#    - treatments["drought"]["synonyms"]
# 3. Construye query optimizada:
query = "Arabidopsis thaliana[Organism] AND " \
        "(drought OR \"water stress\" OR \"water deficit\" " \
        "OR dehydration[All Fields])"

# Resultado: Precisión 3-5x mejor, menos ruido
```

---

## ✅ Validación del KB

### Checklist de calidad
- [ ] Todos los organismos tienen taxon_id válido (NCBI Taxonomy)
- [ ] Todas las estrategias existen en NCBI SRA (67 valores permitidos)
- [ ] Sinónimos son verificables en literatura bioinformática
- [ ] Frecuencias suman correctamente
- [ ] No há duplicados en organismos/tratamientos
- [ ] JSON es válido y parseable
- [ ] Tamaño total < 500 MB

---

## 🚀 Próximos Pasos

1. **Implementar extractores** (Fase 1.1)
2. **Procesar 100k archivos** (Fase 1.1)
3. **Validar calidad del KB** (Fase 1.2)
4. **Integrar con QueryOptimizer** (Fase 2)
5. **Evaluar mejoras** (Fase 3)

---

## 📚 Referencias

- [NCBI SRA XML Schema](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/document.cgi?study_id=phs000001.v3.p1&phd=3346)
- [NCBI Query Syntax](https://www.ncbi.nlm.nih.gov/books/NBK3827/)
- [Plant Ontology](http://purl.obolibrary.org/obo/PO)
- [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy/)

---

**Versión:** 1.0  
**Última actualización:** Febrero 6, 2026
