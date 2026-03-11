# Synonym Builder — Contexto Completo
**Última actualización:** 2026-03-10
**Estado:** ✅ **COMPLETADO Y LISTO PARA PRODUCCIÓN**

---

## 🎯 Objetivo del módulo

Construir un sistema de expansión de sinónimos para el **Query Generator** de PYNER que:
1. Pre-computa sinónimos de términos MeSH + variantes ortográficas reales en un diccionario unificado (O(1)).
2. Implementa **Heurísticas Protectoras Definitivas (Zero-Bias)**:
   - *Shield Genético*: Evita reinterpretaciones de genes mediante validación en el Diccionario Universal de Genes.
   - *Shield de Elementos/Nutrientes*: Impide que términos base químicos colapsen hacia vías/enzimas.
   - *Rastreo de Condición Experimental*: Mediante LLM (`Qwen2.5`) y un diccionario masivo de corpus n-grams, extrae y enriquece estreses y tratamientos.
   - *Detección de Frases Compuestas*: Aísla y respeta el anclaje de un organismo + un proceso ("wheat development").
3. Genera queries NCBI booleanas altamente precisas sin sesgos médicos inducidos (PubMedBERT fue eliminado totalmente).

**Ejemplo:**
```
Input: "Arabidopsis WRKY33 heat stress drought tolerance"

Output NCBI Query:
("A. thaliana" OR ... OR "Arabidopsis") AND
("ATWRKY33" OR ... OR "WRKY33") AND
("heat stress" OR "heat treatment" OR ...) AND
("drought tolerance" OR "drought stress" OR ...)
```

---

## 📋 Pipeline Completado

### ✅ Step 1 — MeSH Official Synonyms
**Script:** `step1_mesh_synonyms.py`
**Output:** `output/step1_mesh_synonyms.json` (13 MB)

**Contenido:**
- 28,157 términos MeSH con sinónimos oficiales curados por NLM
- Formato: `{term_lowercase: {preferred_term, synonyms[], mesh_ui}}`

**Fuente:** `pubmed_explorer/mesh_data/desc2026.xml`

---

### ✅ Step 2A — PubMed Corpus Builder
**Script:** `batch_download_parallel.py`
**Output:** `output/step2_mesh_consolidated_parallel_1000s.json` (9.5 GB)

**Contenido:**
- Por cada término MeSH: hasta 1000 abstracts de PubMed Baseline (26M+ artículos)
- Incluye: Título + Abstract + Keywords

**Configuración:**
- 50 workers paralelos
- 1000 samples por término MeSH

---

### ✅ Step 2B — Text Variant Finder
**Script:** `step2b_find_variants_1000s.py`
**Output:** `output/step2_text_variants_1000s.json` (19 MB)

**Contenido:**
- Variantes ortográficas reales observadas en papers científicos
- Con frecuencia y confidence score

**Algoritmo:**
- Fuzzy matching (Levenshtein, threshold 0.70)
- N-gramas (hasta 4 palabras)
- Regex para variantes de organismos

**Nota:** `step2b_find_variants.py` contiene la lógica core (NO borrar)

---

### ✅ Step 3 — Dictionary Merge (MeSH + Text Variants)
**Script:** `step3_merge_dictionary.py`
**Output:** `output/final_synonym_dictionary.json` (~40 MB)

**Contenido:**
- Diccionario unificado para lookup O(1)
- ~28,157 términos MeSH
- ~250,000+ sinónimos totales:
  - ~205,000 sinónimos MeSH oficiales
  - ~45,000+ variantes ortográficas

**Estructura:**
```json
{
  "rna-seq": {
    "preferred_term": "RNA-Seq",
    "mesh_ui": "D000073888",
    "synonyms": [
      {"term": "Whole Transcriptome Shotgun Sequencing", "source": "mesh", "confidence": "official"},
      {"term": "rnaseq", "source": "text", "confidence": 0.95, "frequency": 523},
      {"term": "rna sequencing", "source": "text", "confidence": 0.85, "frequency": 187}
    ],
    "synonym_count": {"total": 15, "mesh": 1, "text": 14}
  }
}
```

---

### ✅ Query Expander con Multi-Nivel de Análisis Heurístico
**Script:** `query_expander.py`
**Dependencias:** `final_synonym_dictionary.json`, `final_gene_dictionary.json`, `condition_lexicon.json` y Ollama LLM.

**Estrategia híbrida:**

1. **Búsqueda Exacta y Léxica:**
   ```python
   Usuario: "drought"
   → Busca en diccionario MeSH pre-computado
   → Devuelve: ["Droughts", "water deficit", "Desiccation", ...]
   ```

2. **Detección Dinámica de Escudos:**
   - **Genes:** "WRKY33" carga "Gene Shield", cruza genes y prohíbe traducciones abstractas.
   - **Elementos Químicos:** "Zinc" se marca como nutriente, solo agrega "ion" y omite MeSH erróneos (como "Zinc Fingers").
   - **Compuestos Organismo+Proceso:** "Wheat development", se expande el taxon pero el proceso se integra y blinda.
   - **Condiciones Experimentales:** "heat stress" interactúa con el `Qwen2.5` y el _Lexicon de n-grams_ para obtener "heat exposed", "heat treatment", etc.; bloqueando asociaciones médicas directas de MeSH para botánica/microbiología.

**Uso:**
```bash
# Testing interactivo del Expansor Completo (incluye extracciones LLM)
python3 test_query_expander.py "Wheat heat stress grain quality metabolomics"
```

---

## 📁 Estructura de archivos (FINAL)

```
synonym_builder/
├── CONTEXT.md                          ← Este archivo
├── config.py                           ← Configuración global
│
├── Scripts ejecutados:
│   ├── step1_mesh_synonyms.py          ✅ Extrae sinónimos MeSH
│   ├── batch_download_parallel.py      ✅ Descarga abstracts PubMed
│   ├── step2b_find_variants.py         ✅ Variantes de ortografía
│   ├── step3_merge_dictionary.py       ✅ Merge MeSH + Text Variants
│   ├── step4_taxonomy_builder.py       ✅ Extrae el NCBI TaxDump (~1GB)
│   ├── step5_lexicon_builder.py        ✅ Lexicon n-gram de las bases textuales (~2M)
│   └── step6_gene_dictionary.py        ✅ Construcción comprimida NCBI Gene + UniProt
│
├── Módulo de producción:
│   ├── query_expander.py               🚀 Query expander robusto sin PubMedBERT
│   └── test_query_expander.py          🧪 CLI interactivo End-to-End
│
├── gene_data/                          ← Bases NCBI/UniProt en crudo
├── taxdump/                            ← NCBI Taxonomy dump extraído
├── output/
│   ├── step1_mesh_synonyms.json        (13 MB)
│   ├── step2_mesh_consolidated...json  (9.5 GB corpus)
│   ├── condition_lexicon.json          (22 MB - corpus experimental)
│   ├── final_synonym_dictionary.json   (40 MB) ⭐ DICCIONARIO MAIN
│   ├── final_taxonomy_dictionary.json  (170 MB) ⭐ DICCIONARIO TAXONES
│   └── final_gene_dictionary.json      (~120 MB) ⭐ DICCIONARIO GENES
```

---

### **Integración con pyner_miner.sh**

El wrapper principal de Pyner llama internamente a `test_query_expander.py` si seleccionas opciones de Query Generation. Toma la frase original, usa a Ollama para extraer componentes, evalúa cada uno y entrega la cadena final de NCBI.

```bash
# Interactivo (con output extendido)
python3 test_query_expander.py "water stress tomato gene expression"

# Básico sin LLM
python3 test_query_expander.py --no-llm "sequía tomate"
```

---

## 📊 Estadísticas del Sistema

### Cobertura de los Diccionarios
- **Términos MeSH + Texto:** ~250,000+
- **NCBI Taxonomy:** ~1,200,000+ especies y ramas
- **NCBI Genes / UniProt:** ~3,000,000+ genes y alias pre-comprimidos
- **Experimental Corpus N-grams:** ~2,256,780 frases

### Performance
| Operación | Tiempo | RAM |
|-----------|--------|-----|
| Cargar diccionarios base | ~3.5s | 500 MB |
| Lazy-Load Genes (si se necesita) | ~8s | 2 GB |
| Query NLP -> Keywords LLM | ~0.8s | (Llamada Ollama local) |
| Lookup Multi-Fase/Boolean Query | < 1s | - |

### Precisión
- **Zero-Bias:** No se añaden palabras de campos que el usuario no especificó (ni enfermedades no relacionadas, ni funciones hiper-específicas de elementos químicos). Todo se basa en validaciones de bases maestras.

---

## 🔧 Configuración (config.py)

```python
BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / "output"

# Step 1: MeSH
MESH_DESC_FILE = ../pubmed_explorer/mesh_data/desc2026.xml
MESH_INCLUDE_PLURAL = True

# Step 2B: Text Variants
ABSTRACT_MIN_FREQUENCY = 5      # Mínimo 5 apariciones
ABSTRACT_MIN_PERCENTAGE = 0.05  # Mínimo 5% de abstracts
MAX_NGRAM_SIZE = 4              # N-gramas hasta 4 palabras

# Query Expander
SEMANTIC_THRESHOLD = 0.80       # Fallback mínimo 80% similaridad
MAX_FALLBACK_RESULTS = 3        # Máximo 3 sugerencias semánticas
```

---

## 🧩 Integración con PYNER

```
Pyner_PGRLAB/
├── Query_generator/phases/
│   ├── phase1/synonym_builder/   ← ESTE MÓDULO (completo)
│   │   ├── query_expander.py     ← Usar en producción
│   │   └── output/final_synonym_dictionary.json  ← Diccionario principal
│   │
│   ├── phase2/                   ⚠️ DEPRECADO (no usar, datos de SRA)
│   │
│   └── phase3/api/               ← INTEGRAR AQUÍ
│       ├── main.py               ← Entry point
│       └── query_generator.py    ← Modificar para usar query_expander.py
│
├── Fetcher_NCBI/                 ← Descarga papers con boolean query
└── Data_Analyzer/                ← Reanálisis con LLM
```

**Flujo completo actualizado:**
```
1. Input NL (usuario): "arabidopsis drought stress rnaseq"
2. LLM extrae keywords: ["arabidopsis", "drought stress", "rnaseq"]
3. query_expander.py:
   - "arabidopsis" → lookup exacto → "Arabidopsis thaliana" + sinónimos
   - "drought stress" → lookup exacto → "Droughts" + variantes
   - "rnaseq" → lookup exacto → "RNA-Seq" + variantes
4. Boolean query: ("Arabidopsis..." OR ...) AND ("Droughts" OR ...) AND ("RNA-Seq" OR ...)
5. Fetcher_NCBI → papers
6. Data_Analyzer → clasificación LLM
```

---

## ⚠️ Decisiones Importantes

### ❌ PubMedBERT Eliminado (Fallback Semántico DEPRECADO)
**Razón:** Para tópicos de botánica y biología molecular genera "sesgos médicos" inaceptables y reclasificaba términos puramente funcionales hacia enfermedades u órganos de mamíferos. Consumía RAM alta y frenaba el arranque.
**Reemplazo:** Heurísticas LLM guiadas con léxico experimental. Todo bajo el control del `Qwen2.5`.

### ❌ Phase 3 (Antiguo generador simple) REEMPLAZADA
**Razón:** Era muy básica y no comprendía organismos frente a genes. El nuevo test_query_expander.py es la herramienta de interacción en lenguaje natural para la integración final.

---

## 🧪 Testing

### Test del pipeline completo (Lenguaje Natural → Query)
```bash
python3 test_query_expander.py "drought heat stress WRKY33"
```

---

## 📚 Próximos pasos

1. ✅ **Diccionarios base superados (Genes, Taxonomía, N-grams)**
2. ✅ **Query expander implementado y sin sesgos de PubMedBERT**
3. ✅ **pyner_miner.sh actualizado para delegar la interfaz al Expander**
4. 🚀 **Deploy final del CLI**

---

**Estado:** ✅ **LISTO PARA PRODUCCIÓN**
**Maintainer:** PGR Lab
**Última actualización:** 2026-03-10

Para integración con Phase 3, ver ejemplos en `query_expander.py`
