# 🛠️ Plan de Implementación: Phase 3b - NCBI API Integration

## Resumen Ejecutivo

**Objetivo**: Integrar búsqueda real en NCBI con enriquecimiento del KB

**Duración estimada**: 2-3 horas

**Resultado**: API que transforma preguntas naturales en búsquedas NCBI + análisis local

---

## 📋 Componentes a Crear

### 1. NCBI E-utilities Wrapper
**Archivo**: `phase3/api/ncbi_integration.py`

```python
class NCBISearcher:
    """
    Wrapper para NCBI E-utilities API
    
    Métodos:
    - search(query: str) → {count, webenv, query_key}
    - fetch(webenv, query_key, start, count) → [{id, title, organism, ...}]
    - parse_xml_response(xml) → [{study_metadata}]
    """
```

**Responsibilities**:
- Conectar a NCBI E-utilities (eutils.ncbi.nlm.nih.gov)
- Ejecutar búsquedas booleanas
- Parsear respuestas XML
- Manejar rate limiting
- Transformar resultados a JSON

### 2. Query Generator Mejorado
**Archivo**: `phase3/api/query_generator.py`

```python
class NCBIQueryGenerator:
    """
    Genera queries NCBI desde lenguaje natural
    
    Métodos:
    - generate_query(natural_language: str) → str
    - validate_syntax(query: str) → bool
    - add_filters(query: str, filters: dict) → str
    """
```

**Responsibilities**:
- Usar LLM para transformar --> NCBI syntax
- Validar variables booleanas
- Añadir filtros (organism, strategy, etc.)
- Optimizar para cobertura máxima

### 3. Caching Layer
**Archivo**: `phase3/cache/query_cache.db`

```python
class QueryCache:
    """
    SQLite cache para evitar re-queries
    
    Estructura:
    - query_hash(query) → int
    - store(query, results, timestamp)
    - retrieve(query) → cached_results
    - is_valid(query) → boolean (< 7 días)
    """
```

**Responsibilities**:
- Cachear queries frecuentes
- Timestamp de búsqueda
- Evitar rate limiting
- Respuestas rápidas

### 4. Result Enricher (KB Integration)
**Archivo**: `phase3/api/result_enricher.py`

```python
class ResultEnricher:
    """
    Enriquece resultados NCBI con KB metadata
    
    Métodos:
    - enrich_result(ncbi_result, kb) → {ncbi_data + kb_stats}
    - score_relevance(result, kb) → float
    - correlate_organisms(organism, kb) → [{related_organisms}]
    """
```

**Responsibilities**:
- Buscar en KB datos del organismo
- Añadir estadísticas
- Correlacionar con estudios similares
- Generar scores de relevancia

### 5. API Endpoint Mejorado
**Archivo**: `phase3/api/main.py` (modificar)

```python
@app.post("/search/real")
async def search_ncbi_with_enrichment(request: SearchRequest):
    """
    Nueva ruta mejorada:
    1. Generar query NCBI
    2. Ejecutar en NCBI
    3. Enriquecer con KB
    4. Retornar resultados
    """
```

---

## 🔄 Flujo de Ejecución

```
POST /search/real
│
├─ 1. Input Validation
│   └─ natural_language: str
│
├─ 2. Query Generation (LLM)
│   ├─ Input: "CRISPR en arabidopsis"
│   ├─ LLM: qwen2.5:14b
│   └─ Output: "((CRISPR) AND (arabidopsis OR A. thaliana))"
│
├─ 3. Cache Check
│   ├─ hash(query)
│   ├─ Si existe: return cached (< 10ms)
│   └─ Si no existe: continue
│
├─ 4. NCBI Search
│   ├─ POST to eutils API
│   ├─ esearch → get count + webenv
│   └─ efetch → get full results (XML)
│
├─ 5. Parse Results
│   ├─ XML → JSON
│   ├─ Extract: id, title, organism, strategy...
│   └─ Top 100 results
│
├─ 6. Enrich with KB
│   ├─ Para cada resultado:
│   │  ├─ Lookup organism en KB
│   │  ├─ Add statistics
│   │  ├─ Calculate relevance score
│   │  └─ Correlate with similar studies
│   └─ Return enriched results
│
├─ 7. Cache Results
│   ├─ Store: query + enriched_results + timestamp
│   └─ TTL: 7 días
│
└─ 8. Return Response (JSON)
    ├─ Query generated
    ├─ NCBI results count
    ├─ Top 10 enriched results
    ├─ Analysis summary
    ├─ Lookup time
    └─ Reproducibility link
```

---

## 📊 Ejemplo Completo

### INPUT
```json
{
  "query": "¿Estudios sobre CRISPR en arabidopsis con tecnología RNA-Seq?",
  "top_k": 10,
  "filters": {
    "organism": "Arabidopsis thaliana",
    "strategy": "RNA-Seq"
  }
}
```

### PROCESO

**Paso 1: LLM Query Generation**
```
Input:  "Estudios sobre CRISPR en arabidopsis con RNA-Seq"
Output: "((CRISPR OR CRISPR-Cas9 OR gene editing) AND 
         (arabidopsis OR Arabidopsis thaliana)) AND 
         (RNA-Seq OR transcriptome) AND [Strategy]"
```

**Paso 2: NCBI E-utilities**
```
POST https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi
├─ db: sra
├─ term: ((CRISPR)AND(arabidopsis))AND(RNA-Seq)
└─ retmax: 100

Response:
├─ <Count>15742</Count>
├─ <WebEnv>MCID_...</WebEnv>
└─ <QueryKey>1</QueryKey>

Time: 2.3 segundos
```

**Paso 3: Fetch Results**
```
POST https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
├─ webenv: MCID_...
├─ query_key: 1
├─ rettype: xml
└─ retmax: 100

Response: 100 estudios en XML format
Time: 1.2 segundos
```

**Paso 4: KB Enrichment**
```
Para cada resultado:
  organism: "Arabidopsis thaliana"
  └─ Lookup KB:
     ├─ experiments: 6,479,083
     ├─ studies: 351
     ├─ avg_samples: 18,000
     └─ related_studies: 42

SCORE = 0.94 (muy relevante al KB)
```

### OUTPUT
```json
{
  "query_original": "CRISPR en arabidopsis con RNA-Seq",
  "ncbi_query_generated": "((CRISPR)AND(arabidopsis))AND(RNA-Seq)",
  "result_metrics": {
    "ncbi_total_matches": 15742,
    "displayed_top_k": 10,
    "ncbi_search_time": 3.5,
    "enrichment_time": 0.2
  },
  "top_results": [
    {
      "rank": 1,
      "ncbi_id": "SRP123456",
      "title": "CRISPR-based transcriptome editing in A. thaliana",
      "organism": "Arabidopsis thaliana",
      "strategy": "RNA-Seq",
      "platform": "Illumina NextSeq",
      "study_type": "experimental",
      "publication_date": "2024-11-15",
      "kb_enrichment": {
        "organism_experiments": 6479083,
        "organism_studies": 351,
        "similar_crispr_studies": 42,
        "strategy_prevalence": 0.22,
        "relevance_score": 0.94
      },
      "ncbi_link": "https://www.ncbi.nlm.nih.gov/sra/SRP123456",
      "metadata": {
        "author": "...",
        "institution": "...",
        "samples_count": 18
      }
    },
    { ... 9 más ...}
  ],
  "summary_analysis": {
    "trend": "CRISPR usage in plants increasing",
    "top_organism": "Arabidopsis thaliana (98.2%)",
    "most_common_platform": "Illumina (76%)",
    "estimated_total_reads": "12.3 billion",
    "data_freshness": "within last 6 months"
  },
  "reproducibility": {
    "ncbi_search_url": "https://www.ncbi.nlm.nih.gov/sra/advanced?term=...",
    "copy_paste_query": "((CRISPR)AND(arabidopsis))AND(RNA-Seq)",
    "instructions": "Paste above query directly into NCBI SRA search"
  },
  "cache_info": {
    "from_cache": false,
    "cached_at": "2026-02-06T14:32:00Z",
    "ttl_days": 7
  }
}
```

---

## 🔨 Paso a Paso: Implementación

### Paso 1: NCBI E-utilities Wrapper (45 min)

```python
# phase3/api/ncbi_integration.py

import requests
import xml.etree.ElementTree as ET
from typing import Dict, List
from pathlib import Path

class NCBISearcher:
    def __init__(self):
        self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.db = "sra"
        self.email = "la.humada@gmail.com"  # NCBI email
        
    def search(self, query: str, retmax: int = 100) -> Dict:
        """Ejecutar búsqueda en NCBI SRA"""
        url = f"{self.base_url}/esearch.fcgi"
        params = {
            'db': self.db,
            'term': query,
            'retmax': retmax,
            'email': self.email,
            'tool': 'pyner_v3'
        }
        
        response = requests.get(url, params=params)
        response.raise_for_status()
        
        root = ET.fromstring(response.text)
        return {
            'count': int(root.find('Count').text),
            'webenv': root.find('WebEnv').text,
            'query_key': root.find('QueryKey').text
        }
    
    def fetch(self, webenv: str, query_key: str, start: int = 0, count: int = 100) -> List[Dict]:
        """Obtener resultados en detalle"""
        url = f"{self.base_url}/efetch.fcgi"
        params = {
            'db': self.db,
            'webenv': webenv,
            'query_key': query_key,
            'rettype': 'xml',
            'retstart': start,
            'retmax': count,
            'email': self.email,
            'tool': 'pyner_v3'
        }
        
        response = requests.get(url, params=params)
        response.raise_for_status()
        
        return self._parse_sra_xml(response.text)
    
    def _parse_sra_xml(self, xml_text: str) -> List[Dict]:
        """Parsear XML de SRA a formato JSON"""
        root = ET.fromstring(xml_text)
        results = []
        
        for package in root.findall('.//EXPERIMENT_PACKAGE'):
            result = {
                'id': package.find('.//PRIMARY_ID').text,
                'title': package.find('.//TITLE')?.text or "N/A",
                'organism': self._extract_organism(package),
                'strategy': self._extract_strategy(package),
                'platform': self._extract_platform(package),
                'sample_count': len(package.findall('.//SAMPLE')),
                'run_count': len(package.findall('.//RUN'))
            }
            results.append(result)
        
        return results
    
    # ... helper methods ...
```

### Paso 2: Query Generator (30 min)

```python
# phase3/api/query_generator.py

from phase3.api.ollama_integration import OllamaClient

class NCBIQueryGenerator:
    def __init__(self):
        self.ollama = OllamaClient()
    
    def generate_query(self, natural_language: str) -> str:
        """Generar query NCBI desde lenguaje natural"""
        
        prompt = f"""
        Eres experto en búsqueda NCBI SRA. Convierte esta pregunta natural
        a una query booleana NCBI-compatible.
        
        RESTRICCIONES:
        - Usa solo campos válidos: [Organism], [Strategy], [Platform], [All Fields]
        - Operadores: AND, OR, NOT
        - Wildcards: * para búsqueda parcial
        - NO uses comillas simples
        
        PREGUNTA: {natural_language}
        
        Devuelve SOLO la query, sin explicación.
        Ejemplo formato: ((CRISPR) AND (arabidopsis)) AND (RNA-Seq)
        """
        
        response = self.ollama.expand_query(prompt)
        return self._clean_query(response)
    
    def _clean_query(self, query_text: str) -> str:
        """Limpiar y validar query"""
        # Remover comillas problemas
        query_text = query_text.replace("'", "")
        # Validar booleano
        return query_text.strip()
```

### Paso 3: Caching Layer (30 min)

```python
# phase3/cache/cache_manager.py

import sqlite3
import hashlib
from datetime import datetime, timedelta
from pathlib import Path

class QueryCache:
    def __init__(self, db_path: str = "phase3/cache/queries.db"):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()
    
    def _init_db(self):
        """Crear tabla si no existe"""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute('''
                CREATE TABLE IF NOT EXISTS cached_queries (
                    query_hash TEXT PRIMARY KEY,
                    query_text TEXT,
                    ncbi_results JSON,
                    cached_at TIMESTAMP,
                    ttl_days INTEGER DEFAULT 7
                )
            ''')
            conn.commit()
    
    def store(self, query: str, results: str, ttl_days: int = 7):
        """Guardar query + resultados"""
        query_hash = hashlib.md5(query.encode()).hexdigest()
        with sqlite3.connect(self.db_path) as conn:
            conn.execute('''
                INSERT OR REPLACE INTO cached_queries
                VALUES (?, ?, ?, ?, ?)
            ''', (query_hash, query, results, datetime.now(), ttl_days))
            conn.commit()
    
    def retrieve(self, query: str) -> str:
        """Obtener del cache si es válido"""
        query_hash = hashlib.md5(query.encode()).hexdigest()
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute('''
                SELECT ncbi_results, cached_at, ttl_days
                FROM cached_queries
                WHERE query_hash = ?
            ''', (query_hash,))
            
            row = cursor.fetchone()
            if not row:
                return None
            
            results, cached_at, ttl_days = row
            cached_time = datetime.fromisoformat(cached_at)
            
            if datetime.now() - cached_time > timedelta(days=ttl_days):
                return None
            
            return results
```

### Paso 4: Result Enricher (45 min)

```python
# phase3/api/result_enricher.py

import json
from pathlib import Path

class ResultEnricher:
    def __init__(self, kb_path: str = "phase1/output/stage3_knowledge_base.json"):
        self.kb = self._load_kb(kb_path)
    
    def _load_kb(self, kb_path: str) -> Dict:
        """Cargar Knowledge Base"""
        with open(kb_path, 'r') as f:
            return json.load(f)
    
    def enrich_result(self, ncbi_result: Dict) -> Dict:
        """Enriquecer resultado NCBI con datos del KB"""
        
        organism = ncbi_result.get('organism', 'unknown')
        
        # Buscar en KB
        kb_data = self.kb['organisms'].get(organism, {})
        
        return {
            **ncbi_result,
            'kb_enrichment': {
                'organism_experiments': kb_data.get('count', 0),
                'organism_studies': kb_data.get('studies', 0),
                'relevance_score': self._calculate_score(ncbi_result, kb_data)
            }
        }
    
    def _calculate_score(self, ncbi_result: Dict, kb_data: Dict) -> float:
        """Calcular score de relevancia"""
        score = 0.5  # Base
        
        if kb_data.get('count', 0) > 1_000_000:
            score += 0.3  # Organismo bien estudiado
        
        if ncbi_result.get('strategy') == 'RNA-Seq':
            score += 0.1  # Strategy frecuente
        
        return min(score, 1.0)
```

### Paso 5: API Endpoint (30 min)

```python
# Añadir a phase3/api/main.py

@app.post("/search/real")
async def search_ncbi_with_enrichment(request: SearchRequest):
    """
    Phase 3b: Búsqueda real en NCBI + enriquecimiento
    """
    import time
    start_time = time.time()
    
    try:
        # 1. Generar query
        query_gen = NCBIQueryGenerator()
        ncbi_query = query_gen.generate_query(request.query)
        
        # 2. Verificar cache
        cache = QueryCache()
        cached = cache.retrieve(ncbi_query)
        if cached:
            results = json.loads(cached)
        else:
            # 3. Búsqueda NCBI
            searcher = NCBISearcher()
            search_result = searcher.search(ncbi_query, retmax=100)
            
            # 4. Obtener resultados
            ncbi_results = searcher.fetch(
                search_result['webenv'],
                search_result['query_key'],
                count=10
            )
            
            # 5. Enriquecer con KB
            enricher = ResultEnricher()
            enriched = [enricher.enrich_result(r) for r in ncbi_results]
            
            # 6. Cachear
            cache.store(ncbi_query, json.dumps(enriched))
            results = enriched
        
        return {
            "query_original": request.query,
            "ncbi_query_generated": ncbi_query,
            "results": results,
            "execution_time": time.time() - start_time
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
```

---

## ✅ Checklist de Implementación

- [ ] **Paso 1**: NCBI E-utilities Wrapper
  - [ ] search() method
  - [ ] fetch() method
  - [ ] XML parsing
  - [ ] Error handling

- [ ] **Paso 2**: Query Generator
  - [ ] LLM integration
  - [ ] Query validation
  - [ ] Filter application

- [ ] **Paso 3**: Cache Layer
  - [ ] SQLite setup
  - [ ] store() method
  - [ ] retrieve() method
  - [ ] TTL validation

- [ ] **Paso 4**: Result Enricher
  - [ ] KB loading
  - [ ] Score calculation
  - [ ] Data correlation

- [ ] **Paso 5**: API Integration
  - [ ] New endpoint `/search/real`
  - [ ] Pipeline orchestration
  - [ ] Response formatting

- [ ] **Testing**:
  - [ ] Unit tests (each component)
  - [ ] Integration tests (full pipeline)
  - [ ] Real NCBI queries
  - [ ] Performance benchmarks

- [ ] **Documentation**:
  - [ ] API docs (Swagger)
  - [ ] Usage examples
  - [ ] Deployment guide

- [ ] **Deployment**:
  - [ ] GitHub commit
  - [ ] Version bump (3.1)
  - [ ] Production testing

---

## 📊 Performance Expectations

```
Component                  Time      Notes
════════════════════════════════════════════
Query Generation (LLM)     50-200ms  Cached después
NCBI Search (esearch)      1-3s      Network
NCBI Fetch (efetch)        1-2s      Network
XML Parsing                100-300ms Local
KB Enrichment              50-200ms  Local
Caching                    < 5ms     SQLite
────────────────────────────────────
Total (first run)          3-8s      Network bound
Total (cached)             < 100ms   Local only
```

---

## 🚀 Próximos Pasos

1. Implementar los 5 componentes arriba
2. Testear con queries reales
3. Optimizar performance
4. Deployer v3.1
5. Feedback del usuario

¿Comenzamos con Paso 1?

