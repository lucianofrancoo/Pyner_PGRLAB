# 🔄 Flujo de Prompt → LLM → Output

## 1️⃣ ENTRADA (INPUT)

### Datos que llegan a `analyze_paper()`:

```
Title:      "Transcriptome-based meta-analysis of drought stress regulatory genes in tomato."
Abstract:   "Plants possess various molecular defense systems to ward off biotic... [~500 chars]"
User Query: "("Solanum lycopersicum"[Organism] OR "tomato"[All Fields]) AND ("RNA-Seq"[All Fields])"
Full Text:  "Methods: Plant seeds were germinated... Results: 1,234 genes were..." [~10-12KB]
```

---

## 2️⃣ CONSTRUCCIÓN DEL PROMPT

El prompt que se envía a Ollama es una **string gigante** estructurada así:

```
You are a scientific paper classifier. Analyze the following paper and extract structured information.

USER QUERY (for relevance): ("Solanum lycopersicum"[Organism] OR "tomato"[All Fields]) AND ("RNA-Seq"[All Fields])

PAPER TITLE: Transcriptome-based meta-analysis of drought stress regulatory genes in tomato.

PAPER ABSTRACT: Plants possess various molecular defense systems...

PAPER FULL TEXT (Methods & Results sections):
Methods: In this study, we collected samples from Solanum lycopersicum plants subjected to drought stress...
Results: RNA-Seq analysis revealed 1,234 drought-responsive genes...

[DETAILED INSTRUCTIONS FOR EXTRACTION]

YOUR TASK:
1. Evaluate if this paper is RELEVANT to the user query. Score from 0-10 (0=completely irrelevant, 10=perfectly relevant)
2. Extract ALL organisms mentioned...
3. Extract ALL tissues/organs mentioned...
4. Extract ALL experimental conditions/treatments...
5. Extract ALL experimental strategies/techniques used in the study...

RESPONSE FORMAT (JSON only, no explanation):
{
  "relevance_score": <number 0-10>,
  "organisms": ["organism1", "organism2"],
  "tissues": ["tissue1", "tissue2"],
  "conditions": ["condition1", "condition2"],
  "strategies": ["strategy1", "strategy2"]
}

IMPORTANT:
- Return ONLY valid JSON
- Use empty arrays [] if category not found (but try harder to find!)
- Use scientific names for organisms when possible
- For strategies: be comprehensive. Include ALL techniques mentioned in Methods/Results
- Grammar: "qRT-PCR" (not "qrt-pcr"), "RNA-Seq" (not "RNA-seq")

JSON:
```

### Tamaño del Prompt:
- **Sin full text**: ~3-4 KB
- **Con full text**: ~10-12 KB total
- **Timeout dinámico**: ~1 minuto por cada 2KB (267-353 segundos)

---

## 3️⃣ LLAMADA A OLLAMA

```python
requests.post(
    "http://localhost:11434/api/generate",
    json={
        "model": "qwen2.5:14b",
        "prompt": f"""[PROMPT COMPLETO DE ARRIBA]""",
        "stream": False,
        "options": {
            "temperature": 0.1,  # ← Bajo para respuestas consistentes
            "top_p": 0.9,
            "num_predict": 500   # ← Limita respuesta a 500 tokens
        }
    },
    timeout=315  # segundos (calculado dinámicamente)
)
```

---

## 4️⃣ RESPUESTA DE OLLAMA

El LLM devuelve una string con JSON embebido:

```
You are analyzing... [EXPLICACIÓN DEL LLM]

Here's the analysis:

{
  "relevance_score": 10,
  "organisms": ["Solanum lycopersicum", "Solanum tuberosum", "Nicotiana"],
  "tissues": ["root"],
  "conditions": ["drought stress", "water shortage", "heat"],
  "strategies": ["RNA-Seq", "qRT-PCR", "microarray", "statistical analysis"]
}
```

---

## 5️⃣ PARSING DE LA RESPUESTA

El código **extrae el JSON** de la respuesta:

```python
def _parse_response(self, response: str) -> Dict:
    # 1. Busca el primer '{' y último '}'
    start = response.find('{')        # Posición del primer {
    end = response.rfind('}') + 1     # Posición del último }
    
    # 2. Extrae substring JSON
    json_str = response[start:end]
    # json_str = '{"relevance_score": 10, "organisms": [...], ...}'
    
    # 3. Parsea JSON a diccionario Python
    parsed = json.loads(json_str)
    
    # 4. Valida y sanitiza
    result = {
        'relevance_score': int(parsed.get('relevance_score', 0)),  # Asegurar int
        'organisms': parsed.get('organisms', []),                   # Asegurar lista
        'tissues': parsed.get('tissues', []),
        'conditions': parsed.get('conditions', []),
        'strategies': parsed.get('strategies', [])
    }
    
    # 5. Clamp relevance score entre 0-10
    result['relevance_score'] = max(0, min(10, result['relevance_score']))
    
    return result
```

**Output:**
```python
{
    'relevance_score': 10,
    'organisms': ['Solanum lycopersicum', 'Solanum tuberosum', 'Nicotiana'],
    'tissues': ['root'],
    'conditions': ['drought stress', 'water shortage', 'heat'],
    'strategies': ['RNA-Seq', 'qRT-PCR', 'microarray', 'statistical analysis']
}
```

---

## 6️⃣ SALIDA (OUTPUT) - FILA CSV

Se convierte a fila CSV:

```csv
PMID,PMCID,Title,Relevance_Score,Is_Relevant,Organisms,Tissues,Conditions,Strategies,Year,Journal,DOI,Abstract_Preview
41068586,PMC12513158,Transcriptome-based meta-analysis...,10,Yes,"Solanum lycopersicum ; Solanum tuberosum ; Nicotiana",root,"drought stress ; water shortage ; heat","RNA-Seq ; qRT-PCR ; microarray ; statistical analysis",2025,BMC plant biology,10.1186/s12870-025-07348-2,Plants possess various...
```

---

# 🎯 ¿CÓMO SE CALCULA EL SCORE DE RELEVANCIA?

## Criterios que usa el LLM

El Qwen LLM analiza el paper usando estos criterios (del prompt):

### 1️⃣ **Comparación con Query**
```
USER QUERY: ("Solanum lycopersicum"[Organism] OR "tomato") AND (RNA-Seq)

¿Aparecen estos conceptos en el paper?
- Solanum lycopersicum / tomato? ✓ (Title + Abstract + Methods)
- RNA-Seq? ✓ (Title + Abstract + Methods)

Puntuación: 10/10 (Perfectamente relevante)
```

### 2️⃣ **Coincidencia de Contenido**
```
Query busca:  Tomato + RNA-Seq + drought genes
Paper trata:  Tomato + RNA-Seq (completo) + drought stress genes

Overlap = 100% → Score alto
```

### 3️⃣ **Completitud del Match**
```
Query pide 2 cosas:
  1. Solanum lycopersicum / tomato → ✓ Encontrado
  2. RNA-Seq → ✓ Encontrado

2/2 coincidencias → 10/10
```

### 4️⃣ **Escala de Puntuación (de 0-10)**

```
10 = TOTALMENTE RELEVANTE
   - Cumple con TODOS los criterios de la query
   - Organismos, técnicas y condiciones exactas
   - Ejemplo: Paper sobre "tomato drought RNA-Seq"

9 = MUY RELEVANTE
   - Cumple la mayoría de criterios
   - Técnicas similares pero no idénticas
   - Puede tener más organismos además del solicitado

8 = RELEVANTE
   - Cumple criterios básicos
   - Pero con algunas variaciones
   - Busca "tomato RNA-Seq" pero paper tiene "tomato qRT-PCR"

7 = PARCIALMENTE RELEVANTE
   - Cumple 1 de 2 criterios principales
   - "tomato" ✓ pero sin RNA-Seq, solo microarray

5 = POCO RELEVANTE (umbral por defecto)
   - Tangencialmente relacionado
   - Puede mencionar tomato pero en contexto diferente

0 = IRRELEVANTE
   - Sin relación con la query
```

---

# 📊 EJEMPLO COMPLETO DEL FLUJO

## Entrada
```
Title:    "Melatonin Improves Drought Stress Tolerance of Tomato..."
Query:    (tomato OR Solanum) AND (drought OR stress)
Abstract: "Tomato is sensitive to drought... melatonin improved..." [500 chars]
```

## Procesamiento LLM
```
¿Hay "tomato" o "Solanum"? SÍ → +5 puntos
¿Hay "drought" o "stress"? SÍ → +4 puntos
¿Es paper original sobre este tema? SÍ → +1 punto
Score = 9/10 (MUY RELEVANTE)
```

## Output
```csv
...,35204192,PMC8868175,"Melatonin Improves...",9,Yes,"Solanum lycopersicum ; tomato","root ; leaf","drought stress ; water withholding",...
```

---

# ⚙️ PARÁMETROS DE CONFIGURACIÓN

| Parámetro | Valor | Efecto |
|-----------|-------|--------|
| **temperature** | 0.1 | Muy bajo = respuestas consistentes (no creativas) |
| **top_p** | 0.9 | Limita diversidad de tokens |
| **num_predict** | 500 | Max 500 tokens de respuesta |
| **Relevance Threshold** | 5 | Papers con score ≥5 se marcan como "Yes" |
| **Max Timeout** | 353s | Para prompts muy largos (12KB+) |

---

# 🔍 EJEMPLO DE CÁLCULO DE SCORE

### Paper 1: "Transcriptome-based meta-analysis of drought stress genes in tomato"

```
Query: (tomato) AND (RNA-Seq)

Análisis LLM:
├─ Contiene "tomato"? ✓ SÍ (Title + Abstract + Methods + Results)
├─ Contiene "RNA-Seq"? ✓ SÍ (Abstract + Methods)
├─ Is relevancia = 100%? ✓ SÍ
└─ Score = 10/10 (PERFECTO)

Output CSV: Relevance_Score=10, Is_Relevant=Yes
```

### Paper 2: "Melatonin Improves Drought Stress in Tomato"

```
Query: (tomato) AND (RNA-Seq)

Análisis LLM:
├─ Contiene "tomato"? ✓ SÍ (Title + Abstract)
├─ Contiene "RNA-Seq"? ✗ NO (pero menciona análisis biochemical)
├─ Relevancia parcial = 70%
└─ Score = 8/10 (MUY RELEVANTE)

Output CSV: Relevance_Score=8, Is_Relevant=Yes (≥5)
```

### Paper 3: "Bacterial Wilt in Tomato"

```
Query: (tomato) AND (RNA-Seq)

Análisis LLM:
├─ Contiene "tomato"? ✓ SÍ
├─ Contiene "RNA-Seq"? ✗ NO (pero usa metaRNAseq)
├─ Relevancia = 60%
└─ Score = 6/10 (RELEVANTE)

Output CSV: Relevance_Score=6, Is_Relevant=Yes (≥5 pero cercano)
```

---

# 🎪 RESUMEN VISUAL DEL FLUJO COMPLETO

```
┌──────────────────────────────────────────────────────────────┐
│ ENTRADA: Paper + Query                                       │
│ (PMID: 41068586, Title, Abstract, Full Text, User Query)     │
└────────────────────┬─────────────────────────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────────────────────────┐
│ CONSTRUCCIÓN DE PROMPT (~10-12 KB)                           │
│ - Instructions for LLM                                       │
│ - Query exacta del usuario                                   │
│ - Title + Abstract + Full Text completo                      │
│ - Formato JSON esperado                                      │
└────────────────────┬─────────────────────────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────────────────────────┐
│ OLLAMA API CALL                                              │
│ - Model: qwen2.5:14b                                         │
│ - Temperature: 0.1 (determinístico)                          │
│ - Timeout: 315 segundos (dinámico)                           │
└────────────────────┬─────────────────────────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────────────────────────┐
│ OLLAMA RESPONSE (JSON embebido en texto)                     │
│ {                                                            │
│   "relevance_score": 10,                                     │
│   "organisms": ["Solanum lycopersicum", "Solanum tuberosum"],│
│   "tissues": ["root"],                                       │
│   "conditions": ["drought stress", "water shortage"],        │
│   "strategies": ["RNA-Seq", "qRT-PCR", "microarray"]        │
│ }                                                            │
└────────────────────┬─────────────────────────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────────────────────────┐
│ PARSING & VALIDATION                                         │
│ - Extrae JSON de la respuesta                                │
│ - Valida tipos (int, list, etc.)                             │
│ - Clamp score 0-10                                           │
│ - Asegura arrays para estrategias                            │
└────────────────────┬─────────────────────────────────────────┘
                     │
                     ▼
┌──────────────────────────────────────────────────────────────┐
│ OUTPUT: Fila CSV                                             │
│ PMID | PMCID | Title | Score | Is_Relevant | Organisms...   │
│ 41068586|PMC12513158|Transcriptome...|10|Yes|Solanum...     │
└──────────────────────────────────────────────────────────────┘
```
