# 🔄 Prompt Flow → LLM → Output

## 1️⃣ INPUT

### Data reaching `analyze_paper()`:

```
Title:      "Transcriptome-based meta-analysis of drought stress regulatory genes in tomato."
Abstract:   "Plants possess various molecular defense systems to ward off biotic... [~500 chars]"
User Query: "("Solanum lycopersicum"[Organism] OR "tomato"[All Fields]) AND ("RNA-Seq"[All Fields])"
Full Text:  "Methods: Plant seeds were germinated... Results: 1,234 genes were..." [~10-12KB]
```

---

## 2️⃣ PROMPT CONSTRUCTION

The prompt sent to Ollama is a structured **large string**:

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
- Use empty arrays [] if category not found
- Use scientific names for organisms when possible
- For strategies: be comprehensive. Include ALL techniques mentioned in Methods/Results
- Grammar: "qRT-PCR" (not "qrt-pcr"), "RNA-Seq" (not "RNA-seq")

JSON:
```

### Prompt Size:
- **Without full text**: ~3-4 KB
- **With full text**: ~10-12 KB total
- **Dynamic Timeout**: Calculated based on length (~1 min per 2KB)

---

## 3️⃣ OLLAMA API CALL

```python
requests.post(
    "http://localhost:11434/api/generate",
    json={
        "model": "qwen3.5:9b",
        "prompt": f"""[FULL PROMPT FROM ABOVE]""",
        "stream": False,
        "options": {
            "temperature": 0.1,  # ← Low for consistent responses
            "top_p": 0.9
        }
    },
    timeout=315  # seconds (dynamically calculated)
)
```

---

## 4️⃣ OLLAMA RESPONSE

The LLM returns a string with embedded JSON:

```
{
  "relevance_score": 10,
  "organisms": ["Solanum lycopersicum", "Solanum tuberosum", "Nicotiana"],
  "tissues": ["root"],
  "conditions": ["drought stress", "water shortage", "heat"],
  "strategies": ["RNA-Seq", "qRT-PCR", "microarray", "statistical analysis"]
}
```

---

## 5️⃣ RESPONSE PARSING

The code **extracts JSON** from the response:

```python
def _parse_response(self, response: str) -> Dict:
    # 1. Look for first '{' and last '}'
    start = response.find('{')
    end = response.rfind('}') + 1
    
    # 2. Extract JSON substring
    json_str = response[start:end]
    
    # 3. Parse JSON to Python dictionary
    parsed = json.loads(json_str)
    
    # 4. Validate and sanitize
    result = {
        'relevance_score': int(parsed.get('relevance_score', 0)),
        'organisms': parsed.get('organisms', []),
        ...
    }
    
    # 5. Clamp relevance score between 0-10
    result['relevance_score'] = max(0, min(10, result['relevance_score']))
    
    return result
```

---

## 6️⃣ OUTPUT - CSV ROW

Converted to a row in the CSV table:

```csv
PMID,PMCID,Title,Relevance_Score,Is_Relevant,Organisms,Tissues,Conditions,Strategies,Year,Journal,DOI
41068586,PMC12513158,Transcriptome-based meta-analysis...,10,Yes,"Solanum lycopersicum ; Solanum tuberosum ; Nicotiana",root,"drought stress ; water shortage ; heat","RNA-Seq ; qRT-PCR ; microarray",2025,BMC plant biology,10.1186/s12870-025-07348-2
```

---

# 🎯 HOW IS RELEVANCE SCORE CALCULATED?

## LLM Criteria

The Qwen model analyzes the paper based on these criteria:

### 1️⃣ **Comparison with User Query**
Does the paper contain the key concepts (organisms, techniques, conditions) from the query?

### 2️⃣ **Content Overlap**
How much of the paper's scope matches the user's intent?

### 3️⃣ **Scoring Scale (0-10)**

- **10 = TOTALLY RELEVANT**: Perfect match for all criteria.
- **8-9 = HIGHLY RELEVANT**: Matches most criteria with high specificity.
- **5-7 = RELEVANT**: Basic match but may lack some specific techniques or conditions.
- **0-4 = IRRELEVANT / TANGENTIAL**: Does not meet the primary search interest.

Papers with score ≥ 5 are marked as "Yes" in the `Is_Relevant` column.
