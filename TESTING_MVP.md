# 🚀 PYNER: Query Generation MVP (Sesgado - 500K)

## ⚠️ IMPORTANTE
Este es un MVP que funciona **con información sesgada** (solo 500K archivos, falta 85%).
Es una demostración funcional de la generación de queries.
**Phase 3b añadirá NCBI API integration para datos reales**.

---

## 🎯 Qué Funciona Ahora

✅ Generación de queries desde lenguaje natural
✅ Búsqueda vectorial (FAISS) en 211 queries
✅ Mapeo de resultados a Knowledge Base
✅ API REST en localhost:8000
✅ Múltiples endpoints funcionales

---

## 📋 Cómo Usar

### 1️⃣ Verificar que API está corriendo

```bash
# En terminal 1, el API debería estar corriendo:
curl http://localhost:8000/

# Debería responder:
# {"status": "ok", "service": "Pyner Semantic Search Phase 3", ...}
```

### 2️⃣ Test Rápido (Validación)

```bash
python3 phase3/scripts/quick_test.py
```

**Output esperado:**
```
✅ ALL TESTS PASSED!
```

### 3️⃣ Test de Queries (Simple)

```bash
# Ejemplo simple
python3 test_queries.py "virus humanos"

# Modo interactivo
python3 test_queries.py --interactive
```

**Output:**
```
❓ Query: "virus humanos"
✅ API está disponible

📊 Resultados (Top-5):
  1. [0.454] ██... Research studies on Vibrio cholerae
  2. [0.443] ██... Molecular mechanisms of COVID-19 in Homo sapiens
  ...
```

### 4️⃣ Test a través de API (cURL)

```bash
# Búsqueda simple
curl -X POST http://localhost:8000/search \
  -H "Content-Type: application/json" \
  -d '{
    "query": "CRISPR en plantas",
    "top_k": 5,
    "expand": false
  }'

# Con query expansion (Ollama)
curl -X POST http://localhost:8000/search \
  -H "Content-Type: application/json" \
  -d '{
    "query": "expresión génica en plantas",
    "top_k": 5,
    "expand": true
  }'

# Ver estadísticas
curl http://localhost:8000/stats
```

---

## 📊 Pipeline Actual (Sesgado)

```
Natural Language Query
        ↓
  FAISS Search (211 indexed)
        ↓
  Vector Similarity Match
        ↓
  Knowledge Base Lookup (500K)
        ↓
  JSON Response (Sesgado)
```

---

## 🎯 Ejemplos de Queries Funcionales

```bash
python3 test_queries.py "virus que infectan humanos"
python3 test_queries.py "CRISPR en plantas"
python3 test_queries.py "bacteria del suelo"
python3 test_queries.py "expresión génica"
python3 test_queries.py "antibiótesistencia"
```

---

## ⚠️ Limitaciones (Sesgado)

| Limitación | Impacto |
|-----------|---------|
| Solo 500K de 3.6M archivos | 0.014% cobertura |
| Datos históricos | No actualizado |
| Sin NCBI real-time | Falta ~15K+ resultados |
| KB estático | Falsos negativos |

**Solución: Phase 3b (En desarrollo)**

---

## 📁 Archivos Nuevos

```
phase3/
├─ scripts/
│  ├─ quick_test.py           # Teste rápido de componentes
│  ├─ demo_query_generator.py  # Demo completo (más largo)
│  └─ __init__.py
│
└─ (Los que ya existían)
    ├─ api/main.py            # API FastAPI
    ├─ api/ollama_integration.py
    ├─ config.py
    └─ logs/

Raíz del proyecto:
├─ test_queries.py            # Script simple para testear
├─ ARCHITECTURE_DECISION.md    # Por qué hacer Phase 3b
├─ IMPLEMENTATION_ROADMAP.md   # Cómo implementar Phase 3b
└─ RESULTS_500K_FILES.md       # Resultados extraídos
```

---

## 🔄 Flujo de Uso Recomendado

```
1. Verifica API corriendo:
   curl http://localhost:8000/

2. Test rápido:
   python3 phase3/scripts/quick_test.py

3. Prueba queries:
   python3 test_queries.py "tu pregunta"

4. Modo interactivo:
   python3 test_queries.py --interactive

5. Ver detalles técnicos:
   Lee ARCHITECTURE_DECISION.md para entender por qué
   falta Phase 3b
```

---

## 📊 Performance

| Métrica | Valor |
|---------|-------|
| Latencia total | 10-250ms |
| Vector search | < 5ms |
| API response | < 1ms (cached) |
| Throughput | 10-15 queries/sec |

---

## 🚀 Próximo Paso: Phase 3b

Para hacer esto **productivo** (sin sesgos):

```
Phase 3b: NCBI API Integration (3 horas)
├─ E-utilities wrapper (45 min)
├─ Caching layer (30 min)
├─ Result enricher (45 min)
├─ Query generator (30 min)
└─ New /search/real endpoint (30 min)

Resultado: Sistema profesional + real-time NCBI
```

Ver: `IMPLEMENTATION_ROADMAP.md`

---

## 📝 Notas

- Este MVP demuestra que el pipeline funciona
- El sesgado es **aceptable para demostración**
- Los resultados son **científicamente sesgados** (no usar en papers sin Phase 3b)
- Phase 3b es **crítico para producción**

---

## ✅ Checklist: Hoy Hicimos

- [x] Quick test script
- [x] Demo query generator
- [x] Simple test CLI (test_queries.py)
- [x] Validación end-to-end
- [x] Documentación MVP
- [x] Ejemplos funcionales

## ⏳ Próximo: Phase 3b

- [ ] NCBI E-utilities wrapper
- [ ] Caching layer
- [ ] Result enricher
- [ ] Real-time search
- [ ] Production deployment

---

**Status:** MVP Demostración Funcional ✅  
**Fecha:** 2026-02-06  
**Versión:** 3.0-mvp (sesgado)

