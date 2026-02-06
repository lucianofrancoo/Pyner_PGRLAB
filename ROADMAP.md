# 🗺️ Roadmap - Pyner PGRLAB

## Visión a Largo Plazo

Convertir Pyner en **la herramienta más precisa para búsqueda bioinformática** aprovechando datos locales de NCBI SRA.

```
User ─→  Natural Language  ─→  Query Optimizer (KB)  ─→  NCBI SRA  ─→  Precise Results
         "arabidopsis          + Machine Learning          búsqueda      3-5x mejor
          drought"              + NCBI Knowledge Base       optimizada    precisión
```

---

## 📅 Roadmap 2026

### Q1 2026 (Enero-Marzo)

#### 🟡 FASE 1: Knowledge Base Extraction [EN PROGRESO]
```
Week 1-2:  Parser XML robusto                          ⏳ Por iniciar
Week 3-4:  Procesar 100k archivos                      ⏳ Por iniciar
Week 5-6:  Normalización de sinónimos                  ⏳ Por iniciar
Week 7-8:  Validación y documentación                  ⏳ Por iniciar

Dependencia: ✅ Plan documentado
Entregable:  data/kb/ (JSON con índices)
Status:      🟡 Planeado
```

**Tareas Abiertas:**
- [ ] T1.1: Parser XML (BUSCA CONTRIBUYENTE)
- [ ] T1.2: Procesamiento paralelo (BUSCA CONTRIBUYENTE)
- [ ] T1.3: Normalización (BUSCA CONTRIBUYENTE, nivel avanzado)

---

### Q2 2026 (Abril-Junio)

#### 🔴 FASE 2: Query Optimizer [BLOQUEADO POR FASE 1]
```
Week 1-2:  Query Builder determinístico              ⏳ Bloqueado
Week 3-4:  Test suite exhaustivo (50+ casos)         ⏳ Bloqueado
Week 5-6:  Integración con Pyner_search_v0.1.py      ⏳ Bloqueado
Week 7-8:  Validación vs LLM baseline                ⏳ Bloqueado

Dependencia: 🟡 FASE 1 debe terminar
Entregable:  scripts/query_optimizer.py
Status:      ⏳ En cola
```

**Tareas Abiertas:**
- [ ] T2.1: Query Builder (ESPERA FASE 1)
- [ ] T2.2: Test suite (ESPERA FASE 1)

---

### Q3-Q4 2026 (Julio-Diciembre)

#### 🔴 FASE 3: Optimization & Release [BLOQUEADO POR FASE 2]
```
Week 1-2:  Escalabilidad a 1.3M archivos              ⏳ Bloqueado
Week 3-4:  Evaluación exhaustiva de precisión        ⏳ Bloqueado
Week 5-6:  Performance tuning                        ⏳ Bloqueado
Week 7-8:  Release Pyner v0.3 + Documentación        ⏳ Bloqueado

Dependencia: 🔴 FASE 2 debe terminar
Entregable:  Pyner v0.3 Release
Status:      ⏳ En cola
```

**Tareas Abiertas:**
- [ ] T3.1: Escalabilidad (ESPERA FASE 2)
- [ ] T3.2: Evaluación (ESPERA FASE 2)

---

## 🎯 Hitos Principales

| Milestone | Fecha | Status | Blocker |
|-----------|-------|--------|---------|
| 📋 Plan documentado | 2026-02-06 | ✅ Completado | Ninguno |
| 🧬 KB Fase 1 (100k) | 2026-03-31 | ⏳ En progreso | Contribuyentes |
| 🤖 Query Optimizer v1 | 2026-05-31 | ⏳ Espera Fase 1 | Fase 1 |
| 📈 Scalability (1.3M) | 2026-08-31 | ⏳ Espera Fase 2 | Fase 2 |
| 🚀 Release v0.3 | 2026-10-31 | ⏳ Espera Fase 3 | Fase 3 |

---

## 📊 Capacidades por Versión

### v0.1 (Actual)
```
✅ Búsqueda simple en NCBI GEO
✅ LLM genérico para generación de queries
✅ Análisis básico con Ollama
❌ Sinónimos contextualizados
❌ Deduplicación
❌ Exportación en CSV
```

### v0.2 (Actual)
```
✅ Búsqueda en NCBI SRA
✅ Deduplicación por BioProject
✅ Exportación en CSV
✅ Análisis con LLM
❌ Query optimization
❌ Knowledge Base
❌ Precisión mejorada
```

### v0.3 (Target)
```
✅ Todo v0.2 +
✅ Knowledge Base de 1.3M estudios
✅ Query Optimizer determinístico
✅ 3-5x mejor precisión
✅ <500ms para generar queries
✅ Documentación completa
✅ Test suite (70%+ coverage)
✅ API pública estable
```

---

## 🚀 Quick Start para Contribuyentes

### Prioridad 1: NECESITAMOS AYUDA (Fase 1)

**Habilidad:** Python + XML  
**Tiempo:** ~20 horas  
**Impacto:** Crítico

```bash
# Tarea: Implementar Parser XML
# Issue: #TODO-1.1

# Qué necesitamos:
1. Function para parsear experiment.xml
2. Function para parsear sample.xml
3. Function para parsear run.xml
4. Tests para 100+ archivos

# Repo: https://github.com/lucianofrancoo/Pyner_PGRLAB
# Docs: CONTRIBUTING.md
```

---

### Prioridad 2: Procesamiento Paralelo (Fase 1)

**Habilidad:** Python + Multiprocessing/Asyncio  
**Tiempo:** ~15 horas  
**Impacto:** Alto

```bash
# Tarea: Script para procesar 100k archivos en paralelo
# Issue: #TODO-1.2

# Toma el parser XML (Prioridad 1)
# Optimiza para velocidad
# Genera índices JSON
```

---

### Prioridad 3: Normalización (Fase 1)

**Habilidad:** Python + Bioinformatica + Regex  
**Tiempo:** ~25 horas  
**Impacto:** Medio

```bash
# Tarea: Normalizar sinónimos de organismos
# Issue: #TODO-1.3

# Ejemplos:
# "arabidopsis" ─→ "arabidopsis_thaliana"
# "Arabidopsis thaliana" ─→ "arabidopsis_thaliana"
# "ath" ─→ "arabidopsis_thaliana"
```

---

## 💡 Oportunidades Futuras (Post v0.3)

```
v0.4 (2027-Q1):
  - Web UI (Streamlit/Flask)
  - API REST pública
  - Multi-lenguaje (ES/EN)

v0.5 (2027-Q2):
  - Machine learning para ranking de resultados
  - Recomendaciones automáticas
  - Integración con Galaxy

v1.0 (2027-Q3):
  - Soporte para GEO, SRA, ArrayExpress
  - Deep learning para análisis de papers
  - Publicación en PyPI
```

---

## 📊 Métricas de Éxito

```
       Baseline    Target    Mejora
       ────────────────────────────
Query Precision:  65%       85%       +31%
Results Relevance: -40%     -10%      +75%
Generation Time:  3-5s      <0.5s     6-10x
```

---

## 🤝 Cómo Involucrarse

### 1️⃣ Reportar Bugs
```bash
GitHub Issues → Bug Report template
Proporciona: ambiente, pasos, error completo
```

### 2️⃣ Sugerir Features
```bash
GitHub Issues → Feature Request template
Proporciona: problema, solución, impacto
```

### 3️⃣ Contribuir Código
```bash
Fork → Feature branch → Pull Request
Ver: CONTRIBUTING.md para detalles
```

### 4️⃣ Mejorar Docs
```bash
README.md, KB_SCHEMA.md, CONTRIBUTING.md
Mejoras de claridad siempre bienvenidas
```

---

## 📞 Contacto

- **GitHub:** [Pyner_PGRLAB](https://github.com/lucianofrancoo/Pyner_PGRLAB)
- **Issues:** [Issues](https://github.com/lucianofrancoo/Pyner_PGRLAB/issues)
- **Email:** lucianofranco.a@gmail.com
- **Project Lead:** Luciano Franco

---

## 📜 Changelog

## 📜 Changelog

| Versión | Fecha | Cambios |
|---------|-------|---------|
| v0.2 | 2026-02-06 | Scripts iniciales comentados, README completo |
| v0.3 (WIP) | 2026-02-06 | Plan de desarrollo + KB extraction |
| v1.0 (Planned) | 2027-Q3 | Release público estable |

---

## 🙏 Agradecimientos

Especialde a:
- 🧬 **PGRLAB** - Por los datos iniciales
- 🧪 **Comunidad científica** - Por feedback

---

**Última actualización:** Febrero 6, 2026  
**Status General:** 🟡 En planificación  
**Next:**  Fase 1 iniciará cuando contribuyentes se unan 🚀
