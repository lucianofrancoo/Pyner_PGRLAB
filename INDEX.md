# 🚀 Guía de Inicio Rápido - Pyner PGRLAB

**¿Eres nuevo en el proyecto? Empieza aquí.**

---

## ⚡ 5 Minutos para Entender Pyner

### Qué es Pyner?
Una herramienta para buscar estudios bioinformáticos en **NCBI SRA** usando IA local.

**Antes (sin Pyner):**
```
"Quiero estudios de arabidopsis bajo sequía"
→ Búsqueda manual en NCBI
→ 10,000+ resultados irrelevantes
→ 4 horas filtrando
```

**Después (con Pyner):**
```
python pyner.py "arabidopsis drought"
→ Query optimizada automática
→ 500 resultados altamente relevantes
→ 10 minutos analizando
```

---

## 📚 Documentación Esencial

| Documento | Tiempo | Para Quién |
|-----------|--------|-----------|
| **[README.md](README.md)** | 3 min | Todos (introducción) |
| **[ROADMAP.md](ROADMAP.md)** | 5 min | Contribuyentes (entender fases) |
| **[DEVELOPMENT_PLAN.md](DEVELOPMENT_PLAN.md)** | 20 min | Devs (detalles técnicos) |
| **[CONTRIBUTING.md](CONTRIBUTING.md)** | 15 min | Contribuyentes (cómo colaborar) |
| **[docs/KB_SCHEMA.md](docs/KB_SCHEMA.md)** | 15 min | Devs avanzados (estructuras) |

---

## 🎯 Donde Estamos (Febrero 2026)

### ✅ COMPLETADO
- ✅ Scripts iniciales (v0.1 y v0.2)
- ✅ Documentación completa
- ✅ Plan de 3 fases definido
- ✅ Repositorio GitHub activo

### 🟡 EN PROGRESO
- 🟡 Reclutamiento de contribuyentes para Fase 1
- 🟡 Diseño del Parser XML

### ⏳ EN COLA
- ⏳ Fase 1: Extracción de Knowledge Base (100k archivos)
- ⏳ Fase 2: Query Optimizer determinístico
- ⏳ Fase 3: Producción y v0.3 release

---

## 🤝 ¿Cómo Contribuir?

### Opción 1: Reportar un Bug (Fácil, 5 min)

```
1. Ve a: GitHub Issues
2. Clic: "New Issue"
3. Seleccionar: "Bug Report"
4. Completar el template

Ejemplos:
- "Parser XML falla con archivos > 10MB"
- "Query generator genera sintaxis inválida"
```

### Opción 2: Sugerir una Feature (Fácil, 5 min)

```
1. Ve a: GitHub Issues
2. Clic: "New Issue"
3. Seleccionar: "Feature Request"
4. Describir la idea

Ejemplos:
- "Soporte para múltiples organismos simultáneamente"
- "Caché para queries ya ejecutadas"
```

### Opción 3: Contribuir Código (Medio-Difícil, 5-20 horas)

**Para Principiantes:**
```
1. Lee: CONTRIBUTING.md
2. Elige tarea: Busca "Prioridad 1" en DEVELOPMENT_PLAN.md
3. Crea rama: git checkout -b feature/nombre
4. Implementa con tests
5. Pull request
```

**Para Avanzados:**
```
- Optimización de performance
- Integración con herramientas externas
- Machine learning para query ranking
```

---

## 📋 Tareas Disponibles (Prioridad)

### 🔴 URGENTE - Fase 1: Knowledge Base

**T1.1: Parser XML** (⭐⭐ Fácil-Medio)
- Requiere: Python, XML
- Tiempo: 20 horas
- Impacto: Crítico

**T1.2: Procesamiento Paralelo** (⭐⭐⭐ Medio)
- Requiere: Multiprocessing, Performance
- Tiempo: 15 horas
- Impacto: Alto

**T1.3: Normalización de Sinónimos** (⭐⭐⭐⭐ Difícil)
- Requiere: Bioinformatica, Regex
- Tiempo: 25 horas
- Impacto: Muy Alto

👉 **Quieres trabajar en una?** Comenta en [Issue de Fase 1](https://github.com/lucianofrancoo/Pyner_PGRLAB/issues)

---

## 🏗️ Arquitectura Simple

```
User Input
    ↓
[Normalización] ← "arabidopsis" → "arabidopsis_thaliana"
    ↓
[Query Builder] ← Usa Knowledge Base local
    ↓
[NCBI SRA Query] ← Búsqueda booleana optimizada
    ↓
[Resultados]
    ↓
[LLM Analysis] ← Opcional: análisis adicional
    ↓
[CSV Export]
```

**Fase 1:** Construir Knowledge Base  
**Fase 2:** Mejorar Query Builder  
**Fase 3:** Optimizar todo

---

## 💻 Instalación Local

### Requisitos
- Python 3.8+
- Git
- 16GB RAM (ideal 32GB)

### Setup en 5 min

```bash
# 1. Clonar
git clone https://github.com/lucianofrancoo/Pyner_PGRLAB.git
cd Pyner_PGRLAB

# 2. Instalar dependencias
pip install -r requirements.txt
pip install -r requirements-dev.txt

# 3. Descargar LLM (opcional)
ollama pull qwen2.5:14b

# 4. Tests
pytest tests/ -v

# ¡Listo!
```

---

## 📊 Estrutura de Carpetas

```
Pyner_PGRLAB/
├── README.md                  ← Empieza aquí (3 min)
├── ROADMAP.md                 ← Timeline del proyecto (5 min)
├── DEVELOPMENT_PLAN.md        ← Plan técnico detallado (20 min)
├── CONTRIBUTING.md            ← Cómo colaborar (15 min)
├── INDEX.md                   ← Este archivo
│
├── scripts/                   ← Scripts Python
│   ├── Pyner_search_v0.1.py
│   ├── Pyner_v0.1.py
│   └── Pyner_v0.2.py
│
├── docs/
│   └── KB_SCHEMA.md          ← Esquema técnico (15 min)
│
├── data/                      ← Datos (por crear)
│   └── kb/
│       ├── organisms_index.json     (Fase 1)
│       ├── strategies_index.json    (Fase 1)
│       └── ...
│
├── tests/                     ← Tests (por crear)
│   └── test_*.py
│
└── .github/
    └── ISSUE_TEMPLATE/       ← Templates GitHub
```

---

## 🎓 Flujo de Contribución (Git)

```bash
# 1. Fork en GitHub (arriba a la derecha)

# 2. Clonar tu fork
git clone https://github.com/TU_USER/Pyner_PGRLAB.git
cd Pyner_PGRLAB

# 3. Crear rama
git checkout -b feature/nombre-descriptivo

# 4. Hacer cambios + commits
git add .
git commit -m "Agrega feature X"

# 5. Push a tu fork
git push origin feature/nombre-descriptivo

# 6. Pull request en GitHub
# (Botón que aparecerá en tu fork)

# 7. Review y merge ✅
```

Ver detalles en: [CONTRIBUTING.md](CONTRIBUTING.md)

---

## 📞 Preguntas?

| Canal | Usar Para |
|-------|-----------|
| **GitHub Issues** | Bugs, features, preguntas técnicas |
| **GitHub Discussions** | Conversaciones, ideas, building en voz alta |
| **Email** | Privado: lucianofranco.a@gmail.com |

---

## 🏆 Contributors

Agradeciamientos a todos los que colaboran:

- 👤 **Luciano Franco** - Project Lead
- 👥 **Tú?** - Contribuyente (¡pronto!)

---

## 📈 Hitos Próximos

```
✅ Plan completado (Hoy)
→ T1.1: Parser XML (3-4 semanas)
→ T1.2: Procesamiento (1-2 semanas)
→ T1.3: Normalización (2-3 semanas)
→ Fase 2: Query Optimizer (4-5 semanas)
→ Fase 3: Production release (2-4 semanas)

**Objetivo Final:** v0.3 con 3-5x mejor precisión en Q3 2026
```

---

## 🚀 ¿Por Dónde Empiezo?

### Si eres **científico** buscando una herramienta:
1. Lee [README.md](README.md)
2. Instala con `pip install -r requirements.txt`
3. Prueba `python Pyner_v0.2.py`

### Si eres **programador** queriendo contribuir:
1. Lee [README.md](README.md) (3 min)
2. Lee [ROADMAP.md](ROADMAP.md) (5 min)
3. Lee [CONTRIBUTING.md](CONTRIBUTING.md) (15 min)
4. Elige una tarea [T1.1, T1.2 o T1.3]
5. Comenta en GitHub Issues "Quiero trabajar en T1.X"

### Si eres **estudiante** buscando proyecto:
Perfecto. Esto es:
- ✅ Open source
- ✅ Real world problem
- ✅ Bioinformatica
- ✅ Machine learning aplicado
- ✅ Excelente para portafolio

Comenta en issues diciendo que eres estudiante, te guiaremos.

---

## 📚 Materiales de Referencia

- [NCBI SRA Documentation](https://www.ncbi.nlm.nih.gov/sra)
- [BioPython Tutorial](https://biopython.org/wiki/Documentation)
- [Boolean Query Syntax](https://www.ncbi.nlm.nih.gov/books/NBK3827/)
- [Git Learning](https://git-scm.com/book/es/v2)

---

## ✨ Cheatsheet

| Acción | Comando |
|--------|---------|
| Ver estado | `git status` |
| Ver cambios | `git diff` |
| Actualizar | `git pull origin main` |
| Crear rama | `git checkout -b feature/X` |
| Cambiar rama | `git checkout branch-name` |
| Guardar cambios | `git add . && git commit -m "mensaje"` |
| Enviar | `git push origin feature/X` |
| Tests | `pytest tests/ -v` |
| Coverage | `pytest tests/ --cov=scripts` |

---

## 🎯 Resumen

| Meta | Status |
|------|--------|
| Documentación | ✅ Completa |
| Código base | ✅ Funcional |
| Tests | ⏳ Por hacer |
| Knowledge Base | ⏳ Por hacer (Fase 1) |
| Query Optimizer | ⏳ Por hacer (Fase 2) |
| v0.3 Release | ⏳ Q3 2026 |

---

**¡Bienvenido a Pyner! 🚀**

**Próximo paso:** Elige tu rol y comienza. ¿Preguntas? GitHub Issues. ¡Nos vemos!

---

*Última actualización: Febrero 6, 2026*  
*Versión: v0.3 Planning*
