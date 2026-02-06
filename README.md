# Pyner PGRLAB - Pipeline de Búsqueda Bioinformática Automática

**Pyner** es una herramienta para buscar, extraer y analizar datos de secuenciación genómico (RNA-seq, genomics) desde **NCBI SRA (Sequence Read Archive)** usando un LLM local para análisis automático.

## 📋 Descripción General

Este proyecto automatiza el proceso de:
1. **Buscar** estudios en NCBI SRA con palabras clave biológicas
2. **Extraer** metadatos estructurados (tejido, tratamientos, estrategia de secuenciación)
3. **Deduplicar** resultados para evitar procesar el mismo proyecto experimental dos veces
4. **Analizar** automáticamente con un LLM local (Ollama + Qwen2.5)
5. **Exportar** resultados en CSV para análisis posterior

## 🔧 Requisitos

### Software Requerido
- Python 3.8+
- [Ollama](https://ollama.ai/) - Para ejecutar modelos de LLM localmente
- NCBI Entrez Tools (BioPython)

### Instalación

```bash
# Clonar/descargar el proyecto
git clone https://github.com/lucianofrancoo/Pyner_PGRLAB.git
cd Pyner_PGRLAB

# Instalar dependencias Python
pip install biopython ollama pandas

# Descargar y ejecutar Ollama (en otra terminal)
# Visita https://ollama.ai para instrucciones de instalación
ollama pull qwen2.5:14b  # Descargar el modelo
ollama serve             # Mantener ejecutándose en background
```

### Configuración de Credenciales NCBI

Editar los scripts y reemplazar:
```python
Entrez.email = "tu_email@ejemplo.com"      # OBLIGATORIO
Entrez.api_key = "tu_api_key_ncbi"          # Opcional pero recomendado
```

Para obtener API key: https://www.ncbi.nlm.nih.gov/account/

---

## 📁 Scripts del Proyecto

### 1. **Pyner_search_v0.1.py** - Generador de Términos de Búsqueda
**Propósito:** Crear consultas booleanas optimizadas para NCBI SRA

**Uso:**
```bash
python3 Pyner_search_v0.1.py arabidopsis drought
python3 Pyner_search_v0.1.py solanum lycopersicum nitrogen
```

**¿Qué hace?**
- Recibe palabras clave como argumentos
- Envía a LLM local para procesamiento
- LLM genera:
  - Sinónimos biológicos relevantes
  - Consulta booleana compatible con NCBI
  - Frase natural describiendo la búsqueda

**Output:**
```json
{
  "natural_query": "RNA-seq studies of Arabidopsis under drought conditions",
  "esearch_query": "arabidopsis[Organism] AND (drought OR \"water stress\" OR dehydration)"
}
```

**Flujo:** `Keywords` → `LLM` → `JSON con consulta` → Copiar a Pyner_v0.2.py

---

### 2. **Pyner_v0.1.py** - Búsqueda Básica en GEO
**Propósito:** Búsqueda simple en NCBI GEO (Gene Expression Omnibus)

**Uso:**
```bash
python3 Pyner_v0.1.py
```

**¿Qué hace?**
1. Conecta a NCBI usando API de Entrez (BioPython)
2. Busca datasets en GEO con términos hardcodeados
3. Extrae el primer resultado
4. Analiza con LLM local
5. Imprime resultados en terminal

**Datos Extraidos:**
- Título del estudio
- Resumen
- Condiciones experimentales
- Si es serie temporal
- Tejidos estudiados

**Ventajas:**
- Rápido para pruebas
- GEO tiene datos más curados
- Bajo número de resultados

**Limitaciones:**
- Solo procesa primer resultado
- Busca hardcodeada (modificar en código)

---

### 3. **Pyner_v0.2.py** - Búsqueda Avanzada en SRA con Deduplicación
**Propósito:** Búsqueda completa y exportación a CSV (versión de producción)

**Uso:**
```bash
python3 Pyner_v0.2.py
```

**¿Qué hace?**
1. Busca en NCBI SRA (Sequence Read Archive)
2. Para cada resultado:
   - Parsea XML de metadatos
   - Extrae BioProject ID
   - **Si es nuevo BioProject:** procesa y analiza
   - **Si es duplicado:** omite (evita procesamiento redundante)
3. Envía cada estudio a LLM para análisis
4. Exporta resultados a CSV

**Datos Extraidos:**
- Título, organismo, tejido
- Estrategia de secuenciación (RNA-Seq, Genomics, etc.)
- Tipo de librería (PAIRED/SINGLE)
- BioProject y BioSample IDs
- Análisis LLM (condiciones, series temporales, tejidos)

**Output:**
Archivo CSV: `Pyner_SRA_arabidopsis_drought_unique.csv`

Columnas:
```
bioproject | study | organism | tissue | conditions | is_time_series | tissues_studied
```

**Ventajas:**
- Procesa múltiples resultados
- Deduplicación automática
- Exporta datos estructurados
- Práctico para análisis posterior

---

## 🔄 Flujo de Trabajo Típico

### Opción 1: Búsqueda Rápida (Pruebas)
```bash
# Generar consulta de búsqueda
python3 Pyner_search_v0.1.py arabidopsis drought

# Ver primer resultado (GEO)
python3 Pyner_v0.1.py
```

### Opción 2: Búsqueda Completa (Producción)
```bash
# 1. Generar consulta optimizada
python3 Pyner_search_v0.1.py arabidopsis drought

# 2. Copiar el "esearch_query" resultante
# 3. Editar Pyner_v0.2.py y reemplazar "search_term"
nano Pyner_v0.2.py  # O tu editor favorito

# 4. Ejecutar búsqueda completa
python3 Pyner_v0.2.py

# 5. Analizar resultados (CSV)
# Abrir en Excel, R, Python, etc.
```

---

## 📊 Archivos de Entrada/Salida

### Entrada
- **CSV:** `Pyner_SRA_arabidopsis_drought_unique.csv` (datos previos)
- **Variables en código:** Credenciales NCBI, términos de búsqueda

### Salida
- **CSV:** `Pyner_SRA_arabidopsis_drought_unique.csv` (actualizado)
- **Terminal:** Información detallada de procesamiento
- **Datos en RAM:** DataFrame de pandas en scripts

---

## 🤖 Integración con LLM (Ollama)

### Modelo Usado
**Qwen2.5:14b** - Modelo open-source optimizado para instrucciones

### Por qué Ollama?
- ✅ Privacidad: Corre localmente, sin enviar datos a internet
- ✅ Gratuito: Modelos open-source sin API costs
- ✅ Rápido: Una vez descargado, respuestas en segundos
- ✅ Flexible: Puedes cambiar de modelo fácilmente

### Cambiar Modelo
```bash
# Otros modelos disponibles (en Ollama)
ollama pull mistral:7b      # Más rápido
ollama pull llama2:13b      # Alternativa meta
```

Luego en el código, reemplazar `'qwen2.5:14b'` con el nuevo modelo.

---

## 🐛 Troubleshooting

### Error: "Connection refused" (Ollama)
```
❌ Error: [Errno 111] Connection refused
```
**Solución:**
```bash
# Asegúrate que ollama está corriendo
ollama serve
# En otra terminal, ejecuta tu script
```

### Error: "model not found"
```
❌ Error: model 'qwen2.5:14b' not found
```
**Solución:**
```bash
ollama pull qwen2.5:14b
```

### Error: "NCBI request failed"
```
❌ AuthenticationError
```
**Solución:**
- Verificar credenciales NCBI (email y API key)
- Consultar: https://www.ncbi.nlm.nih.gov/tools/primer-blast/

### Resultados vacíos
- Aumentar `retmax` en los scripts
- Verificar que la búsqueda tiene terminología correcta
- Probar con `Pyner_search_v0.1.py` primero

---

## 📈 Próximas Mejoras Planeadas

- [ ] GUI web (Streamlit o Flask)
- [ ] Soporte para múltiples organismos simultáneamente
- [ ] Análisis estadísticos de metadatos
- [ ] Descarga automática de datos crudos
- [ ] Pipelines de QC (control de calidad)

---

## 📝 Notas Técnicas

### Deduplicación (Pyner_v0.2.py)
El script evita procesar múltiples experimentos del mismo BioProject usando:
```python
bioprojects_vistos = set()  # Almacenar IDs únicos

if bioproject in bioprojects_vistos:
    continue  # Saltar si ya fue procesado

bioprojects_vistos.add(bioproject)  # Marcar como visto
```

Esto es crucial porque muchos papers publican múltiples experimentos bajo un mismo BioProject.

### XML Parsing
NCBI devuelve metadatos en XML. El script parsea:
```xml
<Experiment>
  <Study name="..."/>
  <Organism ScientificName="..." />
  <Library_descriptor>
    <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
    <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
  </Library_descriptor>
</Experiment>
```

---

## 👥 Contributing

Sugerencias y pull requests son bienvenidas.

---

## 📄 Licencia

Este proyecto está disponible para uso libre.

---

## 🔗 Referencias Útiles

- [NCBI Entrez E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
- [BioPython Tutorial](https://biopython.org/wiki/Documentation)
- [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra)
- [Ollama Project](https://ollama.ai)
- [Qwen2.5 Model Card](https://huggingface.co/Qwen/Qwen2.5-14B)

---

**Última actualización:** Febrero 2026  
**Versión:** 0.2  
**Autor:** Luciano Franco  
**Proyecto:** PGRLAB
