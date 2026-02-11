# Sistema Interactivo de Generación de Queries NCBI - Guía de Uso

## Descripción General

El sistema Phase 3 ahora incluye un **modo interactivo sin aprendizaje dinámico** que permite:

1. **Generar queries** en español o inglés
2. **Validar** si el query generado es correcto
3. **Revisar** rápidamente el resultado

## Vocabulario Técnico

El sistema usa un único archivo de soporte:

- `Query_generator/phases/phase3/support_dictionary/technical_vocabulary.json`

Ahí están los tejidos, condiciones, alias de organismos y palabras clave de estrategias.

---

## Modos de Uso

### 1. **Modo Interactivo (por defecto)**

```bash
python3 main.py -i "tu consulta en español"
```

**Flujo:**
1. Sistema genera query
2. Muestra términos extraídos y query listo para NCBI

**Ejemplo:**
```bash
$ python3 main.py -i "proteoma de células tumorales"

📝 Input: proteoma de células tumorales
🔎 Extracted:
   Organism:   (none)
   Strategies: (none)
   Tissues:    cell
   Keywords:   tumor

✅ NCBI Query: ...
```

### 2. **Modo Rápido (sin interacción)**

```bash
python3 main.py -q "tu consulta"
```

Genera el query sin preguntar. Útil para automatización.

### 3. **Modo Servidor**

```bash
python3 main.py --server
# Inicia FastAPI en http://0.0.0.0:8000
```

### 4. **Ver Estadísticas**

```bash
python3 main.py --stats
```

Muestra cuántas traducciones y términos tienes almacenados.

---

## Consultas en Español - Ejemplos

```bash
# Búsqueda simple
python3 main.py -i "arabidopsis sequia raices"

# Búsqueda con estrategia
python3 main.py -i "expresión génica en trigo con RNA-Seq"

# Búsqueda proteómica
python3 main.py -i "proteoma hígado estrés hídrico"

# Búsqueda de metabolitos
python3 main.py -i "metaboloma raíces auxinas"
```

---

## Ubicación del Vocabulario Técnico

```
Query_generator/phases/phase3/support_dictionary/
└── technical_vocabulary.json
```

Puedes editar este archivo JSON para:
- Agregar sinónimos, estrategias o alias de organismos
- Extender tejidos/condiciones

---

## Integración con LLM

El sistema usa **Qwen 2.5:14b** (modelo multiidioma vía Ollama) para:

1. **Extracción inicial:** detecta organism, estrategia, tejido, condiciones
2. **Traducción de términos desconocidos:** si no está en el vocabulario técnico
3. **Validación:** verifica que organismos no se mezclen con tejidos incompatibles

---

## Mejoras Futuras

- [ ] CLI para editar el vocabulario técnico directamente
- [ ] Exportar/importar vocabulario técnico
- [ ] Validación de términos por expertos
- [ ] API para actualizar vocabulario técnico
- [ ] Búsqueda en Wikipedia/ontologías automáticas

---

## Troubleshooting

**P: El query sigue siendo incorrecto después de corregir**
R: Asegúrate de que el campo está siendo usado correctamente. El sistema regenera solo si es organismo/estrategia/terminología importante.

**P: Qwen tarda mucho en responder**
R: Normal (2-5 segundos). Usa `-q` para modo rápido sin interacción.

**P: No aparece el vocabulario técnico**
R: Verifícalo con `python3 main.py --stats`

---

## Comandos Rápidos

```bash
# Modo interactivo (default)
python3 main.py -i "tu consulta"

# Modo rápido
python3 main.py -q "tu consulta"

# Con servidor
python3 main.py --server

# Ver estadísticas
python3 main.py --stats

# Sin LLM (solo vocabulario técnico)
python3 main.py -i --no-llm "tu consulta"
```
