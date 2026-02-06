# 🤝 Guía de Contribución - Pyner PGRLAB

¡Gracias por tu interés en contribuir a Pyner! Este documento explica cómo participar en el proyecto.

---

## 📋 Índice

1. [Empezando](#empezando)
2. [Flujo de Trabajo (Git)](#flujo-de-trabajo-git)
3. [Estándares de Código](#estándares-de-código)
4. [Testing](#testing)
5. [Documentación](#documentación)
6. [Proceso de PR](#proceso-de-pr)
7. [Tareas Disponibles](#tareas-disponibles)

---

## 🚀 Empezando

### Prerequisitos
```bash
# Clonar el repositorio
git clone https://github.com/lucianofrancoo/Pyner_PGRLAB.git
cd Pyner_PGRLAB

# Instalar dependencias
pip install -r requirements.txt
pip install -r requirements-dev.txt

# Verifica Ollama está ejecutándose (si necesitas LLM)
ollama serve &

# Tests básicos
pytest tests/ -v
```

### Configurar Git
```bash
git config user.name "Tu Nombre"
git config user.email "tu.email@ejemplo.com"
```

---

## 🔄 Flujo de Trabajo (Git)

### 1. Crear rama para tu feature

```bash
# Actualizar main
git fetch origin
git checkout main
git pull origin main

# Crear rama descriptiva
git checkout -b feature/nombre-del-feature

# Ejemplos buenos de nombres:
# - feature/extract-metadata-phase1
# - bugfix/xml-parser-edge-cases
# - docs/add-kb-schema
# - enhancement/optimize-query-builder
```

### 2. Hacer cambios

```bash
# Editar archivos
nano scripts/mi_script.py

# Ver qué cambió
git status
git diff

# Agregar cambios
git add scripts/mi_script.py

# O agregar todo (cuidado)
git add -A

# Commit con mensaje descriptivo
git commit -m "Agregar parser para SAMPLE_ATTRIBUTES (Fase 1)"
```

### 3. Push a GitHub

```bash
# Enviar branch
git push origin feature/nombre-del-feature

# En GitHub, verás un botón para crear Pull Request
```

### 4. Pull Request

- Ir a https://github.com/lucianofrancoo/Pyner_PGRLAB/pulls
- Clic en "New Pull Request"
- Seleccionar tu branch → main
- Escribir descripción clara (ver template más abajo)
- Submit

### 5. Review y Merge

```bash
# Después que se apruebe:
# ✅ Los maintainers mergearan tu PR
# ✅ Tu branch se eliminará automáticamente

# Actualizar tu local
git fetch origin
git checkout main
git pull origin main

# Opcional: eliminar rama local
git branch -d feature/nombre-del-feature
```

---

## 📝 Estándares de Código

### Style Guide (PEP 8)

```python
# ✅ CORRECTO
def extract_organism_data(xml_file: str) -> dict:
    """
    Extract organism information from XML file.
    
    Args:
        xml_file: Path to experiment.xml file
        
    Returns:
        Dictionary with organism metadata
    """
    organism_data = {}
    
    try:
        tree = ET.parse(xml_file)
        # ... resto del código
    except ET.ParseError as e:
        logger.error(f"Error parsing {xml_file}: {e}")
        raise
    
    return organism_data


# ❌ INCORRECTO (ejemplos a evitar)
def extract_organism_data(f):  # Nombre poco descriptivo
    # Sin docstring
    d = {}  # Variable poco descriptiva
    try:
        t = ET.parse(f)
    except:  # Excepción genérica
        pass  # Silent failure
    return d
```

### Comentarios y Docstrings

```python
# Usar docstrings estilo Google
def function_name(param: type) -> type:
    """
    Una línea de descripción breve.
    
    Descripción más detallada si es necesaria,
    explicando qué hace y por qué.
    
    Args:
        param (type): Descripción del parámetro
        
    Returns:
        type: Descripción del retorno
        
    Raises:
        ValueError: Cuándo y por qué se lanza
        
    Example:
        >>> result = function_name(value)
        >>> print(result)
        expected_output
    """
    pass


# Comments inline para lógica compleja
# Pero NO para código obvio
items = []  # ✅ Necesario si la razón no es obvia
for item in collection:
    items.append(item)  # ❌ Este es obvio, no necesita comentario
```

### Imports y Orden

```python
# Orden: 1) stdlib, 2) terceros, 3) locales

# Stdlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

# Terceros
import numpy as np
import pandas as pd
from biopython import Entrez

# Locales
from scripts.utils import normalize_name
from data.kb import load_knowledge_base
```

---

## ✅ Testing

### Escribir Tests

```python
# tests/test_extract_metadata.py

import pytest
from scripts.extract_metadata import extract_organism_data


class TestExtractOrganism:
    """Test suite para extract_organism_data"""
    
    def test_valid_xml_parsing(self):
        """Should correctly parse valid sample.xml"""
        result = extract_organism_data("tests/fixtures/sample.xml")
        assert result["organism"] == "Bacillus subtilis"
        assert result["taxon_id"] == "645657"
    
    def test_missing_file(self):
        """Should raise FileNotFoundError for missing file"""
        with pytest.raises(FileNotFoundError):
            extract_organism_data("nonexistent.xml")
    
    def test_malformed_xml(self):
        """Should handle malformed XML gracefully"""
        with pytest.raises(ET.ParseError):
            extract_organism_data("tests/fixtures/malformed.xml")
```

### Ejecutar Tests

```bash
# Run all tests
pytest tests/ -v

# Run specific test file
pytest tests/test_extract_metadata.py -v

# Run specific test
pytest tests/test_extract_metadata.py::TestExtractOrganism::test_valid_xml_parsing -v

# With coverage
pytest tests/ --cov=scripts --cov-report=html

# Coverage debe ser mínimo 70%
```

### Coverage Report

```bash
# Generar reporte HTML
pytest tests/ --cov=scripts --cov-report=html

# Abrir reporte
open htmlcov/index.html
```

---

## 📚 Documentación

### README.md

```markdown
# Descripción clara y concisa

## Requisitos
- Lista de dependencias
- Versiones

## Instalación
Paso a paso

## Uso
Ejemplos con código

## API
Si corresponde

## Troubleshooting
Problemas comunes
```

### Docstrings

Todos los módulos, funciones y clases deben tener docstrings:

```python
"""
Module-level docstring.

This module contains functions for extracting metadata from NCBI XML files.
"""

def my_function():
    """One-line summary.
    
    Longer description if needed.
    """
    pass
```

### Comentarios Útiles

```python
# ✅ Útil: explica por qué, no qué
# We use a set here because O(1) lookup is critical for performance
# with 1.3M files
bioprojects_seen = set()

# ❌ Inútil: es obvio del código
for item in items:  # Iterate through items
    process(item)
```

---

## 📋 Proceso de PR

### Template de PR

```markdown
## Descripción
Breve descripción de qué cambia

## Tipo de cambio
- [ ] Bug fix
- [ ] New feature
- [ ] Enhancement
- [ ] Documentation

## Fase del Proyecto (si aplica)
- [ ] Fase 1: Extract metadata
- [ ] Fase 2: Query optimizer
- [ ] Fase 3: Optimization

## Tareas Completadas
- [x] Implementé feature X
- [x] Agregué tests (70%+ coverage)
- [x] Actualicé documentación
- [ ] Validé en local

## Testing
Cómo testear estos cambios:
```bash
pytest tests/test_my_feature.py -v
```

## Checklist
- [ ] Mi código sigue PEP 8
- [ ] He agregado tests para mi cambio
- [ ] Mis tests pasan localmente
- [ ] He actualizado la documentación
- [ ] No hay conflictos con main
```

### Revisión y Feedback

Los maintainers pueden solicitar cambios. Es normal, completamente esperado:

```bash
# 1. Hacer los cambios solicitados
git add .
git commit -m "Address PR feedback: ..."

# 2. Push (el PR se actualizará automáticamente)
git push origin feature/nombre

# 3. Esperar nueva revisión
```

---

## 🎯 Tareas Disponibles

### FASE 1: Extracción de Knowledge Base

#### T1.1 Parser XML robusto
**Dificultad:** ⭐⭐ (Fácil-Medio)  
**Requisitos:** Python, XML parsing  
**Descripción:**
- Crear function que parsee experiment.xml, sample.xml y run.xml
- Manejar errores (XML corrompidos)
- Logging
- Tests para 100+ archivos XML

**Archivos a crear:**
```
scripts/extract_metadata.py  # Main parser
tests/test_extract_metadata.py
tests/fixtures/  # XML samples para testing
```

**Entregable:**
```python
def extract_organism_data(xml_file: str) -> dict
def extract_strategy_data(xml_file: str) -> dict
def extract_attributes_data(xml_file: str) -> dict
```

---

#### T1.2 Proceso paralelo de 100k archivos
**Dificultad:** ⭐⭐⭐ (Medio)  
**Requisitos:** Python multiprocessing/asyncio  
**Descripción:**
- Script que procese primeros 100k archivos en paralelo
- Recolectar índices únicos
- Contar frecuencias
- Exportar JSON

**Archivo a crear:**
```
scripts/build_kb.py
```

**Entregable:**
```bash
python scripts/build_kb.py --max-files 100000
# Génera:
# data/kb/organisms_index.json
# data/kb/strategies_index.json
# data/kb/treatments_index.json
```

---

#### T1.3 Normalización de sinónimos
**Dificultad:** ⭐⭐⭐⭐ (Difícil)  
**Requisitos:** Bioinformatica, regex  
**Descripción:**
- Detectar variantes de organismos (ej: "arabidopsis", "Arabidopsis thaliana", "ath")
- Normalizar names
- Crear mappings
- Validar contra NCBI Taxonomy

**Archivos:**
```
scripts/utils/normalize.py
scripts/utils/taxonomy_mapping.py
tests/test_normalize.py
```

---

### FASE 2: Query Optimizer

#### T2.1 Query Builder determinístico
**Dificultad:** ⭐⭐⭐ (Medio)  
**Requisitos:** NCBI query syntax, lógica booleana  
**Descripción:**
- Implementar QueryOptimizer class
- Mapear conceptos usuario → campos NCBI
- Construir queries booleanas válidas
- Manejar ambigüedades

**Archivo a crear:**
```
scripts/query_optimizer.py
```

---

#### T2.2 Test suite exhaustivo
**Dificultad:** ⭐⭐ (Fácil-Medio)  
**Requisitos:** Pytest, casos SRA reales  
**Descripción:**
- 50+ test cases
- Validar contra baseline LLM
- Casos edge: typos, ambigüedades, términos raros

**Archivo a crear:**
```
tests/test_query_optimizer.py
tests/cases/test_queries.json  # Casos de prueba
```

---

### FASE 3: Refinement & Optimization

#### T3.1 Escalabilidad a 1.3M archivos
**Dificultad:** ⭐⭐⭐⭐ (Difícil)  
**Requisitos:** Performance optimization, memoria  
**Descripción:**
- Procesar todos los 1.3M archivos
- Optimizar performance
- Implementar caché
- Generar KB final completa

---

#### T3.2 Evaluación de precisión
**Dificultad:** ⭐⭐⭐ (Medio)  
**Requisitos:** Estadística, análisis  
**Descripción:**
- Comparar queries nuevas vs antiguas
- Medir precisión/recall
- Documentar resultados
- Generar reporte

---

## 📞 Preguntas?

- **GitHub Issues:** [Issues](https://github.com/lucianofrancoo/Pyner_PGRLAB/issues)
- **Discussions:** [Discussions](https://github.com/lucianofrancoo/Pyner_PGRLAB/discussions)
- **Email:** lucianofranco.a@gmail.com

---

## 🎓 Recursos Útiles

- [Git Documentation](https://git-scm.com/doc)
- [Python PEP 8](https://www.python.org/dev/peps/pep-0008/)
- [Pytest Documentation](https://docs.pytest.org/)
- [NCBI Query Syntax](https://www.ncbi.nlm.nih.gov/books/NBK3827/)
- [XML Parsing Python](https://docs.python.org/3/library/xml.etree.elementtree.html)

---

## 📜 Código de Conducta

Por favor sé respetuoso y constructivo. La comunidad científica se basa en la colaboración.

---

**¡Gracias por contribuir! 🚀**
