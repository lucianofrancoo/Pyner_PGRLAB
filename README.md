# Pyner_LA

Generador de términos de búsqueda optimizados para NCBI SRA (Sequence Read Archive) usando modelos de lenguaje locales.

## Descripción

Pyner_LA es una herramienta de bioinformática que utiliza LLMs locales (Ollama) para generar consultas booleanas optimizadas para búsquedas en NCBI SRA. La herramienta convierte palabras clave en términos de búsqueda compatibles con E-utilities de NCBI, utilizando los campos correctos específicos de SRA como `[Organism]`, `[Title]`, `[Strategy]`, etc.

## Requisitos

- Python 3.x
- [Ollama](https://ollama.ai/) instalado localmente
- Modelo `qwen2.5:14b` descargado en Ollama

## Instalación

1. Clonar el repositorio:
```bash
git clone <url-del-repositorio>
cd Pyner_LA
```

2. Asegurarse de tener Ollama instalado y el modelo descargado:
```bash
ollama pull qwen2.5:14b
```

## Uso

### Pyner_search_v0.1.py

Genera términos de búsqueda a partir de palabras clave:

```bash
python3 Pyner_search_v0.1.py solanum lycopersicum nitrogen
```

El script generará:
- Una consulta en lenguaje natural
- Una consulta optimizada para NCBI SRA ESearch

## Archivos del proyecto

- `Pyner_search_v0.1.py` - Generador de términos de búsqueda para SRA
- `Pyner_v0.1.py` - Versión inicial del sistema
- `Pyner_v0.2.py` - Versión mejorada del sistema
- `Pyner_SRA_arabidopsis_drought_unique.csv` - Datos de ejemplo

## Contribuir

Las contribuciones son bienvenidas. Por favor:
1. Haz un fork del proyecto
2. Crea una rama para tu feature (`git checkout -b feature/nueva-funcionalidad`)
3. Haz commit de tus cambios (`git commit -am 'Añadir nueva funcionalidad'`)
4. Push a la rama (`git push origin feature/nueva-funcionalidad`)
5. Crea un Pull Request

## Licencia

Por definir

## Autores

- Luis A. Humada (lahumada)
