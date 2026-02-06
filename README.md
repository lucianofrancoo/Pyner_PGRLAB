# Pyner_PGRLAB — Project Overview (minimal)

Proyecto Pyner: generador de queries semánticos sobre datos NCBI (MVP sesgado).

Estructura importante:

- `Query_generator/` — Código y documentación para el Query Generator MVP.
  - `Query_generator/docs/` — Documentación y roadmaps específicos.
  - `Query_generator/phases/` — Copias de `phase1/`, `phase2/`, `phase3/` usadas por el proyecto.
  - `Query_generator/run_tests.sh` — Script rápido para validar el demo.

- `scripts_iniciales_beta/` — Scripts antiguos iniciales (archivado).
- `archive_old/` — Requisitos, CSVs y otros archivos auxiliares (archivado).

Estado de avance (resumen):
- Phase 1: Procesamiento de 500K archivos — COMPLETO (KB generado)
- Phase 2: Creación de FAISS index y queries — COMPLETO
- Phase 3: API (FastAPI) + Ollama integration — COMPLETO (MVP)

Nota: Los originales `phase1/`, `phase2/`, `phase3/` han sido movidos/copied a
`Query_generator/phases/` para mantener la raíz ordenada. Usa los archivos dentro
de `Query_generator/` para demos y pruebas.

Para probar rápidamente:

```bash
python3 Query_generator/quick_test.py
# o
python3 Query_generator/test_queries.py --interactive
```

Contacto: mantén este README mínimo como índice; la documentación completa está
en `Query_generator/docs/`.
