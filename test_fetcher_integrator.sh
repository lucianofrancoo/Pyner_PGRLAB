#!/bin/bash
#
# PYNER INTEGRATED WORKFLOW
# ==========================
# Complete workflow: Natural Language → Boolean Query → BioProject/PubMed → SRA → Publications → CSV
#
# Features:
# 1. Natural language to boolean query (Phase 3)
# 2. Database selection (BioProject or PubMed)
# 3. For BioProject: Full fetch of SRA + Cascade PubMed linking
# 4. CSV/JSON output with publication info
#

set -e

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

printf "\n${BLUE}================================================================${NC}\n"
printf "${BLUE}         PYNER - INTEGRATED BOOLEAN FETCHER${NC}\n"
printf "${BLUE}================================================================${NC}\n\n"

# Step 1: Get user input
printf "${GREEN}[1/5] Ingresa tu consulta en lenguaje natural:${NC}\n"
read -r -p "> " USER_INPUT

if [ -z "$USER_INPUT" ]; then
    echo -e "${RED}ERROR: Input vacío. Abortando.${NC}"
    exit 1
fi

# Step 2: Ask user what they want to do
printf "\n${GREEN}[2/5] ¿Qué quieres hacer?${NC}\n"
printf "  ${YELLOW}1)${NC} Buscar Papers (PubMed) ${BLUE}→${NC} Rápido, lista de publicaciones con PMID\n"
printf "  ${YELLOW}2)${NC} Buscar BioProjects ${BLUE}→${NC} Lento, busca proyectos ómicos y papers asociados por cascada BioProject --> BioSample --> BioExperiment\n"
read -r -p "Selección [1-2]: " DB_CHOICE

case $DB_CHOICE in
    1)
        DATABASE="pubmed"
        printf "  ✓ Seleccionado: ${BLUE}Búsqueda de Papers (PubMed)${NC}\n"
        ;;
    2)
        DATABASE="bioproject"
        printf "  ✓ Seleccionado: ${BLUE}Búsqueda de BioProjects con cascada${NC}\n"
        ;;
    *)
        echo -e "${RED}ERROR: Opción inválida. Usando búsqueda de Papers por defecto.${NC}"
        DATABASE="pubmed"
        ;;
esac

# Step 3: Generate boolean query with Phase 3
printf "\n${GREEN}[3/5] Generando query booleano con IA...${NC}\n\n"

cd "$ROOT_DIR/Query_generator/phases/phase3"
GEN_OUTPUT=$(python api/main.py --quick "$USER_INPUT" 2>/dev/null)

echo "$GEN_OUTPUT"

# Extract the query
QUERY=$(printf "%s\n" "$GEN_OUTPUT" | awk '/NCBI Query:/{getline; gsub(/^ +/, ""); print; exit}')

if [ -z "$QUERY" ]; then
    echo -e "${RED}ERROR: No se pudo extraer el query.${NC}"
    exit 1
fi

# Step 4: Confirm query
printf "\n${YELLOW}Query generado: ${BLUE}$QUERY${NC}\n"
read -r -p "¿Continuar con este query? [S/n]: " CONFIRM
CONFIRM=$(echo "$CONFIRM" | tr '[:upper:]' '[:lower:]')

if [ "$CONFIRM" = "n" ] || [ "$CONFIRM" = "no" ]; then
    echo -e "${RED}Cancelado por el usuario.${NC}"
    exit 0
fi

# Step 5a: Ask for max results
if [ "$DATABASE" = "pubmed" ]; then
    printf "\n${GREEN}[4/5] ¿Cuántas publicaciones quieres recuperar?${NC}\n"
    printf "  ${YELLOW}(Presiona Enter para obtener todas las disponibles)${NC}\n"
    read -r -p "Número máximo [default: sin límite]: " MAX_RESULTS
    MAX_RESULTS=${MAX_RESULTS:-0}
else
    printf "\n${GREEN}[4/5] ¿Cuántos BioProjects quieres procesar?${NC}\n"
    printf "  ${YELLOW}(Presiona Enter para obtener todos los disponibles)${NC}\n"
    read -r -p "Número máximo [default: sin límite]: " MAX_RESULTS
    MAX_RESULTS=${MAX_RESULTS:-0}
fi

# Convert 0 to a large number (practical "unlimited")
if [ "$MAX_RESULTS" = "0" ]; then
    MAX_RESULTS=10000
fi

# Step 5b: Execute based on user selection
cd "$ROOT_DIR/Fetcher_NCBI"

if [ "$DATABASE" = "pubmed" ]; then
    # Direct PubMed search (FAST)
    OUTPUT_CSV="$ROOT_DIR/pubmed_results_${TIMESTAMP}.csv"
    OUTPUT_JSON="$ROOT_DIR/pubmed_results_${TIMESTAMP}.json"
    
    printf "\n${GREEN}[5/5] Ejecutando búsqueda directa en PubMed...${NC}\n\n"
    
    python pubmed_boolean_search.py "$QUERY" \
        --max "$MAX_RESULTS" \
        --output-csv "$OUTPUT_CSV" \
        --output-json "$OUTPUT_JSON"
    
    printf "\n${GREEN}✓ BÚSQUEDA COMPLETADA${NC}\n"
    printf "${BLUE}Resultados guardados en:${NC}\n"
    printf "  📄 CSV:  %s\n" "$OUTPUT_CSV"
    printf "  📄 JSON: %s\n" "$OUTPUT_JSON"
    
    # Show preview
    if [ -f "$OUTPUT_CSV" ]; then
        printf "\n${YELLOW}Vista previa (primeras 3 líneas):${NC}\n"
        head -3 "$OUTPUT_CSV" | cut -c1-100
        echo "..."
        
        # Statistics
        total=$(( $(wc -l < "$OUTPUT_CSV") - 1 ))
        with_doi=$(awk -F',' 'NR>1 && $5 != "NA" && $5 != "" {count++} END {print count+0}' "$OUTPUT_CSV")
        
        printf "\n${GREEN}Estadísticas:${NC}\n"
        printf "  Total publicaciones: ${BLUE}%d${NC}\n" "$total"
        printf "  Con DOI: ${GREEN}%d${NC}\n" "$with_doi"
        printf "  Sin DOI: ${YELLOW}%d${NC}\n" "$((total - with_doi))"
    fi

elif [ "$DATABASE" = "bioproject" ]; then
    # BioProject workflow with full integration (SLOW - with cascade)
    OUTPUT_CSV="$ROOT_DIR/results_${TIMESTAMP}.csv"
    OUTPUT_JSON="$ROOT_DIR/results_${TIMESTAMP}.json"
    
    printf "\n${GREEN}[5/5] Ejecutando workflow integrado BioProject → SRA → PubMed...${NC}\n\n"
    
    python boolean_fetcher_integrated.py "$QUERY" \
        --max "$MAX_RESULTS" \
        --output-csv "$OUTPUT_CSV" \
        --output-json "$OUTPUT_JSON"
    
    printf "\n${GREEN}✓ WORKFLOW COMPLETADO${NC}\n"
    printf "${BLUE}Resultados guardados en:${NC}\n"
    printf "  📄 CSV:  %s\n" "$OUTPUT_CSV"
    printf "  📄 JSON: %s\n" "$OUTPUT_JSON"
    
    # Show preview
    if [ -f "$OUTPUT_CSV" ]; then
        printf "\n${YELLOW}Vista previa (primeras 3 líneas):${NC}\n"
        head -3 "$OUTPUT_CSV" | cut -c1-120
        echo "..."
        
        # Statistics
        total=$(( $(wc -l < "$OUTPUT_CSV") - 1 ))
        with_papers=$(awk -F',' 'NR>1 && $10 > 0 {count++} END {print count+0}' "$OUTPUT_CSV")
        without_papers=$(( total - with_papers ))
        
        printf "\n${GREEN}Estadísticas:${NC}\n"
        printf "  Total BioProjects: ${BLUE}%d${NC}\n" "$total"
        printf "  Con publicaciones: ${GREEN}%d${NC}\n" "$with_papers"
        printf "  Sin publicaciones: ${YELLOW}%d${NC}\n" "$without_papers"
    fi
fi

printf "\n${BLUE}================================================================${NC}\n"
