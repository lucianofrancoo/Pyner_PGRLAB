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

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# ================================================================
# ENVIRONMENT VALIDATION
# ================================================================

validate_environment() {
    local errors=0
    
    # Check Python 3
    if ! command -v python3 &> /dev/null; then
        echo -e "${RED}ERROR: python3 not found${NC}"
        echo -e "  Please install Python 3: ${YELLOW}sudo apt install python3${NC}"
        ((errors++))
    fi
    
    # Check BioPython
    if ! python3 -c "import Bio" &> /dev/null; then
        echo -e "${RED}ERROR: BioPython not installed${NC}"
        echo -e "  Install with: ${YELLOW}pip install biopython${NC}"
        echo -e "  Or: ${YELLOW}pip install -r Fetcher_NCBI/requirements.txt${NC}"
        ((errors++))
    fi
    
    # Check required files
    local required_files=(
        "Query_generator/phases/phase3/api/main.py"
        "Fetcher_NCBI/pubmed_boolean_search.py"
        "Fetcher_NCBI/boolean_fetcher_integrated.py"
        "Fetcher_NCBI/config.py"
    )
    
    for file in "${required_files[@]}"; do
        if [ ! -f "$ROOT_DIR/$file" ]; then
            echo -e "${RED}ERROR: Missing required file: $file${NC}"
            ((errors++))
        fi
    done
    
    if [ $errors -gt 0 ]; then
        echo -e "\n${RED}❌ Environment validation failed with $errors error(s)${NC}"
        echo -e "${YELLOW}Please fix the issues above and try again.${NC}\n"
        exit 1
    fi
    
    return 0
}

# Run validation
validate_environment

# Now enable strict error handling after validation
set -eo pipefail

printf "\n${BLUE}================================================================${NC}\n"
printf "${BLUE}         PYNER - INTEGRATED BOOLEAN FETCHER${NC}\n"
printf "${BLUE}================================================================${NC}\n\n"

# Step 1: Get user input
printf "${GREEN}[1/5] Enter your query in natural language:${NC}\n"
read -r -p "> " USER_INPUT

if [ -z "$USER_INPUT" ]; then
    echo -e "${RED}ERROR: Empty input. Aborting.${NC}"
    exit 1
fi

# Step 2: Ask user what they want to do
printf "\n${GREEN}[2/5] What do you want to do?${NC}\n"
printf "  ${YELLOW}1)${NC} Search Papers (PubMed) ${BLUE}→${NC} Fast, list of publications with PMID\n"
printf "  ${YELLOW}2)${NC} Search BioProjects ${BLUE}→${NC} Slow, searches omics projects and associated papers via cascade BioProject → BioSample → BioExperiment\n"
read -r -p "Selection [1-2]: " DB_CHOICE

case $DB_CHOICE in
    1)
        DATABASE="pubmed"
        printf "  ✓ Selected: ${BLUE}Papers Search (PubMed)${NC}\n"
        ;;
    2)
        DATABASE="bioproject"
        printf "  ✓ Selected: ${BLUE}BioProjects Search with cascade${NC}\n"
        ;;
    *)
        echo -e "${RED}ERROR: Invalid option. Using Papers search by default.${NC}"
        DATABASE="pubmed"
        ;;
esac

# Step 3: Generate boolean query with Phase 3
printf "\n${GREEN}[3/5] Generating boolean query with AI...${NC}\n\n"

cd "$ROOT_DIR/Query_generator/phases/phase3"
GEN_OUTPUT=$(python3 api/main.py --quick "$USER_INPUT" 2>/dev/null)

echo "$GEN_OUTPUT"

# Extract the query
QUERY=$(printf "%s\n" "$GEN_OUTPUT" | awk '/NCBI Query:/{getline; gsub(/^ +/, ""); print; exit}')

if [ -z "$QUERY" ]; then
    echo -e "${RED}ERROR: Could not extract query.${NC}"
    exit 1
fi

# Step 4: Confirm query
printf "\n${YELLOW}Generated query: ${BLUE}$QUERY${NC}\n"
read -r -p "Continue with this query? [Y/n]: " CONFIRM
CONFIRM=$(echo "$CONFIRM" | tr '[:upper:]' '[:lower:]')

if [ "$CONFIRM" = "n" ] || [ "$CONFIRM" = "no" ]; then
    echo -e "${RED}Cancelled by user.${NC}"
    exit 0
fi

# Step 5a: Ask for max results
if [ "$DATABASE" = "pubmed" ]; then
    printf "\n${GREEN}[4/5] How many publications do you want to retrieve?${NC}\n"
    printf "  ${YELLOW}(Press Enter to get all available)${NC}\n"
    read -r -p "Maximum number [default: unlimited]: " MAX_RESULTS
    MAX_RESULTS=${MAX_RESULTS:-0}
else
    printf "\n${GREEN}[4/5] How many BioProjects do you want to process?${NC}\n"
    printf "  ${YELLOW}(Press Enter to get all available)${NC}\n"
    read -r -p "Maximum number [default: unlimited]: " MAX_RESULTS
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
    
    printf "\n${GREEN}[5/5] Executing direct PubMed search...${NC}\n\n"
    
    python3 pubmed_boolean_search.py "$QUERY" \
        --max "$MAX_RESULTS" \
        --output-csv "$OUTPUT_CSV" \
        --output-json "$OUTPUT_JSON"
    
    printf "\n${GREEN}✓ SEARCH COMPLETED${NC}\n"
    printf "Results saved in:\n"
    printf "  📄 CSV:  %s\n" "$OUTPUT_CSV"
    printf "  📄 JSON: %s\n" "$OUTPUT_JSON"

elif [ "$DATABASE" = "bioproject" ]; then
    # BioProject workflow with full integration (SLOW - with cascade)
    OUTPUT_CSV="$ROOT_DIR/results_${TIMESTAMP}.csv"
    OUTPUT_JSON="$ROOT_DIR/results_${TIMESTAMP}.json"
    
    printf "\n${GREEN}[5/5] Executing integrated workflow BioProject → SRA → PubMed...${NC}\n\n"
    
    python3 boolean_fetcher_integrated.py "$QUERY" \
        --max "$MAX_RESULTS" \
        --output-csv "$OUTPUT_CSV" \
        --output-json "$OUTPUT_JSON"
    
    printf "\n${GREEN}✓ WORKFLOW COMPLETE${NC}\n"
    printf "Results saved in:\n"
    printf "  📄 CSV:  %s\n" "$OUTPUT_CSV"
    printf "  📄 JSON: %s\n" "$OUTPUT_JSON"
fi

printf "\n${BLUE}================================================================${NC}\n"
