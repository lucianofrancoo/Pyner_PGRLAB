#!/bin/bash
#
# PYNER COMPLETE WORKFLOW WITH ANALYSIS
# ======================================
# Complete workflow: Natural Language → Boolean Query → BioProject/PubMed → AI Analysis → Classified Table
#
# Features:
# 1. Natural language to boolean query (Phase 3)
# 2. Database selection (BioProject or PubMed)
# 3. Fetch papers/projects from NCBI
# 4. AI-powered analysis with Ollama (relevance scoring + structured extraction)
# 5. Final output: Classified table with organisms, tissues, conditions, strategies
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
    
    # Check Ollama (for analysis)
    if ! python3 -c "import sys; sys.path.insert(0, 'Data_Analyzer'); from ollama_client import OllamaClient; exit(0 if OllamaClient().is_available() else 1)" &> /dev/null; then
        echo -e "${YELLOW}WARNING: Ollama not available${NC}"
        echo -e "  Analysis step will be skipped"
        echo -e "  To enable: ${YELLOW}ollama serve${NC} (in another terminal)"
        echo -e "  And: ${YELLOW}ollama pull qwen2.5:14b${NC}"
        OLLAMA_AVAILABLE=false
    else
        OLLAMA_AVAILABLE=true
    fi
    
    # Check required files
    local required_files=(
        "Query_generator/phases/phase3/api/main.py"
        "Fetcher_NCBI/pubmed_boolean_search.py"
        "Fetcher_NCBI/boolean_fetcher_integrated.py"
        "Fetcher_NCBI/config.py"
        "Data_Analyzer/paper_analyzer.py"
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
printf "${BLUE}    PYNER - COMPLETE WORKFLOW WITH AI ANALYSIS${NC}\n"
printf "${BLUE}================================================================${NC}\n\n"

# Step 1: Get user input
printf "${GREEN}[1/6] Enter your query in natural language:${NC}\n"
read -r -p "> " USER_INPUT

if [ -z "$USER_INPUT" ]; then
    echo -e "${RED}ERROR: Empty input. Aborting.${NC}"
    exit 1
fi

# Store original query for analysis
ORIGINAL_QUERY="$USER_INPUT"

# Step 2: Ask user what they want to do
printf "\n${GREEN}[2/6] What do you want to do?${NC}\n"
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
printf "\n${GREEN}[3/6] Generating boolean query with AI...${NC}\n\n"

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
    printf "\n${GREEN}[4/6] How many publications do you want to retrieve?${NC}\n"
    printf "  ${YELLOW}(Press Enter to get all available)${NC}\n"
    read -r -p "Maximum number [default: unlimited]: " MAX_RESULTS
    MAX_RESULTS=${MAX_RESULTS:-0}
else
    printf "\n${GREEN}[4/6] How many BioProjects do you want to process?${NC}\n"
    printf "  ${YELLOW}(Press Enter to get all available)${NC}\n"
    read -r -p "Maximum number [default: unlimited]: " MAX_RESULTS
    MAX_RESULTS=${MAX_RESULTS:-0}
fi

# Convert 0 to a large number (practical "unlimited")
if [ "$MAX_RESULTS" = "0" ]; then
    MAX_RESULTS=10000
fi

# Step 5b: Execute fetcher based on user selection
cd "$ROOT_DIR/Fetcher_NCBI"

if [ "$DATABASE" = "pubmed" ]; then
    # Direct PubMed search (FAST)
    OUTPUT_CSV="$ROOT_DIR/pubmed_results_${TIMESTAMP}.csv"
    OUTPUT_JSON="$ROOT_DIR/pubmed_results_${TIMESTAMP}.json"
    
    printf "\n${GREEN}[5/6] Executing direct PubMed search...${NC}\n\n"
    
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
    
    printf "\n${GREEN}[5/6] Executing integrated workflow BioProject → SRA → PubMed...${NC}\n\n"
    
    python3 boolean_fetcher_integrated.py "$QUERY" \
        --max "$MAX_RESULTS" \
        --output-csv "$OUTPUT_CSV" \
        --output-json "$OUTPUT_JSON"
    
    printf "\n${GREEN}✓ WORKFLOW COMPLETED${NC}\n"
    printf "Results saved in:\n"
    printf "  📄 CSV:  %s\n" "$OUTPUT_CSV"
    printf "  📄 JSON: %s\n" "$OUTPUT_JSON"
fi

# Step 6: AI Analysis with Data_Analyzer (only for PubMed for now)
if [ "$DATABASE" = "pubmed" ] && [ "$OLLAMA_AVAILABLE" = true ]; then
    printf "\n${BLUE}================================================================${NC}\n"
    printf "${GREEN}[6/6] Analyzing papers with AI (Ollama)...${NC}\n"
    printf "${BLUE}================================================================${NC}\n\n"
    
    cd "$ROOT_DIR/Data_Analyzer"
    
    ANALYZED_CSV="$ROOT_DIR/Data_Analyzer/output/classified_papers_${TIMESTAMP}.csv"
    
    # Run analysis (pass the original user query context for better relevance scoring)
    printf "Original query for context: ${YELLOW}$ORIGINAL_QUERY${NC}\n\n"
    
    python3 paper_analyzer.py "$OUTPUT_JSON" "$ANALYZED_CSV"
    
    ANALYSIS_EXIT=$?
    
    if [ $ANALYSIS_EXIT -eq 0 ] && [ -f "$ANALYZED_CSV" ]; then
        printf "\n${BLUE}================================================================${NC}\n"
        printf "${GREEN}✓ COMPLETE WORKFLOW FINISHED${NC}\n"
        printf "${BLUE}================================================================${NC}\n\n"
        printf "${GREEN}Final Results:${NC}\n"
        printf "  📊 Classified Table: ${YELLOW}%s${NC}\n" "$ANALYZED_CSV"
        printf "  📄 Raw CSV:         %s\n" "$OUTPUT_CSV"
        printf "  📄 Raw JSON:        %s\n" "$OUTPUT_JSON"
        
        # Show preview of classified results
        printf "\n${YELLOW}Preview of classified papers (first 2):${NC}\n"
        printf "${BLUE}────────────────────────────────────────────────────────────────${NC}\n"
        head -3 "$ANALYZED_CSV" | cut -c1-150 | while IFS= read -r line; do
            printf "%s...\n" "$line"
        done
        printf "${BLUE}────────────────────────────────────────────────────────────────${NC}\n"
        
        # Statistics
        total_papers=$(( $(wc -l < "$ANALYZED_CSV") - 1 ))
        relevant_papers=$(awk -F',' 'NR>1 && $4=="Yes" {count++} END {print count+0}' "$ANALYZED_CSV")
        
        printf "\n${GREEN}Analysis Summary:${NC}\n"
        printf "  Total papers:    ${BLUE}%d${NC}\n" "$total_papers"
        printf "  Relevant papers: ${GREEN}%d${NC} (%.1f%%)\n" "$relevant_papers" "$(echo "scale=1; $relevant_papers * 100 / $total_papers" | bc 2>/dev/null || echo "N/A")"
        printf "  Non-relevant:    ${YELLOW}%d${NC}\n" "$((total_papers - relevant_papers))"
        
    else
        printf "\n${RED}WARNING: Analysis step failed${NC}\n"
        printf "Check log: Data_Analyzer/logs/analyzer_${TIMESTAMP}.log\n"
        printf "Raw results still available in: $OUTPUT_CSV\n"
    fi
    
elif [ "$DATABASE" = "bioproject" ]; then
    printf "\n${YELLOW}Note: AI analysis is currently only supported for PubMed searches${NC}\n"
    printf "${YELLOW}BioProject analysis will be added in future updates${NC}\n"
    
elif [ "$OLLAMA_AVAILABLE" = false ]; then
    printf "\n${YELLOW}Note: Ollama not available - skipping analysis step${NC}\n"
    printf "To enable analysis, start Ollama: ${YELLOW}ollama serve${NC}\n"
fi

printf "\n${BLUE}================================================================${NC}\n\n"
