#!/bin/bash
#
# PYNER MINER - Complete Mining Workflow
# ========================================
# Integrated pipeline: Natural Language → Boolean Query → Fetch → Analysis (Optional)
#
# Features:
# 1. Lite Mode: Fast fetch only (basic CSV/JSON with PMID, title, abstract)
# 2. Pro Mode: Full fetch + AI analysis (53 metadata columns with comprehensive experimental details)
# 3. Natural language to boolean query conversion (Phase 3)
# 4. Database selection (BioProject or PubMed)
# 5. For BioProject: Full cascade linking (BioProject → SRA → PubMed)
#

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="$ROOT_DIR/output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
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
    if ! python3 -c "import Bio" &> /dev/null 2>&1; then
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
    
    if [ ! -d "$ROOT_DIR/Data_visualization" ]; then
        echo -e "${RED}ERROR: Missing Data_visualization directory${NC}"
        ((errors++))
    fi
    
    for file in "${required_files[@]}"; do
        if [ ! -f "$ROOT_DIR/$file" ]; then
            echo -e "${RED}ERROR: Missing required file: $file${NC}"
            ((errors++))
        fi
    done
    
    # If Pro mode will be selected, check Data Analyzer
    if [ $errors -eq 0 ]; then
        if [ ! -f "$ROOT_DIR/Data_Analyzer/paper_analyzer.py" ]; then
            echo -e "${YELLOW}WARNING: Data Analyzer not found (Pro mode will not be available)${NC}"
        fi
    fi
    
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

printf "\n${CYAN}╔═══════════════════════════════════════════════════════════════╗${NC}\n"
printf "${CYAN}║                                                               ║${NC}\n"
printf "${CYAN}║                    ${BLUE}🔬 PYNER MINER${CYAN}                         ║${NC}\n"
printf "${CYAN}║                                                               ║${NC}\n"
printf "${CYAN}║         ${YELLOW}Intelligent Literature & Data Mining Pipeline${CYAN}        ║${NC}\n"
printf "${CYAN}║                                                               ║${NC}\n"
printf "${CYAN}╚═══════════════════════════════════════════════════════════════╝${NC}\n\n"

# ================================================================
# STEP 1: ANALYSIS MODE SELECTION
# ================================================================

printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
printf "${GREEN}[1/6] Select analysis mode${NC}\n"
printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"

printf "  ${YELLOW}1)${NC} ${BLUE}Lite${NC} - Fast fetch only\n"
printf "     ${CYAN}→${NC} Basic TSV/JSON output (PMID, title, abstract, DOI)\n"
printf "     ${CYAN}→${NC} Quick results for literature review\n\n"

printf "  ${YELLOW}2)${NC} ${BLUE}Pro${NC}  - Full analysis with AI extraction\n"
printf "     ${CYAN}→${NC} 53 metadata columns with comprehensive experimental details\n"
printf "     ${CYAN}→${NC} AI-powered extraction: organisms, conditions, molecules, time courses, etc.\n"
printf "     ${CYAN}→${NC} Relevance scoring and quality metrics\n"
printf "     ${CYAN}→${NC} Requires Ollama with qwen3.5:9b model\n\n"

read -r -p "Selection [1-2]: " ANALYSIS_MODE

case $ANALYSIS_MODE in
    1)
        MODE="lite"
        printf "\n  ✓ Selected: ${BLUE}Lite mode${NC} (fetch only)\n"
        ;;
    2)
        MODE="pro"
        printf "\n  ✓ Selected: ${BLUE}Pro mode${NC} (fetch + full analysis)\n"
        
        # Validate Data Analyzer availability
        if [ ! -f "$ROOT_DIR/Data_Analyzer/paper_analyzer.py" ]; then
            echo -e "${RED}ERROR: Data Analyzer not found${NC}"
            echo -e "Pro mode requires Data_Analyzer/paper_analyzer.py"
            exit 1
        fi
        
        # Check Ollama connection for Pro mode
        if ! python3 -c "import sys; sys.path.insert(0, '$ROOT_DIR/Data_Analyzer'); from config import OLLAMA_BASE_URL; import requests; requests.get(OLLAMA_BASE_URL, timeout=2)" 2>/dev/null; then
            echo -e "${YELLOW}WARNING: Cannot connect to Ollama${NC}"
            read -r -p "Continue anyway? [y/N]: " CONTINUE
            CONTINUE=$(echo "$CONTINUE" | tr '[:upper:]' '[:lower:]')
            if [ "$CONTINUE" != "y" ] && [ "$CONTINUE" != "yes" ]; then
                echo -e "${RED}Cancelled by user.${NC}"
                exit 0
            fi
        fi
        ;;
    *)
        echo -e "${RED}ERROR: Invalid option. Using Lite mode by default.${NC}"
        MODE="lite"
        ;;
esac

# ================================================================
# STEP 2: NATURAL LANGUAGE QUERY INPUT
# ================================================================

printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
printf "${GREEN}[2/6] Enter your research question${NC}\n"
printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"

printf "${CYAN}Examples:${NC}\n"
printf "  - \"Tomato drought stress RNA-Seq\"\n"
printf "  - \"Arabidopsis phosphate starvation transcriptome\"\n"
printf "  - \"Maize cold stress gene expression\"\n\n"

read -r -p "Your query: " USER_INPUT

if [ -z "$USER_INPUT" ]; then
    echo -e "${RED}ERROR: Empty input. Aborting.${NC}"
    exit 1
fi

printf "\n  ✓ Query: ${BLUE}$USER_INPUT${NC}\n"

# ================================================================
# STEP 3: DATABASE SELECTION
# ================================================================

printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
printf "${GREEN}[3/6] Select database to search${NC}\n"
printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"

printf "  ${YELLOW}1)${NC} ${BLUE}PubMed${NC} - Literature search (title/abstract)\n"
printf "     ${CYAN}→${NC} Fast: Search titles and abstracts\n"
printf "     ${CYAN}→${NC} Returns: Papers with PMID, abstract, metadata\n\n"

printf "  ${YELLOW}2)${NC} ${BLUE}PMC${NC} - Full-text literature search\n"
printf "     ${CYAN}→${NC} Searches entire article body (more results)\n"
printf "     ${CYAN}→${NC} Returns: Papers with PMID, abstract, metadata\n\n"

printf "  ${YELLOW}3)${NC} ${BLUE}BioProject${NC} - Omics data search\n"
printf "     ${CYAN}→${NC} Slow: Full cascade linking (BioProject → SRA → PubMed)\n"
printf "     ${CYAN}→${NC} Returns: Experiments + associated publications\n\n"

read -r -p "Selection [1-3]: " DB_CHOICE

case $DB_CHOICE in
    1)
        DATABASE="pubmed"
        SEARCH_DB="pubmed"
        printf "\n  ✓ Selected: ${BLUE}PubMed${NC} (title/abstract search)\n"
        ;;
    2)
        DATABASE="pubmed"
        SEARCH_DB="pmc"
        printf "\n  ✓ Selected: ${BLUE}PMC${NC} (full-text search)\n"
        ;;
    3)
        DATABASE="bioproject"
        SEARCH_DB="pubmed"
        printf "\n  ✓ Selected: ${BLUE}BioProject${NC} (omics data search with cascade)\n"
        ;;
    *)
        echo -e "${RED}ERROR: Invalid option. Using PubMed by default.${NC}"
        DATABASE="pubmed"
        SEARCH_DB="pubmed"
        ;;
esac

# ================================================================
# STEP 4: BOOLEAN QUERY GENERATION
# ================================================================

printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
printf "${GREEN}[4/6] Generating boolean query with AI...${NC}\n"
printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"

cd "$ROOT_DIR/Query_generator/phases/phase3"

# Create a temporary file to capture the output while allowing interaction
TEMP_GEN_OUTPUT=$(mktemp)

# Export SEARCH_DB so the interactive CLI can adapt field tags for PMC
export SEARCH_DB

# Run main.py in interactive mode (default) and tee output to file
# We don't use --quick to allow the new refinement loop
python3 api/main.py "$USER_INPUT" 2>/dev/null | tee "$TEMP_GEN_OUTPUT"
GEN_EXIT_CODE=${PIPESTATUS[0]}

if [ $GEN_EXIT_CODE -ne 0 ]; then
    echo -e "${RED}Generation cancelled or failed.${NC}"
    rm "$TEMP_GEN_OUTPUT"
    exit 1
fi

# Extract the LAST occurrence of the query from the log
QUERY=$(grep "NCBI Query:" "$TEMP_GEN_OUTPUT" -A 1 | tail -n 1 | sed 's/^[[:space:]]*//')
rm "$TEMP_GEN_OUTPUT"

if [ -z "$QUERY" ]; then
    echo -e "${RED}ERROR: Could not extract query.${NC}"
    exit 1
fi

# If PMC selected, ensure field tags are converted (safety net)
if [ "$SEARCH_DB" = "pmc" ]; then
    QUERY=$(echo "$QUERY" | sed 's/\[Organism\]/[all]/g' | sed 's/\[All Fields\]/[all]/g')
fi

# ================================================================
# STEP 5: MAX RESULTS CONFIGURATION
# ================================================================

printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
printf "${GREEN}[5/6] Configure result limits${NC}\n"
printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"

if [ "$DATABASE" = "pubmed" ]; then
    printf "How many publications to retrieve?\n"
    printf "${YELLOW}(Press Enter for unlimited)${NC}\n\n"
    read -r -p "Maximum number [default: unlimited]: " MAX_RESULTS
    MAX_RESULTS=${MAX_RESULTS:-0}
else
    printf "How many BioProjects to process?\n"
    printf "${YELLOW}(Press Enter for unlimited)${NC}\n\n"
    read -r -p "Maximum number [default: unlimited]: " MAX_RESULTS
    MAX_RESULTS=${MAX_RESULTS:-0}
fi

# Convert 0 to a large number (practical "unlimited")
if [ "$MAX_RESULTS" = "0" ]; then
    MAX_RESULTS=10000
    printf "\n  ✓ Limit: ${BLUE}Unlimited${NC} (max 10,000)\n"
else
    printf "\n  ✓ Limit: ${BLUE}%s${NC}\n" "$MAX_RESULTS"
fi

# ================================================================
# STEP 6: EXECUTION
# ================================================================

printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
printf "${GREEN}[6/6] Executing mining pipeline...${NC}\n"
printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"

cd "$ROOT_DIR/Fetcher_NCBI"

if [ "$DATABASE" = "pubmed" ]; then
    # ============================================================
    # PubMed Direct Search
    # ============================================================
    
    FETCH_TSV="$OUTPUT_DIR/pubmed_fetch_${TIMESTAMP}.tsv"
    FETCH_JSON="$OUTPUT_DIR/pubmed_fetch_${TIMESTAMP}.json"
    
    if [ "$SEARCH_DB" = "pmc" ]; then
        printf "${CYAN}[Fetcher]${NC} Searching PMC (full-text)...\n\n"
    else
        printf "${CYAN}[Fetcher]${NC} Searching PubMed...\n\n"
    fi
    
    python3 pubmed_boolean_search.py "$QUERY" \
        --max "$MAX_RESULTS" \
        --db "$SEARCH_DB" \
        --output-tsv "$FETCH_TSV" \
        --output-json "$FETCH_JSON"
    
    printf "\n${GREEN}✓ Fetch completed${NC}\n"
    printf "  📄 TSV:  %s\n" "$FETCH_TSV"
    printf "  📄 JSON: %s\n" "$FETCH_JSON"
    
    # If Pro mode, run Data Analyzer
    if [ "$MODE" = "pro" ]; then
        printf "\n${CYAN}[Data Analyzer]${NC} Running AI analysis with Ollama...\n\n"
        
        cd "$ROOT_DIR/Data_Analyzer"
        ANALYSIS_TSV="$OUTPUT_DIR/classified_papers_${TIMESTAMP}.tsv"
        
        python3 paper_analyzer.py "$FETCH_JSON" "$ANALYSIS_TSV"
        
        printf "\n${GREEN}✓ Analysis completed${NC}\n"
        printf "  📊 Classified TSV: %s\n" "$ANALYSIS_TSV"
    fi

elif [ "$DATABASE" = "bioproject" ]; then
    # ============================================================
    # BioProject Cascade Search
    # ============================================================
    
    FETCH_TSV="$OUTPUT_DIR/bioproject_fetch_${TIMESTAMP}.tsv"
    FETCH_JSON="$OUTPUT_DIR/bioproject_fetch_${TIMESTAMP}.json"
    
    printf "${CYAN}[Fetcher]${NC} Executing cascade workflow (BioProject → SRA → PubMed)...\n\n"
    
    python3 boolean_fetcher_integrated.py "$QUERY" \
        --max "$MAX_RESULTS" \
        --output-tsv "$FETCH_TSV" \
        --output-json "$FETCH_JSON"
    
    printf "\n${GREEN}✓ Fetch completed${NC}\n"
    printf "  📄 TSV:  %s\n" "$FETCH_TSV"
    printf "  📄 JSON: %s\n" "$FETCH_JSON"
    
    # If Pro mode, run Data Analyzer
    if [ "$MODE" = "pro" ]; then
        printf "\n${CYAN}[Data Analyzer]${NC} Running AI analysis with Ollama...\n\n"
        
        cd "$ROOT_DIR/Data_Analyzer"
        ANALYSIS_TSV="$OUTPUT_DIR/classified_papers_${TIMESTAMP}.tsv"
        
        python3 paper_analyzer.py "$FETCH_JSON" "$ANALYSIS_TSV"
        
        printf "\n${GREEN}✓ Analysis completed${NC}\n"
        printf "  📊 Classified TSV: %s\n" "$ANALYSIS_TSV"
    fi
fi

# ================================================================
# STEP 7: DATA VISUALIZATION (Interactive HTML)
# ================================================================

printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
printf "${GREEN}[7/7] Generating interactive visualization...${NC}\n"
printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"

cd "$ROOT_DIR"
python3 Data_visualization/paper_visualizer.py "$FETCH_JSON"

VISUALIZATION_HTML="$OUTPUT_DIR/paper_network_${TIMESTAMP}.html"
if [ ! -f "$VISUALIZATION_HTML" ]; then
    # Fallback if timestamp doesn't match perfectly (unlikely but safe)
    VISUALIZATION_HTML=$(ls -t "$OUTPUT_DIR"/paper_network_*.html 2>/dev/null | head -n 1 || echo "")
fi

# ================================================================
# FINAL SUMMARY
# ================================================================

printf "\n${CYAN}╔═══════════════════════════════════════════════════════════════╗${NC}\n"
printf "${CYAN}║                                                               ║${NC}\n"
printf "${CYAN}║                    ${GREEN}✓ MINING COMPLETE${CYAN}                       ║${NC}\n"
printf "${CYAN}║                                                               ║${NC}\n"
printf "${CYAN}╚═══════════════════════════════════════════════════════════════╝${NC}\n\n"

printf "${BLUE}Query:${NC} %s\n" "$USER_INPUT"
printf "${BLUE}Mode:${NC}  %s\n" "$MODE"
printf "${BLUE}Database:${NC} %s\n\n" "$DATABASE"

printf "${YELLOW}Output files in: %s${NC}\n\n" "$OUTPUT_DIR"

if [ "$MODE" = "lite" ]; then
    printf "  📄 Fetch TSV:  $(basename "$FETCH_TSV")\n"
    printf "  📄 Fetch JSON: $(basename "$FETCH_JSON")\n"
else
    printf "  📄 Fetch TSV:       $(basename "$FETCH_TSV")\n"
    printf "  📄 Fetch JSON:      $(basename "$FETCH_JSON")\n"
    printf "  📊 Classified TSV:  $(basename "$ANALYSIS_TSV")\n"
fi

if [ -f "$VISUALIZATION_HTML" ]; then
    printf "  ✨ Interactive HTML: $(basename "$VISUALIZATION_HTML")\n"
fi

printf "\n${CYAN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n\n"
