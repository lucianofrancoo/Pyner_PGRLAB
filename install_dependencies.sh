#!/bin/bash
#
# PYNER - Quick Installation Script
# ==================================
# Installs minimal dependencies required to run Pyner workflows
#

set -e

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "\n${BLUE}================================================================${NC}"
echo -e "${BLUE}         PYNER - INSTALLATION SCRIPT${NC}"
echo -e "${BLUE}================================================================${NC}\n"

# Check Python 3
echo -e "${YELLOW}[1/3] Checking Python 3...${NC}"
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version)
    echo -e "${GREEN}✓ Found: $PYTHON_VERSION${NC}"
else
    echo -e "${RED}✗ Python 3 not found${NC}"
    echo -e "  Install with: ${YELLOW}sudo apt install python3 python3-pip${NC}"
    exit 1
fi

# Check pip
echo -e "\n${YELLOW}[2/3] Checking pip...${NC}"
if command -v pip3 &> /dev/null || command -v pip &> /dev/null; then
    echo -e "${GREEN}✓ pip is available${NC}"
    PIP_CMD=$(command -v pip3 || command -v pip)
else
    echo -e "${RED}✗ pip not found${NC}"
    echo -e "  Install with: ${YELLOW}sudo apt install python3-pip${NC}"
    exit 1
fi

# Install BioPython
echo -e "\n${YELLOW}[3/3] Installing BioPython...${NC}"
if python3 -c "import Bio" &> /dev/null; then
    echo -e "${GREEN}✓ BioPython already installed${NC}"
else
    echo -e "  Installing biopython..."
    $PIP_CMD install biopython --quiet
    
    if python3 -c "import Bio" &> /dev/null; then
        echo -e "${GREEN}✓ BioPython installed successfully${NC}"
    else
        echo -e "${RED}✗ BioPython installation failed${NC}"
        echo -e "  Try manually: ${YELLOW}pip install biopython${NC}"
        exit 1
    fi
fi

# Success
echo -e "\n${BLUE}================================================================${NC}"
echo -e "${GREEN}✓ Installation Complete!${NC}\n"
echo -e "You can now run the workflow:"
echo -e "  ${YELLOW}bash test_fetcher_integrator.sh${NC}\n"
echo -e "Optional: Configure NCBI API key for faster fetching"
echo -e "  Edit: ${YELLOW}Fetcher_NCBI/config.py${NC}"
echo -e "  Or set: ${YELLOW}export NCBI_API_KEY=your_key${NC}\n"
echo -e "${BLUE}================================================================${NC}\n"
