"""
Data Analyzer - Configuration
==============================
Configuration for paper analysis with Ollama LLM
"""

from pathlib import Path

# ============================================
# PATHS
# ============================================
DATA_ANALYZER_DIR = Path(__file__).parent
ROOT_DIR = DATA_ANALYZER_DIR.parent
OUTPUT_DIR = DATA_ANALYZER_DIR / "output"
LOGS_DIR = DATA_ANALYZER_DIR / "logs"

# Create directories
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

# ============================================
# OLLAMA CONFIGURATION
# ============================================
OLLAMA_BASE_URL = "http://localhost:11434"
OLLAMA_MODEL = "qwen2.5:14b"
OLLAMA_TIMEOUT = 120  # seconds

# ============================================
# ANALYSIS CONFIGURATION
# ============================================
# Relevance threshold (0-10 scale)
RELEVANCE_THRESHOLD = 5

# Maximum abstract length to send to LLM (characters)
MAX_ABSTRACT_LENGTH = 3000

# Use full text from PMC when available (PMCID)
USE_PMC_FULL_TEXT = True

# Maximum full text length to send to LLM (0 = no limit)
# Set to 0 to use complete text from PMC
MAX_FULL_TEXT_LENGTH = 0  # No limit - use complete text

# Batch size for processing papers
BATCH_SIZE = 10

# ============================================
# OUTPUT CONFIGURATION
# ============================================
CSV_DELIMITER = ","
MULTIVALUE_SEPARATOR = " ; "  # Separator for multiple organisms/tissues/conditions

# CSV columns order - EXPANDED with detailed experimental metadata
CSV_COLUMNS = [
    # Identification
    "PMID",
    "PMCID",
    "Title",
    "Year",
    "Journal",
    "DOI",
    
    # Relevance
    "Relevance_Score",
    "Is_Relevant",
    
    # Organism Details
    "Organisms",
    "Species",
    "Strain_Variety",
    "Genotype",
    "Tissues_Organs",
    "Source_Tissue_Origin",
    "Cell_Type",
    
    # Developmental / Temporal Biology
    "Developmental_Stage",
    "Organism_Age",
    "Growth_Phase",
    
    # Sample Collection Conditions
    "Conditions",
    "Environmental_Stress",
    "Temperature_Range",
    "Light_Conditions",
    "Growth_Medium",
    "Sample_Collection_Conditions",
    
    # Molecular Analysis
    "Molecules_Extracted",
    "RNA_Type",
    "DNA_Type",
    "Protein_Type",
    "Other_Molecules",
    
    # Experimental Design
    "Strategies",
    "Measurement_Tools",
    "Detection_Method",
    
    # Time Course
    "Time_Course_Design",
    "Time_Points",
    "Time_Intervals",
    "Time_Duration",
    
    # Replication
    "Sample_Size",
    "Biological_Replicates",
    "Technical_Replicates",
    "Replication_Design",
    
    # Treatment / Comparisons
    "Treatment_Groups",
    "Control_Type",
    "Dose_Range",
    
    # Quality & Contamination
    "Quality_Metrics",
    "Contamination_Check",
    
    # Biological Context
    "Pathway_Focus",
    "Biomarkers_Measured",
    "Disease_Model",
    "Differential_Expression_Threshold",
    
    # Data Availability
    "Normalization_Method",
    "Statistical_Method",
    "Raw_Data_Available",
    
    # Preview
    "Abstract_Preview"
]

# ============================================
# LOGGING
# ============================================
from datetime import datetime
LOG_FILE = LOGS_DIR / f"analyzer_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

print(f"✅ Data Analyzer Config loaded")
print(f"   Ollama: {OLLAMA_BASE_URL}/{OLLAMA_MODEL}")
print(f"   Output: {OUTPUT_DIR}")
