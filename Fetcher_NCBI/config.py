"""
Fetcher_NCBI Configuration
==========================
Manages NCBI Entrez API credentials and fetcher settings.
"""

import os
from pathlib import Path

# ============================================
# NCBI ENTREZ API CONFIGURATION
# ============================================

# NCBI requires email for all Entrez requests
# Override with NCBI_EMAIL environment variable
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "lucianofranco.a@gmail.com")

# API Key increases rate limits from 3/sec to 10/sec
# Get your key at: https://www.ncbi.nlm.nih.gov/account/settings/
# Override with NCBI_API_KEY environment variable
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "8f86973f360600c58b8bb20daf1728a4be08")

# ============================================
# FETCHER SETTINGS
# ============================================

# Maximum results to fetch per query
# None means fetch all results (may take a long time)
MAX_RESULTS = None

# Stop after collecting this many unique BioProjects
MIN_UNIQUE_BIOPROJECTS = 30

# Database to search (bioproject or sra)
DATABASE = "bioproject"

# Rate limiting (seconds between requests)
# With API key: 0.1s (10 req/sec)
# Without API key: 0.34s (3 req/sec)
RATE_LIMIT = 0.1 if NCBI_API_KEY else 0.34

# ============================================
# FILE PATHS
# ============================================

# Base directory for Fetcher_NCBI
BASE_DIR = Path(__file__).parent

# Data storage directory
DATA_DIR = BASE_DIR / "data"
DATA_DIR.mkdir(exist_ok=True)

# Logs directory
LOGS_DIR = BASE_DIR / "logs"
LOGS_DIR.mkdir(exist_ok=True)

# Default output file for search results
DEFAULT_OUTPUT = DATA_DIR / "bioproject_results.json"

# Deduplication cache (tracks processed BioProjects)
DEDUP_CACHE = DATA_DIR / "bioproject_cache.json"

# ============================================
# VALIDATION
# ============================================

def validate_config():
    """
    Validate configuration and warn about missing credentials.
    Returns True if configuration is valid, False otherwise.
    """
    issues = []
    
    if NCBI_EMAIL == "your.email@example.com":
        issues.append("⚠️  NCBI_EMAIL not configured. Set via environment variable or edit config.py")
    
    if not NCBI_API_KEY:
        issues.append("⚠️  NCBI_API_KEY not set. Rate limited to 3 req/sec (vs 10 with key)")
    
    if issues:
        print("\n" + "=" * 70)
        print("CONFIGURATION WARNINGS")
        print("=" * 70)
        for issue in issues:
            print(issue)
        print("=" * 70 + "\n")
        return False
    
    return True


def get_credentials():
    """
    Get NCBI credentials as a dictionary.
    """
    return {
        "email": NCBI_EMAIL,
        "api_key": NCBI_API_KEY
    }


if __name__ == "__main__":
    # Test configuration
    print("Fetcher_NCBI Configuration")
    print("=" * 40)
    print(f"NCBI Email: {NCBI_EMAIL}")
    print(f"API Key: {'✓ Set' if NCBI_API_KEY else '✗ Not Set'}")
    print(f"Max Results: {MAX_RESULTS}")
    print(f"Database: {DATABASE}")
    print(f"Rate Limit: {RATE_LIMIT}s")
    print(f"Data Directory: {DATA_DIR}")
    print(f"Logs Directory: {LOGS_DIR}")
    print("\n")
    validate_config()
