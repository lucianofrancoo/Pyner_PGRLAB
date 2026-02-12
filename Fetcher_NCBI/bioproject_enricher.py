"""
Function to enrich BioProjects with SRA experiment counts and BioSample info
by searching in those databases
"""

import time
from typing import Dict, Optional, List
from Bio import Entrez
from config import NCBI_EMAIL, NCBI_API_KEY, RATE_LIMIT

Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

def enrich_bioproject_with_sra(bioproject_data: Dict) -> Dict:
    """
    Enrich BioProject data with SRA experiment count and BioSample count.
    
    Searches SRA and BioSample databases using the BioProject accession.
    
    Args:
        bioproject_data: Dictionary with at least 'bioproject' key
        
    Returns:
        Same dict with added fields:
        - sra_experiments_count: Number of SRA experiments
        - biosamples_count: Number of BioSamples
        - sra_experiments: List of first 10 SRA IDs (optional, for details)
    """
    
    bp_acc = bioproject_data.get("bioproject", "")
    if not bp_acc:
        return bioproject_data
    
    # Add new fields
    bioproject_data["sra_experiments_count"] = 0
    bioproject_data["biosamples_count"] = 0
    
    try:
        # Try to find SRA experiments
        # The connection is often through project title/accession appearing in run metadata
        time.sleep(RATE_LIMIT)
        
        sra_search = Entrez.esearch(
            db="sra",
            term=f"{bp_acc}",  # Search for BioProject accession in SRA
            rettype="json",
            retmax=1
        )
        sra_result = Entrez.read(sra_search)
        sra_count = int(sra_result.get("esearchresult", {}).get("count", 0))
        bioproject_data["sra_experiments_count"] = sra_count
        
    except Exception as e:
        pass  # Keep default 0
    
    try:
        # Try to find BioSamples
        time.sleep(RATE_LIMIT)
        
        bs_search = Entrez.esearch(
            db="biosample",
            term=f"{bp_acc}",  # Search for BioProject accession in BioSample
            rettype="json",
            retmax=1
        )
        bs_result = Entrez.read(bs_search)
        bs_count = int(bs_result.get("esearchresult", {}).get("count", 0))
        bioproject_data["biosamples_count"] = bs_count
        
    except Exception as e:
        pass  # Keep default 0
    
    return bioproject_data


# Test
if __name__ == "__main__":
    print("Testing SRA/BioSample enrichment...\n")
    
    # Sample BioProject data
    test_data = [
        {"bioproject": "PRJNA1418118", "title": "RNA-Seq analyses of gasa3 mutants"},
        {"bioproject": "PRJNA1357100", "title": "ALKBH10B confers drought tolerance"},
        {"bioproject": "PRJNA1302545", "title": "14-3-3 proteins GRF6 and GRF8"},
    ]
    
    for data in test_data:
        print(f"🔍 {data['bioproject']}")
        enriched = enrich_bioproject_with_sra(data)
        print(f"   SRA experiments: {enriched['sra_experiments_count']}")
        print(f"   BioSamples: {enriched['biosamples_count']}")
        print()
