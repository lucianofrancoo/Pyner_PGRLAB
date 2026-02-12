"""
Test: Enrich BioProjects with SRA experiments and BioSamples info
"""
import sys
sys.path.insert(0, 'Fetcher_NCBI')
from Bio import Entrez
from config import NCBI_EMAIL, NCBI_API_KEY, RATE_LIMIT
import time

Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

def get_related_experiments(bioproject_acc: str, max_results: int = 100) -> dict:
    """Get SRA experiments, BioSamples, and other data for a BioProject."""
    
    info = {
        "bioproject": bioproject_acc,
        "sra_experiments": [],
        "biosamples": [],
        "sra_count": 0,
        "biosample_count": 0
    }
    
    try:
        # Get SRA experiments
        print(f"\n📊 Fetching SRA experiments for {bioproject_acc}...")
        time.sleep(RATE_LIMIT)
        
        sra_search = Entrez.esearch(
            db="sra",
            term=f"{bioproject_acc}[BioProject]",
            retmax=max_results,
            rettype="json"
        )
        sra_result = Entrez.read(sra_search)
        sra_ids = sra_result.get("esearchresult", {}).get("idlist", [])
        info["sra_count"] = int(sra_result.get("esearchresult", {}).get("count", 0))
        
        if sra_ids:
            print(f"  ✅ Found {len(sra_ids[:10])} SRA experiments (total: {info['sra_count']})")
            info["sra_experiments"] = sra_ids[:10]  # Store first 10 IDs
        
        # Get BioSamples
        print(f"📊 Fetching BioSamples for {bioproject_acc}...")
        time.sleep(RATE_LIMIT)
        
        biosample_search = Entrez.esearch(
            db="biosample",
            term=f"{bioproject_acc}[BioProject]",
            retmax=max_results,
            rettype="json"
        )
        biosample_result = Entrez.read(biosample_search)
        biosample_ids = biosample_result.get("esearchresult", {}).get("idlist", [])
        info["biosample_count"] = int(biosample_result.get("esearchresult", {}).get("count", 0))
        
        if biosample_ids:
            print(f"  ✅ Found {len(biosample_ids[:10])} BioSamples (total: {info['biosample_count']})")
            info["biosamples"] = biosample_ids[:10]  # Store first 10 IDs
        
        return info
        
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return info

# Test with sample BioProjects
test_bioprojects = [
    "PRJNA1418118",  # gasa3 mutants
    "PRJNA1357100",  # ALKBH10B
    "PRJNA1302545",  # 14-3-3 proteins
]

print("=" * 70)
print("Testing BioProject enrichment with SRA/BioSample data")
print("=" * 70)

for bp in test_bioprojects:
    info = get_related_experiments(bp, max_results=50)
    print(f"\n📦 {bp}:")
    print(f"   SRA Experiments: {info['sra_count']} total")
    print(f"   BioSamples: {info['biosample_count']} total")
    if info["sra_experiments"]:
        print(f"   First SRA IDs: {info['sra_experiments'][:3]}")
    if info["biosamples"]:
        print(f"   First BioSample IDs: {info['biosamples'][:3]}")
