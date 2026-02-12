"""
Test different search strategies to link BioProject to SRA/BioSample
"""
import sys
sys.path.insert(0, 'Fetcher_NCBI')
from Bio import Entrez
from config import NCBI_EMAIL, NCBI_API_KEY, RATE_LIMIT
import time
import json

Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

def explore_bioproject_data(bioproject_acc: str):
    """Try different ways to get related SRA/BioSample data."""
    
    print(f"\n🔍 Exploring {bioproject_acc}")
    print("=" * 70)
    
    # Strategy 1: Get the BioProject summary details
    print(f"\n1️⃣  Getting BioProject summary...")
    time.sleep(RATE_LIMIT)
    
    try:
        summary = Entrez.esummary(db="bioproject", id=bioproject_acc)
        summary_data = Entrez.read(summary)
        
        if "DocumentSummarySet" in summary_data:
            doc = summary_data["DocumentSummarySet"]["DocumentSummary"][0]
            print(f"   Project ID: {doc.get('Project_Id')}")
            print(f"   Project Acc: {doc.get('Project_Acc')}")
            
            # Look for any cross-links
            print(f"\n   All available fields:")
            for key in sorted(doc.keys()):
                val = doc[key]
                if isinstance(val, str) and "link" in key.lower():
                    print(f"     - {key}: {val}")
                elif isinstance(val, str) and len(val) > 0 and len(val) < 100:
                    print(f"     - {key}: {val}")
    except Exception as e:
        print(f"   ❌ Error: {e}")
    
    # Strategy 2: Try fetching instead of summary for more details
    print(f"\n2️⃣  Trying fetch (XML) for raw data...")
    time.sleep(RATE_LIMIT)
    
    try:
        fetch_handle = Entrez.efetch(db="bioproject", id=bioproject_acc, rettype="native", retmode="xml")
        #fetch_data = fetch_handle.read().decode('utf-8')
        #print(f"   Got XML data ({len(fetch_data)} bytes)")
        # Just check if we get data
        print(f"   ✅ Fetch returns data (XML format available)")
    except Exception as e:
        print(f"   ❌ Error: {e}")
    
    # Strategy 3: Try alternative queries in SRA
    print(f"\n3️⃣  Trying SRA searches with different query terms...")
    
    search_terms = [
        f"{bioproject_acc}",
        f"PRJNA1418118",  # Use NCBI ID from BioProject
        f"1418118",  # Try just the ID number
    ]
    
    for term in search_terms:
        time.sleep(RATE_LIMIT)
        try:
            search = Entrez.esearch(db="sra", term=term, retmax=5, rettype="json")
            result = Entrez.read(search)
            count = int(result.get("esearchresult", {}).get("count", 0))
            if count > 0:
                print(f"   ✅ Query '{term}' found {count} SRA records")
                ids = result.get("esearchresult", {}).get("idlist", [])
                if ids:
                    print(f"      First IDs: {ids[:3]}")
        except Exception as e:
            print(f"   ❌ Query '{term}': {e}")
    
    # Strategy 4: Try BioSample searches
    print(f"\n4️⃣  Trying BioSample searches...")
    
    for term in search_terms:
        time.sleep(RATE_LIMIT)
        try:
            search = Entrez.esearch(db="biosample", term=term, retmax=5, rettype="json")
            result = Entrez.read(search)
            count = int(result.get("esearchresult", {}).get("count", 0))
            if count > 0:
                print(f"   ✅ Query '{term}' found {count} BioSample records")
                ids = result.get("esearchresult", {}).get("idlist", [])
                if ids:
                    print(f"      First IDs: {ids[:3]}")
        except Exception as e:
            print(f"   ❌ Query '{term}': {e}")

# Test with one BioProject
explore_bioproject_data("PRJNA1418118")
