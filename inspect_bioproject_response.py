"""
Inspect raw NCBI BioProject response to see what metadata is actually available.
"""
import sys
sys.path.insert(0, 'Fetcher_NCBI')
from Bio import Entrez
from config import NCBI_EMAIL, NCBI_API_KEY, RATE_LIMIT
import time

# Set email
Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

# Query for a specific BioProject
query = "arabidopsis[Organism] AND drought"
print(f"Searching BioProject with: {query}\n")
print("=" * 80)

try:
    # Search
    search_handle = Entrez.esearch(
        db="bioproject",
        term=query,
        retmax=1
    )
    search_result = Entrez.read(search_handle)
    ids = search_result["IdList"]
    
    if not ids:
        print("No results found")
        sys.exit(1)
    
    bioproject_id = ids[0]
    print(f"\nFetching summary for BioProject ID: {bioproject_id}")
    print("=" * 80)
    
    time.sleep(RATE_LIMIT)
    
    # Fetch summary - get raw XML to see structure
    summary_handle = Entrez.esummary(
        db="bioproject",
        id=bioproject_id,
        rettype="xml"
    )
    summary_result = Entrez.read(summary_handle)
    
    print("\n📋 PARSED SUMMARY STRUCTURE:")
    print("-" * 80)
    
    print(f"Top level keys: {list(summary_result.keys())}")
    
    if "DocumentSummarySet" in summary_result:
        docset = summary_result["DocumentSummarySet"]
        print(f"\nDocumentSummarySet keys: {list(docset.keys())}")
        
        if "DocumentSummary" in docset:
            docs = docset["DocumentSummary"]
            if isinstance(docs, list):
                print(f"DocumentSummary is a list with {len(docs)} items")
                doc = docs[0]
            else:
                print(f"DocumentSummary is a single item")
                doc = docs
            
            print(f"\n📍 First DocumentSummary attributes:")
            print(f"   Attributes: {doc.attributes if hasattr(doc, 'attributes') else 'N/A'}")
            
            print(f"\n📍 First DocumentSummary fields ({len(doc)} fields):")
            for i, key in enumerate(sorted(doc.keys())):
                value = doc[key]
                if isinstance(value, str):
                    val_preview = value[:80].replace('\n', ' ')
                    print(f"  {key}: {val_preview}")
                elif isinstance(value, list):
                    print(f"  {key}: [LIST with {len(value)} items]")
                    if value and isinstance(value[0], dict):
                        print(f"       → First item has keys: {list(value[0].keys())}")
                elif isinstance(value, dict):
                    print(f"  {key}: [DICT with keys: {list(value.keys())}]")
                else:
                    print(f"  {key}: {value}")

except Exception as e:
    print(f"\n❌ Error: {e}")
    import traceback
    traceback.print_exc()
