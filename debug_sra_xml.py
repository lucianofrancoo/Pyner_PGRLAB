#!/usr/bin/env python
"""
Debug script to see what's in the SRA XML for run statistics
"""
import sys
sys.path.insert(0, 'Fetcher_NCBI')

from Bio import Entrez
from config import NCBI_EMAIL, NCBI_API_KEY
import xml.etree.ElementTree as ET

Entrez.email = NCBI_EMAIL
Entrez.api_key = NCBI_API_KEY

# Get one of the experiments
experiment_id = 42185164  # First experiment

print(f"Fetching SRA XML for experiment ID: {experiment_id}...")
print("=" * 70)

try:
    # Fetch using Entrez
    handle = Entrez.efetch(db='sra', id=experiment_id, rettype='xml', retmode='xml')
    xml_data = handle.read()
    handle.close()
    
    root = ET.fromstring(xml_data)
    
    # Find all runs
    runs = root.findall('.//RUN')
    print(f"Found {len(runs)} runs in XML\n")
    
    for i, run in enumerate(runs[:2]):  # Check first 2 runs
        run_acc = run.get('accession', 'N/A')
        print(f"\nRun #{i+1}: {run_acc}")
        print("-" * 40)
        
        # Check Statistics
        stats = run.find('.//Statistics')
        if stats is not None:
            print("✓ Statistics element found")
            # Look for reads
            reads = stats.findall('.//Read')
            print(f"  - Found {len(reads)} Read elements")
            for j, read in enumerate(reads[:3]):
                count = read.get('count', 'N/A')
                bases = read.get('bases', 'N/A')
                print(f"    Read {j}: count={count}, bases={bases}")
        else:
            print("✗ No Statistics element found")
        
        # Check what elements are in the run
        children = [child.tag for child in run]
        print(f"  - Run contains: {children}")
            
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
