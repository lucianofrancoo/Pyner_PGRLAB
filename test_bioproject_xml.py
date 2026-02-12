"""
Extract full XML from BioProject to find embedded links
"""
import sys
sys.path.insert(0, 'Fetcher_NCBI')
from Bio import Entrez
from config import NCBI_EMAIL, NCBI_API_KEY, RATE_LIMIT
import time
import xml.etree.ElementTree as ET

Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

# Get full XML for a BioProject
bp_id = "1418118"  # PRJNA1418118

print(f"Fetching raw XML for BioProject {bp_id}...")
time.sleep(RATE_LIMIT)

try:
    handle = Entrez.efetch(
        db="bioproject",
        id=bp_id,
        rettype="native",
        retmode="xml"
    )
    
    xml_content = handle.read().decode('utf-8')
    
    # Print interesting XML elements
    print(f"\n📄 XML Size: {len(xml_content)} bytes\n")
    
    # Parse XML
    root = ET.fromstring(xml_content)
    
    # Find all unique element tags
    print("🔍 XML Structure - All unique tags:")
    tags = set()
    for elem in root.iter():
        tags.add(elem.tag)
    
    for tag in sorted(tags):
        print(f"   - {tag}")
    
    # Look for any ID references
    print("\n🔗 Looking for data links in XML...")
    xml_str = xml_content
    
    for pattern in ["SRA", "accession", "run", "sample", "experiment", "RefSeq", "BioSample"]:
        if pattern.lower() in xml_str.lower():
            print(f"   ✅ Found '{pattern}' in XML")
            # Find context
            idx = xml_str.lower().find(pattern.lower())
            context = xml_str[max(0, idx-100):min(len(xml_str), idx+200)]
            print(f"      Context: ...{context}...")
    
    # Extract specific fields
    print("\n📋 Key information from XML:")
    
    # Pretty print first 2000 chars
    print(xml_content[:2000])
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
