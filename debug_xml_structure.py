#!/usr/bin/env python
"""
Debug script to inspect the XML structure of SRA runs
"""
import sys
sys.path.insert(0, 'Fetcher_NCBI')
from ncbi_fetcher_sra_fixed import fetch_run_details
import xml.etree.ElementTree as ET

# Test with a specific run
run_accession = "SRR36541090"
print(f"Fetching XML for run: {run_accession}")
print("=" * 80)

# Get the XML
xml_data = fetch_run_details(run_accession)

if xml_data:
    print("\n✓ XML received")
    print("\nParsing XML structure...")
    
    try:
        root = ET.fromstring(xml_data)
        print(f"\nRoot tag: {root.tag}")
        
        # Show all direct children
        print("\n📋 Direct children of root:")
        for child in root:
            print(f"  - {child.tag}")
            # Show attributes
            if child.attrib:
                for key, val in child.attrib.items():
                    print(f"    @{key} = {val}")
        
        # Look for Statistics section
        print("\n🔍 Looking for Statistics...")
        for stats in root.findall('.//Statistics'):
            print(f"\nFound Statistics element:")
            for child in stats:
                print(f"  - {child.tag}: {child.attrib}")
                for subchild in child:
                    print(f"    - {subchild.tag}: {subchild.attrib}")
        
        # Look for any 'size' or 'bases' attribute
        print("\n🔍 Looking for size/bases/bases_count attributes...")
        for elem in root.iter():
            if elem.attrib:
                for key, val in elem.attrib.items():
                    if 'base' in key.lower() or 'size' in key.lower():
                        print(f"  {elem.tag} @{key} = {val}")
        
        # Print full XML for inspection
        print("\n" + "=" * 80)
        print("📄 Full XML structure (first 2000 chars):")
        print("=" * 80)
        print(xml_data[:2000])
        
    except Exception as e:
        print(f"❌ Error parsing XML: {e}")
else:
    print("❌ No XML data received")
