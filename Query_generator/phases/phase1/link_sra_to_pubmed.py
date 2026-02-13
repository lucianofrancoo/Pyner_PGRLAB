#!/usr/bin/env python3
"""
Link SRA Experiments to PubMed Publications
============================================

Enriches SRA data with PubMed links via BioSample identifiers.

Usage:
    python link_sra_to_pubmed.py sra_PRJNA1179470_*.json
    python link_sra_to_pubmed.py --biosamples SAMN44494209 SAMN44494208
    python link_sra_to_pubmed.py --bioproject PRJNA1179470
"""

import sys
import json
import argparse
import os
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict

# Add paths to sys - support both direct execution and module import
__file_dir__ = os.path.dirname(os.path.abspath(__file__))
__project_root__ = os.path.dirname(os.path.dirname(os.path.dirname(__file_dir__)))
__fetcher_path__ = os.path.join(__project_root__, 'Fetcher_NCBI')

sys.path.insert(0, __project_root__)
sys.path.insert(0, __fetcher_path__)

from ncbi_linkout import LinkoutFetcher


def link_sra_json_file(json_file: Path, output_file: Optional[Path] = None) -> Dict:
    """
    Load SRA JSON and link to PubMed.
    
    Args:
        json_file: Path to SRA JSON file
        output_file: Optional output file path
    
    Returns:
        Enriched data with PubMed links
    """
    print(f"\n📂 Loading SRA data: {json_file}")
    
    with open(json_file) as f:
        sra_data = json.load(f)
    
    print(f"   Experiments: {len(sra_data.get('experiments', []))}")
    
    # Extract unique BioSamples
    biosamples = set()
    accessions = set()
    
    for exp in sra_data.get('experiments', []):
        biosample = exp.get('biosample', '')
        if biosample:
            biosamples.add(biosample)
        
        exp_acc = exp.get('exp_accession', '')
        if exp_acc:
            accessions.add(exp_acc)
    
    biosamples = list(biosamples)
    accessions = list(accessions)
    
    print(f"   Unique BioSamples: {len(biosamples)}")
    print(f"   Unique Accessions: {len(accessions)}")
    
    # Link to PubMed
    linkout = LinkoutFetcher()
    
    print(f"\n🔍 Searching PubMed for BioSamples...")
    pubmed_biosample_results = linkout.search_publications_for_biosamples(biosamples[:15])
    
    print(f"\n🔍 Searching PubMed for SRA Accessions...")
    pubmed_accession_results = linkout.search_publications_for_sra_accessions(accessions[:10])
    
    # Combine results
    enriched_data = sra_data.copy()
    enriched_data['pubmed_links'] = {
        'via_biosamples': pubmed_biosample_results,
        'via_accessions': pubmed_accession_results,
        'unique_pmids': list(linkout.stats["unique_pmids"]),
        'summary': {
            'biosamples_queried': len(biosamples),
            'accessions_queried': len(accessions),
            'pubmed_hits': linkout.stats['pubmed_hits'],
            'unique_publications': len(linkout.stats["unique_pmids"])
        }
    }
    
    linkout.print_summary()
    
    # Save results
    if output_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        bioproject = sra_data['metadata'].get('bioproject', 'unknown')
        output_file = Path(f"sra_pubmed_{bioproject}_{timestamp}.json")
    
    print(f"\n💾 Saving enriched data to: {output_file}")
    
    with open(output_file, 'w') as f:
        json.dump(enriched_data, f, indent=2)
    
    return enriched_data


def main():
    """Main entry point."""
    
    parser = argparse.ArgumentParser(
        description='Link SRA experiments to PubMed publications',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s sra_PRJNA1179470_*.json
  %(prog)s --biosamples SAMN44494209 SAMN44494208
  %(prog)s --bioproject PRJNA1179470
        """
    )
    
    parser.add_argument(
        'sra_file',
        nargs='?',
        help='SRA JSON file to enrich'
    )
    
    parser.add_argument(
        '--biosamples',
        nargs='+',
        help='BioSample accessions to search'
    )
    
    parser.add_argument(
        '--bioproject',
        help='BioProject accession to search'
    )
    
    parser.add_argument(
        '--output',
        help='Output file path'
    )
    
    args = parser.parse_args()
    
    linkout = LinkoutFetcher()
    
    # Handle different input types
    if args.sra_file:
        # Input: SRA JSON file
        sra_file = Path(args.sra_file)
        if not sra_file.exists():
            print(f"❌ File not found: {sra_file}")
            return 1
        
        enriched = link_sra_json_file(sra_file, Path(args.output) if args.output else None)
        
    elif args.biosamples:
        # Input: List of BioSamples
        print(f"\n🔍 Searching PubMed for BioSamples...")
        results = linkout.search_publications_for_biosamples(args.biosamples)
        
        output = {
            "metadata": {
                "search_type": "biosamples",
                "biosamples": args.biosamples,
                "fetched_at": datetime.now().isoformat(),
                "unique_pmids": len(linkout.stats["unique_pmids"])
            },
            "results": results
        }
        
        # Save
        if args.output:
            output_file = Path(args.output)
        else:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = Path(f"pubmed_biosamples_{timestamp}.json")
        
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=2)
        
        print(f"\n✓ Results saved to: {output_file}")
        
    elif args.bioproject:
        # Input: BioProject
        print(f"\n🔍 Searching PubMed for BioProject: {args.bioproject}")
        results = linkout.search_publications_for_bioproject(args.bioproject)
        
        output = {
            "metadata": {
                "search_type": "bioproject",
                "bioproject": args.bioproject,
                "fetched_at": datetime.now().isoformat(),
                "publications_found": len(results.get('publications', []))
            },
            "results": results
        }
        
        # Save
        if args.output:
            output_file = Path(args.output)
        else:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_file = Path(f"pubmed_{args.bioproject}_{timestamp}.json")
        
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=2)
        
        print(f"\n✓ Results saved to: {output_file}")
        
        # Display results
        if results.get('publications'):
            print(f"\n✅ Found {len(results['publications'])} publications!")
            for i, pub in enumerate(results['publications'][:3], 1):
                print(f"\n{i}. {pub.get('title', '')[:70]}...")
                print(f"   Authors: {', '.join(pub.get('authors', []))}")
                print(f"   Year: {pub.get('year', '')}")
                print(f"   PMID: {pub.get('pmid', '')}")
        
    else:
        parser.print_help()
        return 1
    
    linkout.print_summary()
    return 0


if __name__ == "__main__":
    sys.exit(main())
