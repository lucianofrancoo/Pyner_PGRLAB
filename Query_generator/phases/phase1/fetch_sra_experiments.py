#!/usr/bin/env python3
"""
Pyner Phase 1 - SRA Data Fetcher Integration
=============================================
Integrates corrected SRA fetcher into Phase 1 workflow.

Usage:
    python fetch_sra_experiments.py PRJNA1179470
    python fetch_sra_experiments.py PRJNA1179470 --max 50 --format csv
    python fetch_sra_experiments.py PRJNA1179470 --output results/sra_data.json
"""

import sys
import argparse
import json
from pathlib import Path
from datetime import datetime

# Add Fetcher_NCBI to path (go up 3 levels to project root, then into Fetcher_NCBI)
fetcher_path = Path(__file__).parent.parent.parent.parent / 'Fetcher_NCBI'
sys.path.insert(0, str(fetcher_path))

from ncbi_fetcher_sra_fixed import SRAFetcher


def main():
    """Main entry point for SRA fetcher."""
    
    parser = argparse.ArgumentParser(
        description='Fetch SRA experiments for a BioProject from NCBI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s PRJNA1179470
  %(prog)s PRJNA1179470 --max 100 --format csv
  %(prog)s PRJNA1179470 --output /tmp/results.json
        """
    )
    
    parser.add_argument(
        'bioproject',
        help='BioProject accession (e.g., PRJNA1179470)'
    )
    
    parser.add_argument(
        '--max',
        type=int,
        default=None,
        help='Maximum experiments to fetch per BioProject (default: all)'
    )
    
    parser.add_argument(
        '--format',
        choices=['json', 'csv'],
        default='json',
        help='Output format (default: json)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Output file path (default: stdout or auto-generated)'
    )
    
    parser.add_argument(
        '--no-save',
        action='store_true',
        help='Don\'t save results to file (just display)'
    )
    
    parser.add_argument(
        '--show-all',
        action='store_true',
        help='Show all results (not just first 3)'
    )
    
    args = parser.parse_args()
    
    # Validate BioProject ID
    if not args.bioproject.startswith('PRJNA'):
        print(f"❌ Invalid BioProject ID: {args.bioproject}")
        print("   BioProject IDs start with 'PRJNA'")
        sys.exit(1)
    
    # Initialize fetcher
    print(f"\n🔍 Fetching SRA experiments for {args.bioproject}...")
    print("=" * 70)
    
    fetcher = SRAFetcher()
    
    # Fetch data
    results = fetcher.fetch_all_by_bioproject(
        args.bioproject,
        max_per_bioproject=args.max
    )
    
    if not results:
        print("❌ No results found")
        return 1
    
    # Display results
    print(f"\n✅ Found {len(results)} experiments\n")
    
    show_count = len(results) if args.show_all else min(3, len(results))
    print(f"📊 Showing {show_count} of {len(results)} experiments:\n")
    
    for i, exp in enumerate(results[:show_count], 1):
        print(f"{i}. {exp.get('exp_accession', '?')}")
        print(f"   Title: {exp.get('title', '?')[:60]}...")
        print(f"   Organism: {exp.get('organism', '?')}")
        print(f"   Strategy: {exp.get('library_strategy', '?')}")
        print(f"   Runs: {exp.get('run_count', 0)}")
        if exp.get('runs'):
            print(f"   Sample runs: {', '.join(exp['runs'][:2])}")
        print()
    
    if len(results) > show_count:
        print(f"   ... and {len(results) - show_count} more\n")
    
    # Save results
    if not args.no_save:
        if args.output:
            output_file = Path(args.output)
        else:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"sra_{args.bioproject}_{timestamp}.{args.format}"
            output_file = Path(filename)
        
        if args.format == 'csv':
            fetcher.save_results_csv(output_file)
        else:
            fetcher.save_results(output_file)
        
        print(f"✅ Results saved to: {output_file}\n")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
