#!/usr/bin/env python3
"""
Fetcher_NCBI - Command Line Interface
======================================
Fetch biological sequencing data from NCBI SRA using boolean queries.

Usage:
    # Basic search
    python main.py --query "arabidopsis[Organism] AND drought"
    
    # Use query from file
    python main.py --query-file query.txt
    
    # Specify max results and output
    python main.py --query "arabidopsis AND stress" --max-results 500 --output results.json
    
    # Integration with Query Generator
    python main.py --from-query-generator "../Query_generator output"
"""

import argparse
import sys
import json
from pathlib import Path

from config import validate_config, MAX_RESULTS, DEFAULT_OUTPUT, get_credentials
from ncbi_fetcher import NCBIFetcher


def load_query_from_file(file_path: Path) -> str:
    """
    Load query from text file.
    Supports both plain text and JSON format.
    """
    with open(file_path, 'r') as f:
        content = f.read().strip()
    
    # Try to parse as JSON (from Query Generator output)
    try:
        data = json.loads(content)
        # Try common field names
        for field in ['query', 'boolean_query', 'optimized_query', 'ncbi_query']:
            if field in data:
                return data[field]
        # If JSON but no recognized field, use entire content as string
        return content
    except json.JSONDecodeError:
        # Plain text query
        return content


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Fetch sequencing data from NCBI SRA",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Direct query
  python main.py -q "arabidopsis[Organism] AND drought"
  
  # From file
  python main.py --query-file my_query.txt
  
  # Specify output and max results
  python main.py -q "stress response" -m 500 -o results.json
  
  # Integration with Query Generator (use JSON output)
  python main.py --query-file ../Query_generator/output/query.json
        """
    )
    
    # Query input options (mutually exclusive)
    query_group = parser.add_mutually_exclusive_group(required=True)
    query_group.add_argument(
        '-q', '--query',
        type=str,
        help='Boolean search query for NCBI SRA'
    )
    query_group.add_argument(
        '--query-file',
        type=Path,
        help='File containing query (plain text or JSON)'
    )
    
    # Output options
    parser.add_argument(
        '-o', '--output',
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f'Output file path (default: {DEFAULT_OUTPUT})'
    )
    
    # Fetch options
    parser.add_argument(
        '-m', '--max-results',
        type=int,
        default=MAX_RESULTS,
        help=f'Maximum results to fetch (default: {MAX_RESULTS})'
    )
    
    parser.add_argument(
        '--no-deduplicate',
        action='store_true',
        help='Disable BioProject deduplication (process all results)'
    )
    
    # Cache management
    parser.add_argument(
        '--clear-cache',
        action='store_true',
        help='Clear BioProject deduplication cache before fetching'
    )
    
    # Credentials
    parser.add_argument(
        '--email',
        type=str,
        help='Override NCBI email from config'
    )
    
    parser.add_argument(
        '--api-key',
        type=str,
        help='Override NCBI API key from config'
    )
    
    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Display header
    print("\n" + "=" * 70)
    print("FETCHER_NCBI - NCBI SRA Data Fetcher")
    print("=" * 70)
    
    # Validate configuration
    validate_config()
    
    # Get credentials
    creds = get_credentials()
    email = args.email if args.email else creds['email']
    api_key = args.api_key if args.api_key else creds['api_key']
    
    # Determine query
    if args.query:
        query = args.query
        print(f"Query: {query}")
    else:
        print(f"Loading query from: {args.query_file}")
        try:
            query = load_query_from_file(args.query_file)
            print(f"Query: {query}")
        except Exception as e:
            print(f"✗ Failed to load query from file: {e}")
            sys.exit(1)
    
    # Initialize fetcher
    print(f"\nInitializing fetcher...")
    print(f"  Email: {email}")
    print(f"  API Key: {'✓' if api_key else '✗'}")
    print(f"  Max Results: {args.max_results}")
    print(f"  Deduplication: {'✓' if not args.no_deduplicate else '✗'}")
    
    fetcher = NCBIFetcher(email=email, api_key=api_key)
    
    # Clear cache if requested
    if args.clear_cache:
        print(f"\nClearing deduplication cache...")
        fetcher.cache.clear()
        print("✓ Cache cleared")
    
    # Fetch data
    print(f"\nStarting fetch operation...\n")
    
    try:
        results = fetcher.fetch_all(
            query=query,
            max_results=args.max_results,
            deduplicate=not args.no_deduplicate
        )
        
        if results:
            # Save results
            fetcher.save_results(args.output)
            
            print(f"\n" + "=" * 70)
            print("SUCCESS")
            print("=" * 70)
            print(f"✓ Fetched {len(results)} unique results")
            print(f"✓ Saved to: {args.output}")
            print(f"✓ Log available in: Fetcher_NCBI/logs/")
            print("=" * 70 + "\n")
            
            return 0
        else:
            print(f"\n" + "=" * 70)
            print("NO RESULTS")
            print("=" * 70)
            print("No results found or all were duplicates")
            print("Try adjusting your query or using --clear-cache")
            print("=" * 70 + "\n")
            return 1
            
    except KeyboardInterrupt:
        print("\n\n✗ Interrupted by user")
        return 130
    except Exception as e:
        print(f"\n✗ Error during fetch: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
