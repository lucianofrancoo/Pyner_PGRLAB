#!/usr/bin/env python3
"""
Paper Analyzer - Main Script
=============================
Analyzes papers from Fetcher output using Ollama LLM

Reads JSON output from Fetcher_NCBI (PubMed results) and:
1. Evaluates relevance of each paper to user query
2. Extracts structured information (organisms, tissues, conditions, strategies)
3. Generates classified table in CSV format

Usage:
    python3 paper_analyzer.py <input_json> [output_csv]
    
Example:
    python3 paper_analyzer.py ../pubmed_results_20260213_113559.json
"""

import argparse
import csv
import json
import logging
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List

from config import (
    OUTPUT_DIR, LOGS_DIR, LOG_FILE, 
    CSV_COLUMNS, CSV_DELIMITER, MULTIVALUE_SEPARATOR,
    RELEVANCE_THRESHOLD, MAX_ABSTRACT_LENGTH, USE_PMC_FULL_TEXT
)
from ollama_client import OllamaClient
from pmc_fetcher import PMCFullTextFetcher

# ============================================
# LOGGING
# ============================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)


class PaperAnalyzer:
    """Main analyzer class"""
    
    def __init__(self, ollama_client: OllamaClient):
        self.ollama = ollama_client
        self.pmc_fetcher = PMCFullTextFetcher() if USE_PMC_FULL_TEXT else None
        self.stats = {
            'total': 0,
            'analyzed': 0,
            'relevant': 0,
            'errors': 0,
            'pmc_full_text': 0,
            'abstract_only': 0
        }
    
    def load_fetcher_output(self, json_path: Path) -> Dict:
        """Load JSON output from Fetcher_NCBI"""
        logger.info(f"Loading fetcher output: {json_path}")
        
        with open(json_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        self.stats['total'] = len(data.get('publications', []))
        logger.info(f"Loaded {self.stats['total']} papers")
        
        return data
    
    def analyze_papers(self, fetcher_data: Dict) -> List[Dict]:
        """Analyze all papers and return classified results"""
        publications = fetcher_data.get('publications', [])
        user_query = fetcher_data.get('metadata', {}).get('query', 'Unknown query')
        
        logger.info(f"User query: {user_query[:200]}...")
        logger.info("=" * 80)
        logger.info("Starting paper analysis with Ollama...")
        logger.info("=" * 80)
        
        results = []
        
        for idx, paper in enumerate(publications, 1):
            logger.info(f"\n[{idx}/{self.stats['total']}] Analyzing: {paper.get('pmid', 'N/A')}")
            logger.info(f"   Title: {paper.get('title', 'No title')[:80]}...")
            
            try:
                classified = self._analyze_single_paper(paper, user_query)
                results.append(classified)
                self.stats['analyzed'] += 1
                
                if classified['Is_Relevant'] == 'Yes':
                    self.stats['relevant'] += 1
                
                logger.info(f"   ✓ Relevance: {classified['Relevance_Score']}/10")
                logger.info(f"   ✓ Organisms: {classified['Organisms'][:60]}")
                logger.info(f"   ✓ Tissues: {classified['Tissues'][:60]}")
                logger.info(f"   ✓ Conditions: {classified['Conditions'][:60]}")
                
            except Exception as e:
                logger.error(f"   ✗ Error analyzing paper: {e}")
                self.stats['errors'] += 1
                # Add empty result
                results.append(self._create_error_result(paper))
        
        return results
    
    def _analyze_single_paper(self, paper: Dict, user_query: str) -> Dict:
        """Analyze single paper with Ollama"""
        title = paper.get('title', '')
        abstract = paper.get('abstract', '')[:MAX_ABSTRACT_LENGTH]
        pmcid = paper.get('pmcid', None)
        
        # Try to fetch full text from PMC if PMCID available
        full_text = None
        if self.pmc_fetcher and pmcid:
            logger.info(f"   → Fetching full text from PMC: {pmcid}")
            sections = self.pmc_fetcher.fetch_full_text(pmcid)
            if sections and 'full_text_preview' in sections:
                full_text = sections['full_text_preview']
                self.stats['pmc_full_text'] += 1
                logger.info(f"   ✓ Using PMC full text ({len(full_text)} chars)")
            else:
                logger.warning(f"   ⚠ PMC fetch failed, using abstract only")
                self.stats['abstract_only'] += 1
        else:
            self.stats['abstract_only'] += 1
        
        # Call Ollama (with full text if available)
        analysis = self.ollama.analyze_paper(title, abstract, user_query, full_text=full_text)
        
        # Build result row
        result = {
            'PMID': paper.get('pmid', 'N/A'),
            'Title': title,
            'Relevance_Score': analysis['relevance_score'],
            'Is_Relevant': 'Yes' if analysis['relevance_score'] >= RELEVANCE_THRESHOLD else 'No',
            'Organisms': MULTIVALUE_SEPARATOR.join(analysis['organisms']) if analysis['organisms'] else 'N/A',
            'Tissues': MULTIVALUE_SEPARATOR.join(analysis['tissues']) if analysis['tissues'] else 'N/A',
            'Conditions': MULTIVALUE_SEPARATOR.join(analysis['conditions']) if analysis['conditions'] else 'N/A',
            'Strategies': MULTIVALUE_SEPARATOR.join(analysis['strategies']) if analysis['strategies'] else 'N/A',
            'Year': paper.get('year', 'N/A'),
            'Journal': paper.get('journal', 'N/A'),
            'DOI': paper.get('doi', 'N/A'),
            'Abstract_Preview': abstract[:200] + '...' if len(abstract) > 200 else abstract
        }
        
        return result
    
    def _create_error_result(self, paper: Dict) -> Dict:
        """Create result for papers that failed analysis"""
        return {
            'PMID': paper.get('pmid', 'N/A'),
            'Title': paper.get('title', 'N/A'),
            'Relevance_Score': 0,
            'Is_Relevant': 'Error',
            'Organisms': 'ERROR',
            'Tissues': 'ERROR',
            'Conditions': 'ERROR',
            'Strategies': 'ERROR',
            'Year': paper.get('year', 'N/A'),
            'Journal': paper.get('journal', 'N/A'),
            'DOI': paper.get('doi', 'N/A'),
            'Abstract_Preview': 'Analysis failed'
        }
    
    def save_results(self, results: List[Dict], output_path: Path):
        """Save classified results to CSV"""
        logger.info(f"\nSaving results to: {output_path}")
        
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, delimiter=CSV_DELIMITER)
            writer.writeheader()
            writer.writerows(results)
        
        logger.info(f"✓ Saved {len(results)} classified papers")
    
    def print_statistics(self):
        """Print analysis statistics"""
        print("\n" + "=" * 80)
        print("ANALYSIS STATISTICS")
        print("=" * 80)
        print(f"Total papers:        {self.stats['total']}")
        print(f"Successfully analyzed: {self.stats['analyzed']}")
        print(f"Relevant papers:     {self.stats['relevant']} ({self.stats['relevant']/max(1, self.stats['total'])*100:.1f}%)")
        print(f"Analysis errors:     {self.stats['errors']}")
        if self.pmc_fetcher:
            print(f"PMC full text used: {self.stats['pmc_full_text']}")
            print(f"Abstract only:       {self.stats['abstract_only']}")
        print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze papers from Fetcher output using Ollama LLM",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 paper_analyzer.py ../pubmed_results_20260213_113559.json
  python3 paper_analyzer.py ../pubmed_results_20260213_113559.json my_analysis.csv
  python3 paper_analyzer.py --test-connection
        """
    )
    
    parser.add_argument('input_json', nargs='?', help='Input JSON from Fetcher_NCBI')
    parser.add_argument('output_csv', nargs='?', help='Output CSV file (optional, auto-generated if not provided)')
    parser.add_argument('--test-connection', action='store_true', help='Test Ollama connection and exit')
    
    args = parser.parse_args()
    
    # Test connection mode
    if args.test_connection:
        from ollama_client import test_ollama_connection
        sys.exit(0 if test_ollama_connection() else 1)
    
    # Validate input
    if not args.input_json:
        parser.print_help()
        sys.exit(1)
    
    input_path = Path(args.input_json)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}")
        sys.exit(1)
    
    # Generate output path if not provided
    if args.output_csv:
        output_path = Path(args.output_csv)
    else:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_path = OUTPUT_DIR / f"classified_papers_{timestamp}.csv"
    
    # Initialize Ollama client
    print("\n" + "=" * 80)
    print("PAPER ANALYZER - Powered by Ollama")
    print("=" * 80)
    
    ollama = OllamaClient()
    if not ollama.is_available():
        print("ERROR: Ollama is not available")
        print("Make sure Ollama is running: ollama serve")
        print("And that the model is pulled: ollama pull qwen2.5:14b")
        sys.exit(1)
    
    print(f"✓ Ollama connected")
    print(f"✓ Model: {ollama.model}")
    
    # Run analysis
    analyzer = PaperAnalyzer(ollama)
    
    try:
        # Load data
        fetcher_data = analyzer.load_fetcher_output(input_path)
        
        # Analyze papers
        results = analyzer.analyze_papers(fetcher_data)
        
        # Save results
        analyzer.save_results(results, output_path)
        
        # Print statistics
        analyzer.print_statistics()
        
        print(f"\n✓ ANALYSIS COMPLETE")
        print(f"  Input:  {input_path}")
        print(f"  Output: {output_path}")
        print(f"  Log:    {LOG_FILE}")
        print("=" * 80 + "\n")
        
    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.exception("Fatal error during analysis")
        print(f"\nERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
