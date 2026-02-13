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
import re
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
    
    # Common experimental techniques to look for in post-processing
    TECHNIQUE_KEYWORDS = {
        'qRT-PCR': [r'\bq(?:rt|uantitative\s+real.?time\s+pcr)\b', r'\bqrt[- ]?pcr\b', r'\bqpcr\b'],
        'RT-PCR': [r'\brt[- ]?pcr\b', r'\breverse\s+transcription\s+pcr\b'],
        'RNA-Seq': [r'\brna[- ]?seq(?:uencing)?\b', r'\brnaseq\b'],
        'Microarray': [r'\bmicroarray\b', r'\bgene\s+chip\b'],
        'Western blot': [r'\bwestern\s+blot(?:ting)?\b'],
        'Northern blot': [r'\bnorthern\s+blot(?:ting)?\b'],
        'PCR': [r'\bpcr\b'],
        'ChIP-Seq': [r'\bchip[- ]?seq\b', r'\bchip\s+sequencing\b'],
        'Sequencing': [r'\bsequencing\b', r'\bnext[- ]?gen(?:eration)?\s+sequencing\b', r'\bngs\b'],
        'Microscopy': [r'\bmicroscop(?:y|ic)\b', r'\bconfocal\b', r'\bfluorescence\b'],
        'Histological staining': [r'\bhistolog(?:y|ical)\b', r'\bistochemical?\s+staining\b'],
        'ELISA': [r'\belisa\b', r'\bimmuno[- ]?assay\b'],
        'Flow cytometry': [r'\bflow\s+cytometr(?:y|ic)\b'],
        'Bioinformatics': [r'\bbioinformatics\b', r'\bcomputational\s+analysis\b'],
        'Mass spectrometry': [r'\bmass\s+spectrom(?:etry|etric)\b'],
        'Proteomics': [r'\bproteomics?\b'],
        'Metabolomics': [r'\bmetabolomics?\b'],
        'Phenotyping': [r'\bphenotyping\b', r'\bmorphological\s+analysis\b'],
        'Gas exchange': [r'\bgas\s+exchange\b'],
        'Enzyme assay': [r'\benzyme\s+assay\b'],
        'Antioxidant assay': [r'\bantioxidant\s+assay\b'],
    }
    
    def __init__(self, ollama_client: OllamaClient):
        self.ollama = ollama_client
        self.pmc_fetcher = PMCFullTextFetcher() if USE_PMC_FULL_TEXT else None
        self.stats = {
            'total': 0,
            'analyzed': 0,
            'relevant': 0,
            'errors': 0,
            'pmc_full_text': 0,
            'abstract_only': 0,
            'techniques_enhanced': 0  # Papers where we found missing techniques
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
                    
            except Exception as e:
                logger.error(f"   ✗ Error analyzing paper: {e}")
                self.stats['errors'] += 1
                # Create and append error result when analysis fails
                results.append(self._create_error_result(paper))
                continue  # Skip logging below for error rows
            
            # Log successful results (outside try/except to avoid duplicates)
            try:
                logger.info(f"   ✓ Relevance: {classified['Relevance_Score']}/10")
                logger.info(f"   ✓ Organisms: {classified['Organisms'][:60]}")
                logger.info(f"   ✓ Tissues_Organs: {classified['Tissues_Organs'][:60]}")
                logger.info(f"   ✓ Conditions: {classified['Conditions'][:60]}")
            except Exception as log_err:
                logger.warning(f"   ⚠ Could not log all fields: {log_err}")
        
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
        
        # Post-process: search for additional techniques in text
        # Even if LLM found some, we want to make sure we don't miss any
        text_to_search = full_text if full_text else abstract
        found_techniques_regex = self._extract_techniques_regex(text_to_search)
        
        if found_techniques_regex:
            # Merge LLM results with regex results
            llm_strategies = set(analysis['strategies']) if analysis['strategies'] else set()
            # Remove empty/placeholder values
            llm_strategies = {s for s in llm_strategies if s and s not in ['N/A', 'n/a', '']}
            
            all_strategies = llm_strategies.union(set(found_techniques_regex))
            
            if len(all_strategies) > len(llm_strategies):
                # Found new techniques via regex
                logger.info(f"   ℹ Enhanced techniques: {found_techniques_regex}")
                analysis['strategies'] = sorted(list(all_strategies))
                self.stats['techniques_enhanced'] += 1
        
        # Helper function to format array fields
        def format_array(items):
            """Convert list to semicolon-separated string"""
            if isinstance(items, list):
                items = [str(i).strip() for i in items if i and str(i).strip() not in ['N/A', 'n/a', '']]
                return MULTIVALUE_SEPARATOR.join(items) if items else 'not described'
            return str(items) if items else 'not described'
        
        # Helper function to format scalar fields
        def format_scalar(value):
            """Convert scalar to string, or 'not described' if empty"""
            if not value or value == 'N/A' or value == 'n/a':
                return 'not described'
            return str(value).strip()
        
        # Build result row with EXPANDED METADATA
        result = {
            # Identification
            'PMID': paper.get('pmid', 'N/A'),
            'PMCID': paper.get('pmcid', 'N/A'),
            'Title': title,
            'Year': paper.get('year', 'N/A'),
            'Journal': paper.get('journal', 'N/A'),
            'DOI': paper.get('doi', 'N/A'),
            
            # Relevance
            'Relevance_Score': analysis['relevance_score'],
            'Is_Relevant': 'Yes' if analysis['relevance_score'] >= RELEVANCE_THRESHOLD else 'No',
            
            # Organism Details (EXPANDED)
            'Organisms': format_array(analysis.get('organisms', [])),
            'Species': format_scalar(analysis.get('species', '')),
            'Strain_Variety': format_scalar(analysis.get('strain_variety', '')),
            'Genotype': format_scalar(analysis.get('genotype', '')),
            'Tissues_Organs': format_array(analysis.get('tissues_organs', [])),
            'Source_Tissue_Origin': format_scalar(analysis.get('source_tissue_origin', '')),
            'Cell_Type': format_scalar(analysis.get('cell_type', '')),
            
            # Developmental Biology (EXPANDED)
            'Developmental_Stage': format_scalar(analysis.get('developmental_stage', '')),
            'Organism_Age': format_scalar(analysis.get('organism_age', '')),
            'Growth_Phase': format_scalar(analysis.get('growth_phase', '')),
            
            # Sample Collection (EXPANDED)
            'Conditions': format_array(analysis.get('conditions', [])),
            'Environmental_Stress': format_scalar(analysis.get('environmental_stress', '')),
            'Temperature_Range': format_scalar(analysis.get('temperature_range', '')),
            'Light_Conditions': format_scalar(analysis.get('light_conditions', '')),
            'Growth_Medium': format_scalar(analysis.get('growth_medium', '')),
            'Sample_Collection_Conditions': format_scalar(analysis.get('sample_collection_conditions', '')),
            
            # Molecular Analysis (EXPANDED)
            'Molecules_Extracted': format_array(analysis.get('molecules_extracted', [])),
            'RNA_Type': format_scalar(analysis.get('rna_type', '')),
            'DNA_Type': format_scalar(analysis.get('dna_type', '')),
            'Protein_Type': format_scalar(analysis.get('protein_type', '')),
            'Other_Molecules': format_scalar(analysis.get('other_molecules', '')),
            
            # Experimental Techniques (EXPANDED)
            'Strategies': format_array(analysis.get('strategies', [])),
            'Measurement_Tools': format_scalar(analysis.get('measurement_tools', '')),
            'Detection_Method': format_scalar(analysis.get('detection_method', '')),
            
            # Time Course (EXPANDED)
            'Time_Course_Design': format_scalar(analysis.get('time_course_design', '')),
            'Time_Points': format_scalar(analysis.get('time_points', '')),
            'Time_Intervals': format_scalar(analysis.get('time_intervals', '')),
            'Time_Duration': format_scalar(analysis.get('time_duration', '')),
            
            # Replication (EXPANDED)
            'Sample_Size': format_scalar(analysis.get('sample_size', '')),
            'Biological_Replicates': format_scalar(analysis.get('biological_replicates', '')),
            'Technical_Replicates': format_scalar(analysis.get('technical_replicates', '')),
            'Replication_Design': format_scalar(analysis.get('replication_design', '')),
            
            # Treatment / Comparisons (EXPANDED)
            'Treatment_Groups': format_scalar(analysis.get('treatment_groups', '')),
            'Control_Type': format_scalar(analysis.get('control_type', '')),
            'Dose_Range': format_scalar(analysis.get('dose_range', '')),
            
            # Quality & Contamination (EXPANDED)
            'Quality_Metrics': format_scalar(analysis.get('quality_metrics', '')),
            'Contamination_Check': format_scalar(analysis.get('contamination_check', '')),
            
            # Biological Context (EXPANDED)
            'Pathway_Focus': format_scalar(analysis.get('pathway_focus', '')),
            'Biomarkers_Measured': format_scalar(analysis.get('biomarkers_measured', '')),
            'Disease_Model': format_scalar(analysis.get('disease_model', '')),
            
            # Data Analysis (EXPANDED)
            'Normalization_Method': format_scalar(analysis.get('normalization_method', '')),
            'Statistical_Method': format_scalar(analysis.get('statistical_method', '')),
            'Differential_Expression_Threshold': format_scalar(analysis.get('differential_expression_threshold', '')),
            
            # Data Availability (EXPANDED)
            'Raw_Data_Available': format_scalar(analysis.get('raw_data_available', '')),
            
            # Preview
            'Abstract_Preview': abstract[:200] + '...' if len(abstract) > 200 else abstract
        }
        
        return result
    
    def _create_error_result(self, paper: Dict) -> Dict:
        """Create result for papers that failed analysis"""
        # Create error row with all columns set to error/default values
        error_dict = {
            'PMID': paper.get('pmid', 'N/A'),
            'PMCID': paper.get('pmcid', 'N/A'),
            'Title': paper.get('title', 'N/A'),
            'Year': paper.get('year', 'N/A'),
            'Journal': paper.get('journal', 'N/A'),
            'DOI': paper.get('doi', 'N/A'),
            'Relevance_Score': 0,
            'Is_Relevant': 'Error',
            'Organisms': 'ERROR',
            'Species': 'ERROR',
            'Strain_Variety': 'ERROR',
            'Genotype': 'ERROR',
            'Tissues_Organs': 'ERROR',
            'Source_Tissue_Origin': 'ERROR',
            'Cell_Type': 'ERROR',
            'Developmental_Stage': 'ERROR',
            'Organism_Age': 'ERROR',
            'Growth_Phase': 'ERROR',
            'Conditions': 'ERROR',
            'Environmental_Stress': 'ERROR',
            'Temperature_Range': 'ERROR',
            'Light_Conditions': 'ERROR',
            'Growth_Medium': 'ERROR',
            'Sample_Collection_Conditions': 'ERROR',
            'Molecules_Extracted': 'ERROR',
            'RNA_Type': 'ERROR',
            'DNA_Type': 'ERROR',
            'Protein_Type': 'ERROR',
            'Other_Molecules': 'ERROR',
            'Strategies': 'ERROR',
            'Measurement_Tools': 'ERROR',
            'Detection_Method': 'ERROR',
            'Time_Course_Design': 'ERROR',
            'Time_Points': 'ERROR',
            'Time_Intervals': 'ERROR',
            'Time_Duration': 'ERROR',
            'Sample_Size': 'ERROR',
            'Biological_Replicates': 'ERROR',
            'Technical_Replicates': 'ERROR',
            'Replication_Design': 'ERROR',
            'Treatment_Groups': 'ERROR',
            'Control_Type': 'ERROR',
            'Dose_Range': 'ERROR',
            'Quality_Metrics': 'ERROR',
            'Contamination_Check': 'ERROR',
            'Pathway_Focus': 'ERROR',
            'Biomarkers_Measured': 'ERROR',
            'Disease_Model': 'ERROR',
            'Normalization_Method': 'ERROR',
            'Statistical_Method': 'ERROR',
            'Differential_Expression_Threshold': 'ERROR',
            'Raw_Data_Available': 'ERROR',
            'Abstract_Preview': 'Analysis failed'
        }
        return error_dict
    
    def _extract_techniques_regex(self, text: str) -> List[str]:
        """
        Extract experimental techniques from text using regex patterns
        Used as fallback when LLM doesn't find techniques
        """
        found_techniques = set()
        text_lower = text.lower()
        
        try:
            for technique, patterns in self.TECHNIQUE_KEYWORDS.items():
                for pattern in patterns:
                    try:
                        if re.search(pattern, text_lower, re.IGNORECASE):
                            found_techniques.add(technique)
                            break  # Found this technique, move to next
                    except re.error as e:
                        logger.warning(f"Regex error in pattern '{pattern}' for '{technique}': {e}")
                        continue
        except Exception as e:
            logger.error(f"Error extracting techniques: {e}")
            return []
        
        # Sort for consistency
        return sorted(list(found_techniques)) if found_techniques else []
    
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
        print(f"Techniques enhanced: {self.stats['techniques_enhanced']} (found via text search)")
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
