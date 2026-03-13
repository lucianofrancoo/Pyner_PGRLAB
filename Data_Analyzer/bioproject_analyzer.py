#!/usr/bin/env python3
"""
BioProject Analyzer - Specialized Module for Dataset Usability
==============================================================
Transforms raw BioProject/SRA cascade search JSON into a normalized, 
experiment/run-centric TSV, evaluating data usability using Ollama.

Usage:
    python3 bioproject_analyzer.py <input_json> [output_tsv]
"""

import argparse
import csv
import json
import logging
import sys
import re
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set, Any

from da_config import OUTPUT_DIR, LOGS_DIR
from ollama_client import OllamaClient

# Setup logging
log_file_path = LOGS_DIR / f"bioproject_analyzer_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(log_file_path),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

BIOPROJECT_CSV_COLUMNS = [
    'User_Query', 'BioProject_ID', 'BioSample_ID', 'Experiment_ID', 'Run_ID',
    'Title', 'Organism', 'Project_Type', 'Library_Strategy', 'Library_Source',
    'Library_Selection', 'Layout', 'Platform', 'Instrument', 
    'Tissue', 'Condition', 'Treatment', 'Timepoint', 'Genotype',
    'Total_Samples', 'Total_Runs', 'Estimated_Volume_GB', 
    'Linked_PMIDs', 'Has_Published_Paper', 'Metadata_Completeness_Score', 
    'Relevance_Score', 'Recommended_For_Query', 'Experimental_Design_Summary', 
    'Primary_Use_Case', 'Major_Limitations', 'Confidence_Explanation', 'Flags'
]

class BioProjectNormalizer:
    """Normalizes the BioProject JSON structure into an Experiment/Run-centric list."""
    
    @staticmethod
    def _parse_size_to_gb(size_str: str) -> float:
        """Parse sizes like '1.2 GB', '500 MB', '12412412 bytes' into GB."""
        if not size_str or size_str == 'N/A':
            return 0.0
            
        s = str(size_str).upper().strip()
        try:
            if 'GB' in s or 'GIGABYTE' in s:
                return float(re.sub(r'[^0-9.]', '', s))
            elif 'MB' in s or 'MEGABYTE' in s:
                return float(re.sub(r'[^0-9.]', '', s)) / 1024.0
            elif 'KB' in s or 'KILOBYTE' in s:
                return float(re.sub(r'[^0-9.]', '', s)) / (1024.0 * 1024.0)
            elif 'B' in s or 'BYTE' in s:
                return float(re.sub(r'[^0-9.]', '', s)) / (1024.0 * 1024.0 * 1024.0)
            else:
                return float(s) / (1024.0 * 1024.0 * 1024.0)
        except:
            return 0.0

    @staticmethod
    def _extract_biological_attribute(attrs: Dict, keywords: List[str]) -> str:
        """Search BioSample attributes using keyword presence"""
        if not attrs:
            return 'not described'
        
        # Exact/Partial matches on keys
        for key, value in attrs.items():
            key_lower = key.lower()
            for kw in keywords:
                if kw in key_lower:
                    return str(value)
        return 'not described'

    def normalize(self, bp_data: Dict, user_query: str) -> List[Dict]:
        """Flatten a BioProject into a list of normalized experiment records."""
        rows = []
        bp_id = bp_data.get('bioproject', 'Unknown')
        title = bp_data.get('title', 'not described')
        organism = bp_data.get('organism', 'not described')
        project_type = bp_data.get('project_type', 'not described')
        
        pmids = bp_data.get('pmids', [])
        linked_pmids = "; ".join(pmids) if pmids else "not described"
        # Three levels of evidence as requested, here we do explicit_ncbi_link or missing
        has_paper = "explicit_ncbi_link" if pmids else "No"
        
        hierarchy = bp_data.get('sra_hierarchy', {})
        total_samples = bp_data.get('biosamples_count', 0)
        total_runs = bp_data.get('sra_runs_count', 0)
        
        # If no hierarchy, yield a Project-level row
        if not hierarchy:
            row = self._create_base_row(user_query, bp_id, title, organism, project_type, linked_pmids, has_paper)
            row['Total_Samples'] = total_samples
            row['Total_Runs'] = total_runs
            row['Flags'] = 'project_only_no_sra'
            rows.append(row)
            return rows

        # Iterate through Biosamples and Experiments
        for bs_id, bs_data in hierarchy.items():
            experiments = bs_data.get('experiments', [])
            if not experiments:
                row = self._create_base_row(user_query, bp_id, title, organism, project_type, linked_pmids, has_paper)
                row['BioSample_ID'] = bs_id
                row['Flags'] = 'sample_only_no_exp'
                rows.append(row)
                continue
                
            for exp in experiments:
                exp_id = exp.get('experiment_id', 'Unknown')
                metadata = exp.get('metadata', {})
                attrs = exp.get('sample_attributes', {})
                runs = exp.get('runs', [])
                
                # Biologic metadata from BioSample
                tissue = self._extract_biological_attribute(attrs, ['tissue', 'cell_type', 'cell type', 'organ', 'source'])
                condition = self._extract_biological_attribute(attrs, ['condition', 'stress', 'disease', 'phenotype', 'state'])
                treatment = self._extract_biological_attribute(attrs, ['treatment', 'drug', 'agent', 'exposure', 'stimulus'])
                timepoint = self._extract_biological_attribute(attrs, ['time', 'age', 'stage', 'phase', 'day', 'hour'])
                genotype = self._extract_biological_attribute(attrs, ['genotype', 'strain', 'cultivar', 'breed', 'mutant'])
                
                # Technical metadata
                strategy = metadata.get('library_strategy', 'not described')
                source = metadata.get('library_source', 'not described')
                selection = metadata.get('library_selection', 'not described')
                layout = metadata.get('library_layout', 'not described')
                instrument = metadata.get('instrument', 'not described')
                platform = "not described" # Sometimes it is part of instrument
                if instrument != 'not described':
                    if 'Illumina' in instrument: platform = 'Illumina'
                    elif 'PacBio' in instrument: platform = 'PacBio'
                    elif 'Nanopore' in instrument or 'MinION' in instrument: platform = 'Oxford Nanopore'
                    elif 'Ion Torrent' in instrument: platform = 'Ion Torrent'
                    elif 'BGI' in instrument: platform = 'BGISEQ'
                
                # Completeness heuristics (just a simple 0-10 approximation)
                fields = [tissue, condition, treatment, genotype, strategy, instrument]
                completeness = sum(1 for f in fields if f != 'not described')
                score_compl = int((completeness / len(fields)) * 10)
                
                if not runs:
                    row = self._create_base_row(user_query, bp_id, title, organism, project_type, linked_pmids, has_paper)
                    row.update({
                        'BioSample_ID': bs_id, 'Experiment_ID': exp_id,
                        'Library_Strategy': strategy, 'Library_Source': source,
                        'Library_Selection': selection, 'Layout': layout,
                        'Platform': platform, 'Instrument': instrument,
                        'Tissue': tissue, 'Condition': condition,
                        'Treatment': treatment, 'Timepoint': timepoint, 'Genotype': genotype,
                        'Metadata_Completeness_Score': score_compl,
                        'Flags': 'experiment_only_no_runs'
                    })
                    rows.append(row)
                    continue

                for r in runs:
                    run_id = r.get('accession', 'Unknown')
                    size_gb = self._parse_size_to_gb(r.get('size', '0'))
                    
                    row = self._create_base_row(user_query, bp_id, title, organism, project_type, linked_pmids, has_paper)
                    row.update({
                        'BioSample_ID': bs_id, 'Experiment_ID': exp_id, 'Run_ID': run_id,
                        'Library_Strategy': strategy, 'Library_Source': source,
                        'Library_Selection': selection, 'Layout': layout,
                        'Platform': platform, 'Instrument': instrument,
                        'Tissue': tissue, 'Condition': condition,
                        'Treatment': treatment, 'Timepoint': timepoint, 'Genotype': genotype,
                        'Estimated_Volume_GB': round(size_gb, 2),
                        'Metadata_Completeness_Score': score_compl,
                        'Flags': 'run_level_data'
                    })
                    rows.append(row)
        
        return rows

    def _create_base_row(self, user_query, bp_id, title, organism, project_type, linked_pmids, has_paper) -> Dict:
        return {
            'User_Query': user_query,
            'BioProject_ID': bp_id,
            'BioSample_ID': 'not described',
            'Experiment_ID': 'not described',
            'Run_ID': 'not described',
            'Title': title,
            'Organism': organism,
            'Project_Type': project_type,
            'Library_Strategy': 'not described',
            'Library_Source': 'not described',
            'Library_Selection': 'not described',
            'Layout': 'not described',
            'Platform': 'not described',
            'Instrument': 'not described',
            'Tissue': 'not described',
            'Condition': 'not described',
            'Treatment': 'not described',
            'Timepoint': 'not described',
            'Genotype': 'not described',
            'Total_Samples': 0,
            'Total_Runs': 0,
            'Estimated_Volume_GB': 0.0,
            'Linked_PMIDs': linked_pmids,
            'Has_Published_Paper': has_paper,
            'Metadata_Completeness_Score': 0,
            'Relevance_Score': 0,
            'Recommended_For_Query': 'not described',
            'Experimental_Design_Summary': 'not described',
            'Primary_Use_Case': 'not described',
            'Major_Limitations': 'not described',
            'Confidence_Explanation': 'not described',
            'Flags': ''
        }

class BioProjectAnalyzer:
    """Main analyzer for BioProject Data"""
    
    def __init__(self, ollama_client: OllamaClient):
        self.ollama = ollama_client
        self.normalizer = BioProjectNormalizer()
        
    def analyze_dataset(self, fetcher_data: Dict) -> List[Dict]:
        """Main pipeline to analyze the whole dataset"""
        results_raw = fetcher_data.get('results', [])
        user_query = fetcher_data.get('metadata', {}).get('query', 'Unknown query')
        
        logger.info(f"User query: {user_query}")
        logger.info("=" * 80)
        logger.info("Starting BioProject Analysis...")
        
        final_tsv_rows = []
        
        for idx, bp in enumerate(results_raw, 1):
            bp_id = bp.get('bioproject', 'Unknown')
            logger.info(f"\n[{idx}/{len(results_raw)}] Processing Project: {bp_id}")
            
            # 1. Deterministic heuristic check
            prefilter_score = self._heuristic_prefilter(bp, user_query)
            
            # 2. Invoke LLM if reasonably good (or even if bad, to get a conclusion. Here we ask LLM always for now, but pass prefilter_score as a flag maybe)
            # Actually, to save time we could skip LLM if prefilter_score is 0, but since this is phase 3, we analyze.
            logger.info(f"   → Analyzing usability with Qwen2.5 (Prefilter Baseline: {prefilter_score}/10)")
            llm_result = self.ollama.analyze_bioproject(bp, user_query)
            
            # Combine prefilter insight
            if prefilter_score == 0 and llm_result['relevance_score'] > 5:
                # LLM might be hallucinating relevance, flag it
                llm_result['confidence_explanation'] += f" [Warning: LLM score contradicts low deterministic match]"
            
            # 3. Normalize to tabular form
            normalized_rows = self.normalizer.normalize(bp, user_query)
            
            # 4. Inject LLM evaluation into every row of this project
            for row in normalized_rows:
                # Fill project-level totals if not already filled
                if row['Total_Samples'] == 0:
                    row['Total_Samples'] = bp.get('biosamples_count', 0)
                if row['Total_Runs'] == 0:
                    row['Total_Runs'] = bp.get('sra_runs_count', 0)
                
                row['Relevance_Score'] = llm_result.get('relevance_score', 0)
                row['Recommended_For_Query'] = llm_result.get('recommended_for_query', 'not described')
                row['Experimental_Design_Summary'] = llm_result.get('experimental_design_summary', 'not described')
                row['Primary_Use_Case'] = llm_result.get('primary_use_case', 'not described')
                row['Major_Limitations'] = llm_result.get('major_limitations', 'not described')
                row['Confidence_Explanation'] = llm_result.get('confidence_explanation', 'not described')
                
                # Append to TSV
                final_tsv_rows.append(row)
                
            logger.info(f"   ✓ Generated {len(normalized_rows)} normalized records.")
            logger.info(f"   ✓ Usability Score: {llm_result.get('relevance_score')}/10 ({llm_result.get('recommended_for_query')})")

        return final_tsv_rows

    def _heuristic_prefilter(self, bp_data: Dict, query: str) -> int:
        """Calculate a basic rule-based score. 0 to 10."""
        score = 0
        q = query.lower()
        
        # Checking Organism
        org = bp_data.get('organism', '').lower()
        if org and org in q: score += 3
        
        # Checking data type / strategy broadly
        desc = bp_data.get('description', '').lower()
        title = bp_data.get('title', '').lower()
        
        # Example: RNA-Seq/Transcriptome bonus
        if 'rna' in q or 'transcripto' in q:
            if 'rna-seq' in desc or 'transcripto' in desc or 'rna-seq' in title or 'transcripto' in title:
                score += 3
                
        # If it has SRA data, boost score, otherwise empty project
        if bp_data.get('sra_experiments_count', 0) > 0:
            score += 2
            
        # If it has published paper, it's more usable
        if bp_data.get('pmids', []):
            score += 2
            
        return min(score, 10)

def save_results(results: List[Dict], output_path: Path):
    """Save classified results to TSV"""
    logger.info(f"\nSaving BioProject TSV results to: {output_path}")
    
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=BIOPROJECT_CSV_COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    
    logger.info(f"✓ Saved {len(results)} normalized experimental records")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze BioProject datasets from Fetcher output using Ollama LLM"
    )
    
    parser.add_argument('input_json', nargs='?', help='Input JSON from BioProject Fetcher')
    parser.add_argument('output_tsv', nargs='?', help='Output TSV file')
    
    args = parser.parse_args()
    
    if not args.input_json:
        parser.print_help()
        sys.exit(1)
        
    input_path = Path(args.input_json)
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}")
        sys.exit(1)
        
    if args.output_tsv:
        output_path = Path(args.output_tsv)
    else:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_path = OUTPUT_DIR / f"bioproject_analysis_{timestamp}.tsv"
        
    logger.info("=" * 80)
    logger.info("BIOPROJECT DATASET ANALYZER - Powered by Ollama")
    logger.info("=" * 80)
    
    ollama = OllamaClient()
    if not ollama.is_available():
        logger.error("ERROR: Ollama is not available. Ensure ollama serve is running.")
        sys.exit(1)
        
    try:
        with open(input_path, 'r', encoding='utf-8') as f:
            fetcher_data = json.load(f)
            
        analyzer = BioProjectAnalyzer(ollama)
        tsv_results = analyzer.analyze_dataset(fetcher_data)
        save_results(tsv_results, output_path)
        
        logger.info("\n✓ BIOPROJECT ANALYSIS COMPLETE")
        logger.info(f"  Input:  {input_path}")
        logger.info(f"  Output: {output_path}")
    except Exception as e:
        logger.exception(f"Fatal error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
