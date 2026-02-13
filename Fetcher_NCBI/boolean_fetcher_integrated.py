#!/usr/bin/env python3
"""
Boolean Search → BioProject Fetcher → SRA Experiments → PubMed Linking
=====================================================================

Complete workflow:
1. Search BioProjects using boolean query
2. For each BioProject:
   - Fetch SRA experiments and BioSamples
   - Search for associated publications in PubMed
   - If not found, cascade search through BioSamples and SRA accessions
3. Export enriched results to CSV

Usage:
    python boolean_fetcher_integrated.py "Arabidopsis phosphate" --output results.csv
"""

import sys
import json
import time
import csv
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

from Bio import Entrez

from config import NCBI_EMAIL, NCBI_API_KEY, RATE_LIMIT
from ncbi_fetcher_sra_fixed import SRAFetcher, extract_sra_experiment_metadata
from ncbi_linkout import LinkoutFetcher, extract_pubmed_metadata

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

Entrez.email = NCBI_EMAIL
if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY
    logger.info("✓ Using NCBI API key")


class BooleanFetcherIntegrated:
    """
    Orchestrates complete workflow: Boolean→BioProject→SRA→PubMed
    """
    
    def __init__(self):
        self.sra_fetcher = SRAFetcher()
        self.linkout_fetcher = LinkoutFetcher()
        self.results = []
        
    def search_bioproject_boolean(self, query: str, retmax: int = 100) -> List[Dict]:
        """
        Search BioProject database using boolean query.
        
        Args:
            query: Boolean search query (e.g., "Arabidopsis AND phosphate")
            retmax: Maximum results to return
            
        Returns:
            List of BioProject records with IDs and metadata
        """
        logger.info(f"\n{'='*70}")
        logger.info(f"🔍 BOOLEAN SEARCH on BioProject")
        logger.info(f"{'='*70}")
        logger.info(f"Query: {query}")
        
        bioprojects = []
        
        try:
            time.sleep(RATE_LIMIT)
            
            # Search BioProject
            search_result = Entrez.esearch(
                db="bioproject",
                term=query,
                retmax=retmax
            )
            search_data = Entrez.read(search_result)
            bioproject_ids = search_data.get("IdList", [])
            
            logger.info(f"✓ Found {len(bioproject_ids)} BioProjects")
            
            # Fetch details for each BioProject
            if bioproject_ids:
                time.sleep(RATE_LIMIT)
                
                fetch_result = Entrez.efetch(
                    db="bioproject",
                    id=",".join(bioproject_ids[:retmax]),
                    rettype="xml"
                )
                
                from xml.etree import ElementTree as ET
                xml_text = fetch_result.read()
                root = ET.fromstring(xml_text)
                
                # Parse each DocumentSummary (not Package)
                for doc in root.findall(".//DocumentSummary"):
                    bp_data = {
                        'bioproject': '',
                        'title': '',
                        'description': '',
                        'submission_date': '',
                        'organism': '',
                        'project_type': '',
                    }
                    
                    # Extract BioProject ID
                    archive_id = doc.find(".//ProjectID/ArchiveID")
                    if archive_id is not None and 'accession' in archive_id.attrib:
                        bp_data['bioproject'] = archive_id.attrib['accession']
                    
                    # Extract Title
                    title = doc.findtext(".//Title", "")
                    if title:
                        bp_data['title'] = title
                    
                    # Extract Description
                    desc = doc.findtext(".//Description", "")
                    if desc:
                        bp_data['description'] = desc
                    
                    # Extract Organism Name
                    organism = doc.findtext(".//Name", "")
                    if organism:
                        bp_data['organism'] = organism
                    
                    # Extract Project Type
                    proj_type = doc.findtext(".//ProjectType", "")
                    if proj_type:
                        bp_data['project_type'] = proj_type
                    
                    # Extract Submission/Registration Date
                    reg_date = doc.findtext(".//RegistrationDate", "")
                    if reg_date:
                        bp_data['submission_date'] = reg_date
                    
                    if bp_data['bioproject']:
                        bioprojects.append(bp_data)
                        logger.info(f"   - {bp_data['bioproject']}: {bp_data['title'][:60]}")
            
        except Exception as e:
            logger.error(f"Error searching BioProject: {e}")
        
        return bioprojects
    
    
    def fetch_sra_for_bioproject(self, bioproject_id: str) -> Tuple[List, List]:
        """
        Fetch all SRA experiments and BioSamples for a BioProject with full metadata.
        
        Returns:
            (experiments_list, biosamples_dict, sra_runs_list)
        """
        logger.info(f"   📊 Fetching SRA data for {bioproject_id}...")
        
        experiments = []
        biosamples_dict = {}  # {biosample_id: full metadata dict}
        sra_runs = []  # Collect all SRR codes
        
        try:
            # Use SRAFetcher to get all experiments
            experiments = self.sra_fetcher.fetch_all_by_bioproject(
                bioproject_id,
                max_per_bioproject=None
            )
            
            # Extract BioSamples and SRA runs with full metadata
            for exp in experiments:
                # Collect SRA runs (SRR codes)
                if 'runs' in exp and exp['runs']:
                    sra_runs.extend(exp['runs'])
                
                # Build biosamples dict with complete metadata
                if 'biosample' in exp and exp['biosample']:
                    bs_id = exp['biosample']
                    # Store biosample with all experiment details
                    if bs_id not in biosamples_dict:
                        biosamples_dict[bs_id] = {
                            'experiments': []  # List of experiments for this sample
                        }
                    
                    # Add experiment info to this sample
                    exp_info = {
                        'exp_accession': exp.get('exp_accession', ''),
                        'title': exp.get('exp_title', ''),
                        'library_name': exp.get('library_name', ''),
                        'library_strategy': exp.get('library_strategy', ''),
                        'library_source': exp.get('library_source', ''),
                        'library_selection': exp.get('library_selection', ''),
                        'library_layout': exp.get('library_layout', ''),
                        'instrument': exp.get('instrument', ''),
                        'runs': exp.get('runs', []),
                        'runs_details': exp.get('runs_details', []),  # Run details with spots and size
                        'sample_attributes': exp.get('sample_attributes', {})
                    }
                    biosamples_dict[bs_id]['experiments'].append(exp_info)
            
            logger.info(f"      ✓ {len(experiments)} experiments, {len(biosamples_dict)} biosamples, {len(set(sra_runs))} SRA runs")
            
        except Exception as e:
            logger.error(f"      Error fetching SRA: {e}")
        
        return experiments, biosamples_dict, sra_runs
    
    
    def search_pubmed_publications(
        self,
        bioproject_id: str,
        biosamples_dict: dict,
        sra_runs: list,
        cascade: bool = True
    ) -> Tuple[List[Dict], str]:
        """
        Search for PubMed publications with cascade strategy.
        
        Returns:
            (publications_list, search_strategy_label)
        """
        logger.info(f"   🔗 Searching PubMed for {bioproject_id}...")
        
        publications = []
        search_stat = "NA"
        
        # Level 1: Direct BioProject search
        logger.info(f"      Level 1: Searching BIOPROJECT directly...")
        try:
            result = self.linkout_fetcher.search_publications_for_bioproject(bioproject_id)
            pubs = result.get('publications', [])
            if pubs:
                publications = pubs
                search_stat = "direct"
                logger.info(f"         ✓ Found {len(pubs)} publications (direct)")
                return publications, search_stat
            logger.info(f"         ✗ No results for direct bioproject search")
        except Exception as e:
            logger.info(f"         ✗ Error: {e}")
        
        # If cascade enabled, try BioSamples and SRA runs
        if cascade and not publications:
            # Level 2: BioSample search
            if biosamples_dict:
                biosamples_list = list(biosamples_dict.keys())[:10]  # Limit to 10
                logger.info(f"      Level 2: Searching {len(biosamples_list)} BIOSAMPLES...")
                try:
                    bs_results = self.linkout_fetcher.search_publications_for_biosamples(biosamples_list)
                    
                    # Flatten results from all biosamples
                    pubs = []
                    for bs_id, bs_pubs in bs_results.items():
                        pubs.extend(bs_pubs)
                    
                    if pubs:
                        publications = pubs
                        search_stat = "biosamples"
                        logger.info(f"         ✓ Found {len(pubs)} publications (biosamples)")
                        return publications, search_stat
                    logger.info(f"         ✗ No biosamples found in PubMed")
                except Exception as e:
                    logger.info(f"         ✗ Error: {e}")
            
            # Level 3: SRA runs search
            if sra_runs:
                sra_runs_list = list(set(sra_runs))[:10]  # Limit to 10, remove duplicates
                logger.info(f"      Level 3: Searching {len(sra_runs_list)} SRA RUNS...")
                try:
                    acc_results = self.linkout_fetcher.search_publications_for_sra_accessions(sra_runs_list)
                    
                    # Flatten results from all accessions
                    pubs = []
                    for acc_id, acc_pubs in acc_results.items():
                        pubs.extend(acc_pubs)
                    
                    if pubs:
                        publications = pubs
                        search_stat = "sra_runs"
                        logger.info(f"         ✓ Found {len(pubs)} publications (SRA runs)")
                        return publications, search_stat
                    logger.info(f"         ✗ No SRA runs found in PubMed")
                except Exception as e:
                    logger.info(f"         ✗ Error: {e}")
        
        if not publications:
            logger.info(f"      ⊘ No publications found")
        
        return publications, search_stat
    
    
    def process_bioproject(self, bp_data: Dict) -> Dict:
        """
        Complete workflow for a single BioProject:
        1. Fetch SRA experiments and BioSamples with full metadata
        2. Build hierarchical structure with titles
        3. Search PubMed with cascade strategy
        4. Enrich results
        
        Returns:
            Enriched BioProject data with publication info and SRA hierarchy
        """
        bioproject_id = bp_data['bioproject']
        logger.info(f"\n{'='*70}")
        logger.info(f"📋 Processing: {bioproject_id}")
        logger.info(f"{'='*70}")
        
        result = bp_data.copy()
        
        # Fetch SRA data with full metadata
        experiments, biosamples_dict, sra_runs = self.fetch_sra_for_bioproject(bioproject_id)
        
        result['sra_experiments_count'] = len(experiments)
        result['biosamples_count'] = len(biosamples_dict)
        result['sra_runs_count'] = len(set(sra_runs))  # Count unique runs
        
        # Store IDs lists
        result['sra_experiments'] = [exp.get('exp_accession', '') for exp in experiments if exp.get('exp_accession')]
        result['biosamples'] = sorted(list(biosamples_dict.keys()))
        result['sra_runs'] = sorted(list(set(sra_runs)))  # Unique SRR codes
        
        # Build hierarchical SRA structure with all metadata and titles
        result['sra_hierarchy'] = self.build_hierarchical_sra_structure(experiments, biosamples_dict)
        
        # Log hierarchy summary
        logger.info(f"   🏗️ SRA Hierarchy built:")
        for bs_id, bs_data in result['sra_hierarchy'].items():
            exp_count = len(bs_data.get('experiments', []))
            logger.info(f"      • {bs_id}: {exp_count} experiment(s)")
        
        # PubMed search disabled - only generate SRA hierarchy data
        # To enable PubMed cascade search, uncomment the code below
        
        # publications, search_method = self.search_pubmed_publications(
        #     bioproject_id,
        #     biosamples_dict,
        #     sra_runs,
        #     cascade=True
        # )
        
        # No PubMed search - set all publication fields to NA
        result['publications_found'] = 0
        result['search_method'] = "NA"
        result['dois'] = "NA"
        result['pmids'] = "NA"
        result['papers_summary'] = "NA"
        
        return result
    
    
    def run_workflow(self, boolean_query: str, max_bioproject: int = 50) -> List[Dict]:
        """
        Execute complete workflow:
        1. Boolean search
        2. For each BioProject: fetch SRA + PubMed search
        3. Return enriched results
        """
        logger.info(f"\n{'*'*70}")
        logger.info(f"🚀 BOOLEAN FETCHER INTEGRATED WORKFLOW")
        logger.info(f"{'*'*70}")
        logger.info(f"Query: {boolean_query}")
        logger.info(f"Max BioProjects: {max_bioproject}")
        
        # Step 1: Boolean search
        bioprojects = self.search_bioproject_boolean(boolean_query, retmax=max_bioproject)
        
        if not bioprojects:
            logger.warning("No BioProjects found")
            return []
        
        # Step 2: Process each BioProject
        logger.info(f"\n{'='*70}")
        logger.info(f"⚙️ PROCESSING {len(bioprojects)} BioProjects")
        logger.info(f"{'='*70}")
        
        self.results = []
        for i, bp_data in enumerate(bioprojects, 1):
            logger.info(f"\n[{i}/{len(bioprojects)}]")
            
            try:
                enriched = self.process_bioproject(bp_data)
                self.results.append(enriched)
            except Exception as e:
                logger.error(f"Error processing {bp_data['bioproject']}: {e}")
                bp_data['error'] = str(e)
                self.results.append(bp_data)
            
            time.sleep(RATE_LIMIT)
        
        return self.results
    
    
    def save_results_csv(self, output_file: Path):
        """Save results to CSV format."""
        if not self.results:
            logger.warning("No results to save")
            return
        
        logger.info(f"\n{'='*70}")
        logger.info(f"💾 SAVING RESULTS")
        logger.info(f"{'='*70}")
        
        # Define CSV columns
        fieldnames = [
            'bioproject',
            'title',
            'submission_date',
            'organism',
            'project_type',
            'description',
            'sra_experiments_count',
            'biosamples_count',
            'sra_runs_count',
            'sra_experiments',
            'biosamples',
            'sra_runs',
            'publications_found',
            'search_method',
            'papers_summary',
            'dois',
            'pmids'
        ]
        
        try:
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames, restval='NA')
                writer.writeheader()
                
                for result in self.results:
                    # Format lists as semicolon-separated strings
                    row = {}
                    for k in fieldnames:
                        value = result.get(k, 'NA')
                        if isinstance(value, list):
                            row[k] = "; ".join(value) if value else 'NA'
                        else:
                            row[k] = value if value else 'NA'
                    writer.writerow(row)
            
            logger.info(f"✓ Saved {len(self.results)} results to: {output_file}")
            logger.info(f"  - BioProjects with papers: {sum(1 for r in self.results if r.get('publications_found', 0) > 0)}")
            logger.info(f"  - BioProjects without papers: {sum(1 for r in self.results if r.get('publications_found', 0) == 0)}")
            
        except Exception as e:
            logger.error(f"Error saving CSV: {e}")
    
    
    def save_results_json(self, output_file: Path):
        """Save results to JSON format with hierarchical SRA structure."""
        if not self.results:
            logger.warning("No results to save")
            return
        
        try:
            # Transform results to include hierarchical SRA data
            enhanced_results = []
            for result in self.results:
                enhanced = result.copy()
                
                # If we have SRA hierarchy info, structure it properly
                if 'sra_hierarchy' in result:
                    enhanced['sra_hierarchy'] = result['sra_hierarchy']
                
                enhanced_results.append(enhanced)
            
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump({
                    'metadata': {
                        'total_results': len(enhanced_results),
                        'date': datetime.now().isoformat(),
                        'with_publications': sum(1 for r in enhanced_results if r.get('publications_found', 0) > 0),
                        'without_publications': sum(1 for r in enhanced_results if r.get('publications_found', 0) == 0)
                    },
                    'results': enhanced_results
                }, f, indent=2, default=str)
            
            logger.info(f"✓ Saved JSON to: {output_file}")
        except Exception as e:
            logger.error(f"Error saving JSON: {e}")
    
    
    def build_hierarchical_sra_structure(self, experiments: List, biosamples_dict: Dict) -> Dict:
        """
        Build hierarchical SRA structure showing:
        BioProject → BioSample → Experiment (with titles) → Runs (with spots & size)
        
        Returns:
            Dictionary with organized hierarchy and full metadata
        """
        hierarchy = {}
        
        for bs_id, bs_data in biosamples_dict.items():
            hierarchy[bs_id] = {
                'sample_id': bs_id,
                'experiments': []
            }
            
            # Add all experiments for this sample
            for exp_info in bs_data.get('experiments', []):
                # Build runs list with detailed metadata (spots, size)
                runs_with_details = []
                runs_details = exp_info.get('runs_details', [])
                
                if runs_details:
                    # Use detailed run information if available
                    for run_detail in runs_details:
                        runs_with_details.append({
                            'accession': run_detail.get('accession', ''),
                            'spots': run_detail.get('spots', 0),
                            'size': run_detail.get('size', '0 MB')
                        })
                else:
                    # Fallback: just list the run accessions
                    for run_acc in exp_info.get('runs', []):
                        runs_with_details.append({
                            'accession': run_acc,
                            'spots': 0,
                            'size': 'N/A'
                        })
                
                exp_entry = {
                    'experiment_id': exp_info.get('exp_accession', ''),
                    'title': exp_info.get('title', ''),
                    'metadata': {
                        'library_name': exp_info.get('library_name', ''),
                        'library_strategy': exp_info.get('library_strategy', ''),
                        'library_source': exp_info.get('library_source', ''),
                        'library_selection': exp_info.get('library_selection', ''),
                        'library_layout': exp_info.get('library_layout', ''),
                        'instrument': exp_info.get('instrument', '')
                    },
                    'sample_attributes': exp_info.get('sample_attributes', {}),
                    'runs': runs_with_details  # Now includes spots and size!
                }
                hierarchy[bs_id]['experiments'].append(exp_entry)
        
        return hierarchy


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Boolean search → BioProject → SRA → PubMed integrated fetcher"
    )
    parser.add_argument("query", help="Boolean search query (e.g., 'Arabidopsis phosphate')")
    parser.add_argument("--max", type=int, default=50, help="Max BioProjects to process (default: 50)")
    parser.add_argument("--output-csv", type=Path, help="Output CSV file")
    parser.add_argument("--output-json", type=Path, help="Output JSON file")
    
    args = parser.parse_args()
    
    # Create output files if not specified
    if not args.output_csv and not args.output_json:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        args.output_csv = Path(f"boolean_results_{timestamp}.csv")
    
    # Run workflow
    fetcher = BooleanFetcherIntegrated()
    results = fetcher.run_workflow(args.query, max_bioproject=args.max)
    
    # Save results
    if args.output_csv:
        fetcher.save_results_csv(args.output_csv)
    
    if args.output_json:
        fetcher.save_results_json(args.output_json)
    
    logger.info(f"\n{'*'*70}")
    logger.info(f"✓ WORKFLOW COMPLETE")
    logger.info(f"{'*'*70}\n")


if __name__ == "__main__":
    main()
