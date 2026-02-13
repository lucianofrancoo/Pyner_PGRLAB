"""
NCBI SRA Fetcher - FIXED VERSION
=================================
Properly fetches SRA experiment metadata using XML (not JSON).
Critical fix: Use efetch with rettype="xml" instead of esummary for SRA.
"""

import time
import json
import logging
import csv
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Optional
import xml.etree.ElementTree as ET
from Bio import Entrez
from config import (
    NCBI_EMAIL, NCBI_API_KEY, RATE_LIMIT, LOGS_DIR
)


# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(log_file: Optional[str] = None) -> logging.Logger:
    """Configure logging for the fetcher."""
    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = LOGS_DIR / f"sra_fetcher_{timestamp}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


# ============================================
# METADATA EXTRACTION FROM SRA XML
# ============================================

def extract_sra_experiment_metadata(xml_string: str, sra_id: str) -> Optional[Dict]:
    """
    Extract metadata from SRA XML for a single experiment.
    
    Args:
        xml_string: Complete SRA XML response
        sra_id: SRA runs ID
    
    Returns:
        Dictionary with extracted metadata or None if parsing fails
    """
    try:
        root = ET.fromstring(xml_string)
    except ET.ParseError as e:
        logging.error(f"XML parsing error: {e}")
        return None
    
    # Find EXPERIMENT element
    exp = root.find('.//EXPERIMENT')
    if exp is None:
        return None
    
    # Extract basic experiment info (accession is an attribute!)
    exp_accession = exp.get('accession', '').strip()
    exp_title = exp.findtext('TITLE', '').strip()
    
    # Extract run info with detailed metadata
    # Runs are peer elements of EXPERIMENT, not children
    runs = []
    runs_details = []
    
    for run in root.findall('.//RUN'):
        run_acc = run.get('accession', '')
        if run_acc:
            runs.append(run_acc)
            
            # Extract run details: spots and size
            run_info = {
                'accession': run_acc,
                'spots': 0,
                'size': '0 MB'
            }
            
            # Try to get spots and base_count from Statistics
            base_count_total = 0
            stats = run.find('.//Statistics')
            if stats is not None:
                # Look for Read elements with spot counts
                for read_elem in stats.findall('.//Read'):
                    spot_count = read_elem.get('count', '')
                    base_count = read_elem.get('bases', '')
                    
                    if spot_count:
                        try:
                            run_info['spots'] += int(spot_count)
                        except ValueError:
                            pass
                    
                    if base_count:
                        try:
                            base_count_total += int(base_count)
                        except ValueError:
                            pass
            
            # Convert base_count to human readable size (MB)
            if base_count_total > 0:
                size_mb = base_count_total / (1024 ** 2)
                run_info['size'] = f"{size_mb:.2f} MB"
            
            runs_details.append(run_info)
    
    # Extract study info
    study = root.find('.//STUDY')
    if study is not None:
        study_accession = study.get('accession', '')
        # BioProject is an external ID 
        bioproject = study.findtext('.//EXTERNAL_ID[@namespace="BioProject"]', '')
    else:
        study_accession = ''
        bioproject = ''
    
    # Extract library descriptor
    lib_desc = exp.find('.//LIBRARY_DESCRIPTOR')
    library_strategy = lib_desc.findtext('LIBRARY_STRATEGY', '').strip() if lib_desc is not None else ''
    library_source = lib_desc.findtext('LIBRARY_SOURCE', '').strip() if lib_desc is not None else ''
    library_selection = lib_desc.findtext('LIBRARY_SELECTION', '').strip() if lib_desc is not None else ''
    library_name = lib_desc.findtext('LIBRARY_NAME', '').strip() if lib_desc is not None else ''
    
    # Determine library layout
    if lib_desc is not None:
        library_layout = (
            "PAIRED" if lib_desc.find('LIBRARY_LAYOUT/PAIRED') is not None 
            else "SINGLE"
        )
    else:
        library_layout = "UNKNOWN"
    
    # Extract instruments (sequencing platform)
    instruments = []
    for instr in root.findall('.//INSTRUMENT'):
        instr_text = instr.tag if instr.tag else ''
        if instr.tag and instr.tag.upper() not in ['INSTRUMENT']:
            instruments.append(instr.tag)
    # Alternative: look for instrument in Submission block
    if not instruments:
        for submission in root.findall('.//SUBMISSION/ORGANISM'):
            instr = submission.findtext('INSTRUMENT', '')
            if instr:
                instruments.append(instr)
    # Try common locations
    instrument_name = ''
    for instr_elem in root.findall('.//*'):
        if instr_elem.tag in ['INSTRUMENT_MODEL', 'INSTRUMENT', 'MODEL']:
            if instr_elem.text:
                instrument_name = instr_elem.text.strip()
                break
    
    # Extract sample info
    sample = root.find('.//SAMPLE')
    if sample is not None:
        biosample = sample.findtext('.//EXTERNAL_ID[@namespace="BioSample"]', '')
        organism = sample.findtext('.//SCIENTIFIC_NAME', '').strip()
    else:
        biosample = ''
        organism = ''
    
    # Extract sample attributes (metadata) if available
    sample_attributes = {}
    if sample is not None:
        for attr in sample.findall('.//SAMPLE_ATTRIBUTE'):
            tag = attr.findtext('TAG', '')
            value = attr.findtext('VALUE', '')
            if tag and value:
                sample_attributes[tag] = value
    
    return {
        "sra_id": sra_id,
        "exp_accession": exp_accession,
        "exp_title": exp_title,
        "study_accession": study_accession,
        "bioproject": bioproject,
        "organism": organism,
        "library_strategy": library_strategy,
        "library_source": library_source,
        "library_selection": library_selection,
        "library_name": library_name,
        "library_layout": library_layout,
        "instrument": instrument_name if instrument_name else ("_".join(instruments) if instruments else ""),
        "biosample": biosample,
        "sample_attributes": sample_attributes,
        "runs": runs,
        "runs_details": runs_details,  # Detailed run info with spots and size
        "run_count": len(runs),
        "fetched_at": datetime.now().isoformat()
    }


# ============================================
# SRA FETCHER CLASS
# ============================================

class SRAFetcher:
    """Main class for fetching SRA experiment data from NCBI."""
    
    def __init__(self, email: str = NCBI_EMAIL, api_key: str = NCBI_API_KEY):
        """
        Initialize SRA Fetcher.
        
        Args:
            email: Email for NCBI Entrez (required)
            api_key: Optional API key for higher rate limits
        """
        self.logger = setup_logging()
        
        # Configure Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            self.logger.info("✓ Using NCBI API key (10 req/sec)")
        else:
            self.logger.warning("⚠ No API key - rate limited to 3 req/sec")
        
        # Storage
        self.results: List[Dict] = []
        self.bioproject_experiments: Dict[str, List[str]] = {}  # Map BioProject → SRA IDs
        
        # Statistics
        self.stats = {
            "bioproject_searched": "",
            "total_sra_experiments": 0,
            "successfully_fetched": 0,
            "fetch_errors": 0,
            "unique_bioprojects": set()
        }
    
    def search_sra_by_bioproject(self, bioproject_id: str, retmax: int = 1000) -> List[str]:
        """
        Search SRA database for experiments linked to a BioProject.
        
        Args:
            bioproject_id: BioProject accession (e.g., "PRJNA1179470")
            retmax: Maximum results to retrieve
        
        Returns:
            List of SRA experiment IDs
        """
        self.logger.info(f"\n{'='*70}")
        self.logger.info(f"Searching SRA for BioProject: {bioproject_id}")
        self.logger.info(f"{'='*70}")
        
        self.stats["bioproject_searched"] = bioproject_id
        
        try:
            # CRITICAL: Do NOT use rettype="json" - use XML (default)
            handle = Entrez.esearch(
                db='sra',
                term=bioproject_id,
                retmax=retmax
            )
            search_result = Entrez.read(handle)
            handle.close()
            
            count = int(search_result.get('Count', 0))
            sra_ids = search_result.get('IdList', [])
            
            self.logger.info(f"✅ Found {count} SRA experiments")
            self.logger.info(f"   Retrieved {len(sra_ids)} IDs")
            
            self.stats["total_sra_experiments"] = count
            return sra_ids
            
        except Exception as e:
            self.logger.error(f"Search failed: {e}")
            return []
    
    def fetch_sra_metadata(self, sra_id: str) -> Optional[Dict]:
        """
        Fetch metadata for a single SRA experiment.
        
        CRITICAL: Use efetch with rettype="xml" instead of esummary.
        
        Args:
            sra_id: SRA experiment accession
        
        Returns:
            Metadata dictionary or None if failed
        """
        try:
            time.sleep(RATE_LIMIT)
            
            # FIXED: Use efetch with rettype="xml" (not esummary with rettype="json")
            handle = Entrez.efetch(
                db='sra',
                id=sra_id,
                rettype='xml'
            )
            xml_data = handle.read()
            handle.close()
            
            # Parse and extract metadata
            metadata = extract_sra_experiment_metadata(xml_data.decode('utf-8'), sra_id)
            
            if metadata:
                return metadata
            else:
                self.logger.warning(f"No valid data extracted from {sra_id}")
                return None
                
        except Exception as e:
            self.logger.error(f"Failed to fetch {sra_id}: {e}")
            self.stats["fetch_errors"] += 1
            return None
    
    def fetch_all_by_bioproject(self, bioproject_id: str, max_per_bioproject: Optional[int] = None) -> List[Dict]:
        """
        Main method: search and fetch all SRA experiments for a BioProject.
        
        Args:
            bioproject_id: BioProject accession (e.g., "PRJNA1179470")
            max_per_bioproject: Max experiments to fetch per BioProject (None = all)
        
        Returns:
            List of metadata dictionaries
        """
        # Reset results
        self.results = []
        self.stats["successfully_fetched"] = 0
        self.stats["fetch_errors"] = 0
        
        # Search SRA
        sra_ids = self.search_sra_by_bioproject(bioproject_id)
        
        if not sra_ids:
            self.logger.warning("No SRA experiments found")
            return []
        
        # Limit if specified
        if max_per_bioproject:
            sra_ids = sra_ids[:max_per_bioproject]
        
        # Fetch each experiment
        total = len(sra_ids)
        for idx, sra_id in enumerate(sra_ids, 1):
            self.logger.info(f"\n[{idx}/{total}] Fetching {sra_id}...")
            
            metadata = self.fetch_sra_metadata(sra_id)
            
            if metadata:
                title = metadata.get('exp_title', 'No title')
                self.logger.info(f"   ✓ {title[:80]}...")
                self.logger.info(f"   Strategy: {metadata.get('library_strategy', '?')}")
                self.logger.info(f"   Runs: {metadata.get('run_count', 0)}")
                
                self.results.append(metadata)
                self.stats["successfully_fetched"] += 1
                
                # Track BioProject
                if bioproject_id not in self.stats["unique_bioprojects"]:
                    self.stats["unique_bioprojects"].add(bioproject_id)
        
        self._print_summary()
        return self.results
    
    def save_results(self, output_file: Path):
        """Save results to JSON file."""
        output_file = Path(output_file)
        
        data = {
            "metadata": {
                "fetched_at": datetime.now().isoformat(),
                "bioproject": self.stats["bioproject_searched"],
                "total_experiments": self.stats["total_sra_experiments"],
                "successfully_fetched": self.stats["successfully_fetched"],
                "fetch_errors": self.stats["fetch_errors"],
                "results_count": len(self.results)
            },
            "experiments": self.results
        }
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.logger.info(f"\n✓ Results saved to: {output_file}")
    
    def save_results_csv(self, output_file: Path):
        """Save results to CSV file."""
        output_file = Path(output_file)
        
        if not self.results:
            self.logger.warning("No results to save")
            return
        
        columns = [
            "sra_id",
            "exp_accession",
            "title",
            "study_accession",
            "organism",
            "library_strategy",
            "library_source",
            "library_selection",
            "biosample",
            "run_count",
            "runs",
            "fetched_at"
        ]
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            
            for row in self.results:
                # Convert runs list to comma-separated string
                row_copy = row.copy()
                row_copy['runs'] = ','.join(row_copy.get('runs', []))
                writer.writerow(row_copy)
        
        self.logger.info(f"\n✓ Results saved to: {output_file}")
    
    def _print_summary(self):
        """Print summary statistics."""
        self.logger.info("\n" + "=" * 70)
        self.logger.info("FETCH SUMMARY")
        self.logger.info("=" * 70)
        self.logger.info(f"BioProject: {self.stats['bioproject_searched']}")
        self.logger.info(f"Total SRA experiments found: {self.stats['total_sra_experiments']}")
        self.logger.info(f"Successfully fetched: {self.stats['successfully_fetched']}")
        self.logger.info(f"Fetch errors: {self.stats['fetch_errors']}")
        self.logger.info(f"Results stored: {len(self.results)}")
        self.logger.info("=" * 70)


# ============================================
# CONVENIENCE FUNCTIONS
# ============================================

def fetch_sra_experiments(bioproject_id: str, output_file: Optional[Path] = None,
                          max_per_bioproject: Optional[int] = None) -> List[Dict]:
    """
    Convenience function to fetch SRA experiments for a BioProject.
    
    Args:
        bioproject_id: BioProject accession (e.g., "PRJNA1179470")
        output_file: Optional output file path
        max_per_bioproject: Max experiments to fetch
    
    Returns:
        List of metadata dictionaries
    """
    fetcher = SRAFetcher()
    results = fetcher.fetch_all_by_bioproject(bioproject_id, max_per_bioproject)
    
    if results and output_file:
        if Path(output_file).suffix.lower() == ".csv":
            fetcher.save_results_csv(output_file)
        else:
            fetcher.save_results(output_file)
    
    return results


if __name__ == "__main__":
    # Test with a known BioProject that has SRA data
    bioproject = "PRJNA1179470"
    print(f"\n🔍 Testing SRA Fetcher with: {bioproject}\n")
    
    fetcher = SRAFetcher()
    results = fetcher.fetch_all_by_bioproject(bioproject, max_per_bioproject=5)
    
    print(f"\n✅ Fetched {len(results)} experiments")
    
    if results:
        print("\n📊 Sample results:")
        for exp in results[:2]:
            print(f"\n  Experiment: {exp.get('exp_accession', '?')}")
            print(f"  Title: {exp.get('title', '?')[:60]}...")
            print(f"  Strategy: {exp.get('library_strategy', '?')}")
            print(f"  Organism: {exp.get('organism', '?')}")
            print(f"  Runs: {exp.get('run_count', 0)}")
