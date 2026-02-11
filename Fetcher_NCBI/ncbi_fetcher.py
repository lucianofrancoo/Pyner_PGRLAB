"""
NCBI Fetcher Core Module
=========================
Handles querying NCBI SRA database using Entrez API.
Includes deduplication, metadata extraction, and result storage.
"""

import time
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Optional
import xml.etree.ElementTree as ET
from Bio import Entrez

from config import (
    NCBI_EMAIL, NCBI_API_KEY, MAX_RESULTS, DATABASE,
    RATE_LIMIT, DEFAULT_OUTPUT, DEDUP_CACHE, LOGS_DIR
)


# ============================================
# LOGGING SETUP
# ============================================

def setup_logging(log_file: Optional[str] = None) -> logging.Logger:
    """
    Configure logging for the fetcher.
    """
    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = LOGS_DIR / f"fetcher_{timestamp}.log"
    
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
# DEDUPLICATION CACHE
# ============================================

class BioProjectCache:
    """
    Tracks processed BioProjects to avoid duplicates.
    """
    def __init__(self, cache_file: Path = DEDUP_CACHE):
        self.cache_file = cache_file
        self.seen: Set[str] = self._load_cache()
    
    def _load_cache(self) -> Set[str]:
        """Load existing cache from disk."""
        if self.cache_file.exists():
            try:
                with open(self.cache_file, 'r') as f:
                    data = json.load(f)
                    return set(data.get("bioprojects", []))
            except Exception:
                return set()
        return set()
    
    def save(self):
        """Save cache to disk."""
        with open(self.cache_file, 'w') as f:
            json.dump({
                "bioprojects": list(self.seen),
                "last_updated": datetime.now().isoformat()
            }, f, indent=2)
    
    def is_seen(self, bioproject: str) -> bool:
        """Check if BioProject has been processed."""
        return bioproject in self.seen
    
    def add(self, bioproject: str):
        """Mark BioProject as processed."""
        self.seen.add(bioproject)
    
    def clear(self):
        """Clear the cache."""
        self.seen.clear()
        self.save()


# ============================================
# METADATA EXTRACTION
# ============================================

def extract_metadata(xml_string: str) -> Optional[Dict]:
    """
    Extract metadata from SRA XML response.
    
    Args:
        xml_string: ExpXml string from Entrez.esummary
    
    Returns:
        Dictionary with extracted metadata or None if parsing fails
    """
    # Wrap XML in root element for valid parsing
    wrapped_xml = f"<Root>{xml_string}</Root>"
    
    try:
        root = ET.fromstring(wrapped_xml)
    except ET.ParseError as e:
        logging.error(f"XML parsing error: {e}")
        return None
    
    # Extract BioProject (unique identifier for biological study)
    bioproject = root.findtext(".//Bioproject", "")
    if not bioproject:
        return None  # Skip if no BioProject
    
    # Extract study information
    title = root.findtext(".//Title", "")
    study_elem = root.find(".//Study")
    study_name = study_elem.attrib.get("name", "") if study_elem is not None else ""
    
    # Extract organism
    organism_elem = root.find(".//Organism")
    organism = organism_elem.attrib.get("ScientificName", "") if organism_elem is not None else ""
    
    # Extract library information
    library_strategy = root.findtext(".//Library_descriptor/LIBRARY_STRATEGY", "")
    library_source = root.findtext(".//Library_descriptor/LIBRARY_SOURCE", "")
    library_selection = root.findtext(".//Library_descriptor/LIBRARY_SELECTION", "")
    
    # Determine library layout (paired vs single-end)
    library_layout = (
        "PAIRED" 
        if root.find(".//Library_descriptor/LIBRARY_LAYOUT/PAIRED") is not None 
        else "SINGLE"
    )
    
    # Extract tissue/sample information
    library_name = root.findtext(".//LIBRARY_NAME", "")
    biosample = root.findtext(".//Biosample", "")
    
    metadata = {
        "bioproject": bioproject,
        "title": title,
        "study_name": study_name,
        "organism": organism,
        "library_strategy": library_strategy,
        "library_source": library_source,
        "library_selection": library_selection,
        "library_layout": library_layout,
        "library_name": library_name,
        "biosample": biosample
    }
    
    return metadata


# ============================================
# NCBI FETCHER CLASS
# ============================================

class NCBIFetcher:
    """
    Main class for fetching data from NCBI SRA.
    """
    
    def __init__(self, email: str = NCBI_EMAIL, api_key: str = NCBI_API_KEY):
        """
        Initialize NCBI Fetcher.
        
        Args:
            email: Email for NCBI Entrez (required)
            api_key: Optional API key for higher rate limits
        """
        self.logger = setup_logging()
        
        # Configure Entrez
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            self.logger.info("Using NCBI API key (10 req/sec)")
        else:
            self.logger.warning("No API key - rate limited to 3 req/sec")
        
        # Initialize deduplication cache
        self.cache = BioProjectCache()
        
        # Storage for results
        self.results: List[Dict] = []
        
        # Statistics
        self.stats = {
            "total_found": 0,
            "processed": 0,
            "skipped_duplicate": 0,
            "skipped_error": 0,
            "unique_bioprojects": 0
        }
    
    def search(self, query: str, max_results: int = MAX_RESULTS) -> List[str]:
        """
        Search NCBI SRA with boolean query.
        
        Args:
            query: Boolean search query
            max_results: Maximum results to retrieve
        
        Returns:
            List of study IDs
        """
        self.logger.info(f"Searching NCBI SRA: {query}")
        self.logger.info(f"Max results: {max_results}")
        
        try:
            handle = Entrez.esearch(
                db=DATABASE,
                term=query,
                retmax=max_results
            )
            search_results = Entrez.read(handle)
            handle.close()
            
            total_count = int(search_results.get("Count", 0))
            id_list = search_results.get("IdList", [])
            
            self.stats["total_found"] = total_count
            
            self.logger.info(f"Found {total_count} total studies")
            self.logger.info(f"Retrieved {len(id_list)} IDs to process")
            
            return id_list
            
        except Exception as e:
            self.logger.error(f"Search failed: {e}")
            return []
    
    def fetch_metadata(self, study_id: str) -> Optional[Dict]:
        """
        Fetch metadata for a single study.
        
        Args:
            study_id: NCBI SRA study ID
        
        Returns:
            Metadata dictionary or None if failed
        """
        try:
            # Add rate limiting
            time.sleep(RATE_LIMIT)
            
            # Fetch summary
            handle = Entrez.esummary(db=DATABASE, id=study_id, rettype="docsum")
            study = Entrez.read(handle)[0]
            handle.close()
            
            # Check for ExpXml field
            if "ExpXml" not in study:
                self.logger.warning(f"Study {study_id} has no ExpXml - skipping")
                return None
            
            # Extract metadata from XML
            metadata = extract_metadata(study["ExpXml"])
            
            if metadata:
                # Add study ID
                metadata["sra_id"] = study_id
                # Add fetch timestamp
                metadata["fetched_at"] = datetime.now().isoformat()
            
            return metadata
            
        except Exception as e:
            self.logger.error(f"Failed to fetch study {study_id}: {e}")
            return None
    
    def fetch_all(self, query: str, max_results: int = MAX_RESULTS, 
                  deduplicate: bool = True) -> List[Dict]:
        """
        Main method: search and fetch all results with deduplication.
        
        Args:
            query: Boolean search query
            max_results: Maximum results to fetch
            deduplicate: Whether to skip duplicate BioProjects
        
        Returns:
            List of metadata dictionaries
        """
        self.logger.info("=" * 70)
        self.logger.info("NCBI FETCHER - Starting fetch operation")
        self.logger.info("=" * 70)
        
        # Reset results and stats
        self.results = []
        self.stats["processed"] = 0
        self.stats["skipped_duplicate"] = 0
        self.stats["skipped_error"] = 0
        
        # Step 1: Search
        id_list = self.search(query, max_results)
        
        if not id_list:
            self.logger.warning("No results found or search failed")
            return []
        
        # Step 2: Process each study
        for idx, study_id in enumerate(id_list, 1):
            self.logger.info(f"\n[{idx}/{len(id_list)}] Processing study {study_id}")
            
            metadata = self.fetch_metadata(study_id)
            
            if not metadata:
                self.stats["skipped_error"] += 1
                continue
            
            bioproject = metadata.get("bioproject", "")
            
            # Check deduplication
            if deduplicate and self.cache.is_seen(bioproject):
                self.logger.info(f"BioProject {bioproject} already seen - skipping")
                self.stats["skipped_duplicate"] += 1
                continue
            
            # New BioProject - add to results
            self.logger.info(f"✓ New BioProject: {bioproject}")
            self.logger.info(f"  Title: {metadata.get('title', '')[:60]}...")
            self.logger.info(f"  Organism: {metadata.get('organism', '')}")
            self.logger.info(f"  Strategy: {metadata.get('library_strategy', '')}")
            
            self.results.append(metadata)
            self.cache.add(bioproject)
            self.stats["processed"] += 1
        
        self.stats["unique_bioprojects"] = len(self.results)
        
        # Save cache
        self.cache.save()
        
        self._print_summary()
        
        return self.results
    
    def save_results(self, output_file: Path = DEFAULT_OUTPUT):
        """
        Save results to JSON file.
        
        Args:
            output_file: Path to output file
        """
        output_file = Path(output_file)
        
        data = {
            "metadata": {
                "fetched_at": datetime.now().isoformat(),
                "total_results": len(self.results),
                "statistics": self.stats
            },
            "results": self.results
        }
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.logger.info(f"\n✓ Results saved to: {output_file}")
        self.logger.info(f"  Total records: {len(self.results)}")
    
    def _print_summary(self):
        """Print summary statistics."""
        self.logger.info("\n" + "=" * 70)
        self.logger.info("FETCH SUMMARY")
        self.logger.info("=" * 70)
        self.logger.info(f"Total found in database: {self.stats['total_found']}")
        self.logger.info(f"Studies processed successfully: {self.stats['processed']}")
        self.logger.info(f"Skipped (duplicate BioProject): {self.stats['skipped_duplicate']}")
        self.logger.info(f"Skipped (errors): {self.stats['skipped_error']}")
        self.logger.info(f"Unique BioProjects collected: {self.stats['unique_bioprojects']}")
        self.logger.info("=" * 70)


# ============================================
# CONVENIENCE FUNCTIONS
# ============================================

def fetch_from_query(query: str, output_file: Optional[Path] = None, 
                     max_results: int = MAX_RESULTS) -> List[Dict]:
    """
    Convenience function to fetch data with one call.
    
    Args:
        query: Boolean search query
        output_file: Optional output file path
        max_results: Maximum results to fetch
    
    Returns:
        List of metadata dictionaries
    """
    fetcher = NCBIFetcher()
    results = fetcher.fetch_all(query, max_results=max_results)
    
    if results and output_file:
        fetcher.save_results(output_file)
    
    return results


if __name__ == "__main__":
    # Test with simple query
    test_query = 'arabidopsis[Organism] AND drought'
    print(f"Testing fetcher with query: {test_query}")
    results = fetch_from_query(test_query, max_results=10)
    print(f"\nFetched {len(results)} unique results")
