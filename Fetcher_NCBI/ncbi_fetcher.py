"""
NCBI Fetcher Core Module
=========================
Handles querying NCBI databases using Entrez API.
Includes deduplication, metadata extraction, and result storage.
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
    NCBI_EMAIL, NCBI_API_KEY, MAX_RESULTS, DATABASE,
    RATE_LIMIT, DEFAULT_OUTPUT, DEDUP_CACHE, LOGS_DIR,
    MIN_UNIQUE_BIOPROJECTS
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
    Tracks processed BioProjects to avoid duplicates WITHIN CURRENT SEARCH ONLY.
    Does not persist between searches.
    """
    def __init__(self, cache_file: Path = DEDUP_CACHE):
        self.cache_file = cache_file
        # Start fresh for each search - don't load previous cache
        self.seen: Set[str] = set()
    
    def _load_cache(self) -> Set[str]:
        """Load existing cache from disk."""
        # Disabled - only track duplicates in current search
        return set()
    
    def save(self):
        """Save cache to disk. Disabled for per-search deduplication."""
        # Don't persist cache between searches
        pass
    
    def is_seen(self, bioproject: str) -> bool:
        """Check if BioProject has been processed in current search."""
        return bioproject in self.seen
    
    def add(self, bioproject: str):
        """Mark BioProject as processed in current search."""
        self.seen.add(bioproject)
    
    def clear(self):
        """Clear the cache."""
        self.seen.clear()


# ============================================
# METADATA EXTRACTION
# ============================================

def extract_sra_metadata(xml_string: str) -> Optional[Dict]:
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
    
    return {
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
def extract_bioproject_metadata(summary: Dict) -> Optional[Dict]:
    """
    Extract metadata from BioProject summary response.
    
    NOTE: NCBI BioProject esummary has limited metadata available.
    These fields are populated when available:
    - Organism fields: Organism_Name, Organism_Label, Organism_Strain
    - Date: Registration_Date (no Public_Date or Submission_Date available)
    
    The following are NOT available in BioProject esummary:
    - Total_Studies, Total_Runs (need SRA query)
    - Publications, Experimental_Design (not in API response)
    """
    try:
        bioproject = summary.get("Project_Acc", "")
        title = summary.get("Project_Title", "")
        description = summary.get("Project_Description", "")
        
        # Organismos - NCBI BioProject has individual organism fields, not a list
        organism_name = summary.get("Organism_Name", "")
        organism_label = summary.get("Organism_Label", "")
        organisms = organism_name or organism_label or ""
        
        # Strain (available per BioProject)
        strain = summary.get("Organism_Strain", "")
        
        # Cultivar - NOT available in BioProject esummary
        cultivar = ""
        
        # Dates - only Registration_Date is available, not Submission/Public dates
        registration_date = summary.get("Registration_Date", "")
        
        # Conteos - NOT available in BioProject esummary (would need SRA query)
        # These would require separate SRA database query
        total_studies = ""
        total_runs = ""
        
        # Publicaciones - NOT available in BioProject esummary
        publications = ""
        
        # Diseño experimental - NOT available in BioProject esummary
        experimental_design = ""
        
        if not bioproject:
            return None
        
        return {
            "bioproject": bioproject,
            "title": title,
            "description": description,
            "organisms": organisms,
            "strain": strain,
            "cultivar": cultivar,
            "submission_date": registration_date,  # Only available date from NCBI
            "public_date": "",  # Not available in BioProject esummary
            "total_studies": total_studies,  # Empty - not in BioProject esummary
            "total_runs": total_runs,  # Empty - not in BioProject esummary
            "publications": publications,  # Empty - not in BioProject esummary
            "experimental_design": experimental_design  # Empty - not in BioProject esummary
        }
    except Exception as e:
        print(f"Error extracting metadata: {e}")
        return None


# ============================================
# NCBI FETCHER CLASS
# ============================================

class NCBIFetcher:
    """
    Main class for fetching data from NCBI databases.
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
    
    def search(self, query: str, max_results: Optional[int] = MAX_RESULTS, batch_size: int = 500) -> List[str]:
        """
        Search NCBI database with boolean query.

        Args:
            query: Boolean search query
            max_results: Maximum results to retrieve

        Returns:
            List of record IDs
        """
        self.logger.info(f"Searching NCBI {DATABASE}: {query}")
        if max_results is None:
            self.logger.info("Max results: all")
        else:
            self.logger.info(f"Max results: {max_results}")
        
        try:
            # First request: get total count
            handle = Entrez.esearch(
                db=DATABASE,
                term=query,
                retmax=0
            )
            search_results = Entrez.read(handle)
            handle.close()

            total_count = int(search_results.get("Count", 0))
            self.stats["total_found"] = total_count

            if total_count == 0:
                self.logger.info("No records found")
                return []

            target_count = total_count if max_results is None else min(max_results, total_count)
            self.logger.info(f"Found {total_count} total records")
            self.logger.info(f"Retrieving {target_count} IDs to process")

            id_list: List[str] = []
            retstart = 0
            while retstart < target_count:
                retmax = min(batch_size, target_count - retstart)

                handle = Entrez.esearch(
                    db=DATABASE,
                    term=query,
                    retstart=retstart,
                    retmax=retmax
                )
                page_results = Entrez.read(handle)
                handle.close()

                id_list.extend(page_results.get("IdList", []))
                retstart += retmax

                self.logger.info(f"Retrieved {len(id_list)}/{target_count} IDs")
                time.sleep(RATE_LIMIT)

            return id_list
            
        except Exception as e:
            self.logger.error(f"Search failed: {e}")
            return []
    
    def fetch_metadata(self, study_id: str) -> Optional[Dict]:
        """
        Fetch metadata for a single record.
        
        Args:
                study_id: NCBI record ID
        
        Returns:
            Metadata dictionary or None if failed
        """
        try:
            # Add rate limiting
            time.sleep(RATE_LIMIT)
            
            # Fetch summary
            handle = Entrez.esummary(db=DATABASE, id=study_id, rettype="docsum")
            summary = Entrez.read(handle)
            handle.close()

            if DATABASE == "bioproject":
                docset = summary.get("DocumentSummarySet", {})
                docs = docset.get("DocumentSummary", [])
                if not docs:
                    self.logger.warning(f"Record {study_id} has no DocumentSummary - skipping")
                    return None
                metadata = extract_bioproject_metadata(docs[0])
            else:
                study = summary[0] if isinstance(summary, list) else summary
                if "ExpXml" not in study:
                    self.logger.warning(f"Record {study_id} has no ExpXml - skipping")
                    return None
                metadata = extract_sra_metadata(study["ExpXml"])
            
            if metadata:
                # Add study ID
                metadata["sra_id"] = study_id
                # Add fetch timestamp
                metadata["fetched_at"] = datetime.now().isoformat()
            
            return metadata
            
        except Exception as e:
            self.logger.error(f"Failed to fetch record {study_id}: {e}")
            return None
    
    def fetch_all(self, query: str, max_results: Optional[int] = MAX_RESULTS, 
                  deduplicate: bool = True, min_unique: int = MIN_UNIQUE_BIOPROJECTS) -> List[Dict]:
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
        
        # Step 2: Process each record
        total_ids = len(id_list)
        for idx, study_id in enumerate(id_list, 1):
            self.logger.info(f"\n[{idx}/{total_ids}] Processing record {study_id}")

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
            if DATABASE == "sra":
                self.logger.info(f"  Strategy: {metadata.get('library_strategy', '')}")

            self.results.append(metadata)
            self.cache.add(bioproject)
            self.stats["processed"] += 1

            if min_unique and len(self.results) >= min_unique:
                self.logger.info(f"Reached {min_unique} unique BioProjects. Stopping early.")
                break
        
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

    def save_results_csv(self, output_file: Path):
        """
        Save results to CSV file.

        Args:
            output_file: Path to output file
        """
        output_file = Path(output_file)

        if not self.results:
            self.logger.warning("No results to save")
            return

        if DATABASE == "bioproject":
            columns = [
                "bioproject",
                "title",
                "submission_date",
                "description",
                "fetched_at"
            ]
        else:
            columns = [
                "sra_id",
                "bioproject",
                "title",
                "study_name",
                "organisms",
                "library_strategy",
                "library_source",
                "library_selection",
                "library_layout",
                "library_name",
                "biosample",
                "fetched_at"
            ]

        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(self.results)

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
                     max_results: Optional[int] = MAX_RESULTS,
                     min_unique: int = MIN_UNIQUE_BIOPROJECTS) -> List[Dict]:
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
    results = fetcher.fetch_all(query, max_results=max_results, min_unique=min_unique)
    
    if results and output_file:
        if Path(output_file).suffix.lower() == ".csv":
            fetcher.save_results_csv(output_file)
        else:
            fetcher.save_results(output_file)
    
    return results


if __name__ == "__main__":
    # Test with simple query
    test_query = 'arabidopsis[Organism] AND drought'
    print(f"Testing fetcher with query: {test_query}")
    results = fetch_from_query(test_query, max_results=10)
    print(f"\nFetched {len(results)} unique results")
