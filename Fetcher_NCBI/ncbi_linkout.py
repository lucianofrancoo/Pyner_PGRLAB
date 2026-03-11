"""
NCBI LinkOut: BioSample → PubMed Linking
=========================================
Links SRA experiments to PubMed/PMC publications via BioSample identifiers.

Strategy:
1. Extract BioSamples from SRA experiments
2. Search PubMed for BioSample accessions (SAMN...)
3. Search PubMed for SRA experiment accessions (SRX...)
4. Map BioProject → Publications

Usage:
    from Fetcher_NCBI.ncbi_linkout import LinkoutFetcher
    
    linkout = LinkoutFetcher()
    pubmed_hits = linkout.search_publications_for_biosamples(
        ['SAMN44494209', 'SAMN44494208', ...]
    )
"""

import time
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Set, Optional, Tuple
from Bio import Entrez
from config import NCBI_EMAIL, NCBI_API_KEY, RATE_LIMIT, LOGS_DIR


# ============================================
# LOGGING
# ============================================

def setup_logging(log_file: Optional[str] = None) -> logging.Logger:
    """Configure logging for linkout fetcher."""
    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = LOGS_DIR / f"linkout_{timestamp}.log"
    
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
# PUBMED PARSER
# ============================================

def extract_pubmed_metadata(pubmed_id: str) -> Optional[Dict]:
    """
    Fetch metadata for a PubMed record.
    
    Args:
        pubmed_id: PubMed ID
    
    Returns:
        Dictionary with: pmid, title, authors, year, doi, url
    """
    try:
        time.sleep(RATE_LIMIT)
        
        handle = Entrez.efetch(
            db='pubmed',
            id=pubmed_id,
            rettype='xml'
        )
        xml_data = handle.read()
        handle.close()
        
        # Parse XML (simplified - look for key fields)
        import xml.etree.ElementTree as ET
        root = ET.fromstring(xml_data)
        
        # Extract article info
        article = root.find('.//Article')
        if article is None:
            return None
        
        # Title (use itertext() to get all text including sub-elements like <i>, <sub>, etc.)
        title_elem = article.find('.//ArticleTitle')
        title = ''.join(title_elem.itertext()) if title_elem is not None else ''
        
        # Abstract (use itertext() to get all text including sub-elements like <sub>, <sup>, etc.)
        abstract_elem = article.find('.//AbstractText')
        abstract = ''.join(abstract_elem.itertext()) if abstract_elem is not None else ''
        
        # Year
        year_elem = article.find('.//PubDate/Year')
        year = year_elem.text if year_elem is not None else ''
        
        # DOI
        doi_elem = article.find('.//ELocationID[@EIdType="doi"]')
        doi = doi_elem.text if doi_elem is not None else ''
        
        # PMCID (PubMed Central ID)
        pmc_elem = root.find('.//ArticleId[@IdType="pmc"]')
        pmcid = pmc_elem.text if pmc_elem is not None else ''
        
        # Journal name
        journal_elem = article.find('.//Journal/Title')
        journal = journal_elem.text if journal_elem is not None else ''
        
        # Publication types (e.g., "Review", "Research Article", etc.)
        pub_types = []
        for pt in article.findall('.//PublicationType'):
            if pt.text:
                pub_types.append(pt.text)
        publication_type = "; ".join(pub_types) if pub_types else ''
        
        # Authors
        authors = []
        for author in article.findall('.//Author'):
            last_name = author.findtext('LastName', '')
            initials = author.findtext('Initials', '')
            if last_name:
                authors.append(f"{last_name} {initials}".strip())
        
        # Use the pubmed_id parameter directly (more reliable than parsing XML)
        url = f"https://www.ncbi.nlm.nih.gov/pubmed/{pubmed_id}"
        
        return {
            'pmid': pubmed_id,
            'title': title,
            'abstract': abstract,  # Full abstract without truncation
            'authors': authors[:3],  # First 3 authors
            'year': year,
            'journal': journal,
            'publication_type': publication_type,
            'doi': doi,
            'pmcid': pmcid,
            'url': url,
            'fetched_at': datetime.now().isoformat()
        }
        
    except Exception as e:
        logging.error(f"Failed to fetch PubMed {pubmed_id}: {e}")
        return None


# ============================================
# LINKOUT FETCHER
# ============================================

class LinkoutFetcher:
    """Links SRA data to PubMed publications."""
    
    def __init__(self, email: str = NCBI_EMAIL, api_key: str = NCBI_API_KEY):
        """Initialize linkout fetcher."""
        self.logger = setup_logging()
        
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key
            self.logger.info("✓ Using NCBI API key")
        
        self.results: Dict[str, Dict] = {}  # identifier → pubmed hits
        self.stats = {
            "queries": 0,
            "biosamples_queried": 0,
            "accessions_queried": 0,
            "pubmed_hits": 0,
            "unique_pmids": set()
        }
    
    def search_pubmed_for_identifier(self, identifier: str, search_field: str = "All Fields") -> List[str]:
        """
        Search PubMed for a specific identifier (SAMN, SRX, etc).
        
        Args:
            identifier: SAMN44494209, SRX26886995, etc.
            search_field: PubMed search field
        
        Returns:
            List of PubMed IDs
        """
        try:
            time.sleep(RATE_LIMIT)
            
            # Search with specific identifier
            query = f'"{identifier}"[{search_field}]'
            self.logger.info(f"Searching PubMed: {query}")
            
            handle = Entrez.esearch(
                db='pubmed',
                term=query,
                retmax=100
            )
            result = Entrez.read(handle)
            handle.close()
            
            pmids = result.get('IdList', [])
            
            if pmids:
                self.logger.info(f"  ✓ Found {len(pmids)} PubMed records for {identifier}")
            else:
                self.logger.info(f"  × No PubMed records for {identifier}")
            
            self.stats["queries"] += 1
            return pmids
            
        except Exception as e:
            self.logger.error(f"Search failed for {identifier}: {e}")
            return []
    
    def search_publications_for_biosamples(self, biosamples: List[str]) -> Dict[str, List[Dict]]:
        """
        Search PubMed for all BioSamples from an experiment.
        
        Args:
            biosamples: List of SAMN... accessions
        
        Returns:
            Dictionary: biosample → list of PubMed records
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("🔍 LINKING BIOSAMPLES TO PUBMED")
        self.logger.info("="*70)
        
        results = {}
        
        for biosample in biosamples:
            self.logger.info(f"\n[{biosamples.index(biosample)+1}/{len(biosamples)}] {biosample}")
            
            # Search for BioSample
            pmids = self.search_pubmed_for_identifier(biosample, search_field="All Fields")
            
            if pmids:
                # Fetch metadata for each hit
                publications = []
                for pmid in pmids[:5]:  # Limit to first 5 hits
                    pub = extract_pubmed_metadata(pmid)
                    if pub:
                        publications.append(pub)
                        self.stats["unique_pmids"].add(pmid)
                
                results[biosample] = publications
                self.stats["pubmed_hits"] += len(publications)
                self.stats["biosamples_queried"] += 1
        
        return results
    
    def search_publications_for_sra_accessions(self, sra_accessions: List[str]) -> Dict[str, List[Dict]]:
        """
        Search PubMed for SRA experiment accessions (SRX...).
        
        Args:
            sra_accessions: List of SRX... accessions
        
        Returns:
            Dictionary: accession → list of PubMed records
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("🔍 LINKING SRA ACCESSIONS TO PUBMED")
        self.logger.info("="*70)
        
        results = {}
        
        for acc in sra_accessions:
            self.logger.info(f"\n[{sra_accessions.index(acc)+1}/{len(sra_accessions)}] {acc}")
            
            # Search for SRA accession
            pmids = self.search_pubmed_for_identifier(acc, search_field="All Fields")
            
            if pmids:
                # Fetch metadata
                publications = []
                for pmid in pmids[:3]:  # Limit to first 3 hits
                    pub = extract_pubmed_metadata(pmid)
                    if pub:
                        publications.append(pub)
                        self.stats["unique_pmids"].add(pmid)
                
                results[acc] = publications
                self.stats["pubmed_hits"] += len(publications)
                self.stats["accessions_queried"] += 1
        
        return results
    
    def search_publications_for_bioproject(self, bioproject_id: str) -> Dict[str, List[Dict]]:
        """
        Search PubMed for a BioProject identifier directly.
        
        Args:
            bioproject_id: PRJNA... accession
        
        Returns:
            Dictionary with PubMed results
        """
        self.logger.info("\n" + "="*70)
        self.logger.info(f"🔍 SEARCHING PUBMED FOR BIOPROJECT: {bioproject_id}")
        self.logger.info("="*70)
        
        try:
            time.sleep(RATE_LIMIT)
            
            query = f'"{bioproject_id}"[All Fields]'
            self.logger.info(f"Query: {query}")
            
            handle = Entrez.esearch(db='pubmed', term=query, retmax=100)
            result = Entrez.read(handle)
            handle.close()
            
            pmids = result.get('IdList', [])
            self.logger.info(f"Found {len(pmids)} PubMed records")
            
            publications = []
            for pmid in pmids[:10]:  # First 10 hits
                pub = extract_pubmed_metadata(pmid)
                if pub:
                    publications.append(pub)
                    self.stats["unique_pmids"].add(pmid)
            
            return {
                "bioproject": bioproject_id,
                "publications": publications
            }
            
        except Exception as e:
            self.logger.error(f"Search failed: {e}")
            return {"bioproject": bioproject_id, "publications": []}
    
    def search_publications_by_boolean_query(
        self,
        query: str,
        max_results: int = 100,
        db: str = 'pubmed',
        progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
    ) -> List[Dict]:
        """
        Direct boolean search in PubMed or PMC.
        
        Args:
            query: Boolean query string (e.g., "Arabidopsis AND phosphate")
            max_results: Maximum number of publications to retrieve
            db: Database to search ('pubmed' or 'pmc')
        
        Returns:
            List of publication dictionaries with metadata
        """
        db_label = "PMC (full-text)" if db == 'pmc' else "PubMed"
        
        # PMC uses different field tags than PubMed (see PMC User Guide)
        # [All Fields] → [all], [Organism] not supported → [all]
        search_query = query
        if db == 'pmc':
            import re as re_mod
            search_query = re_mod.sub(r'\[Organism\]', '[all]', search_query)
            search_query = re_mod.sub(r'\[All Fields\]', '[all]', search_query)
            if search_query != query:
                self.logger.info("ℹ️  Adjusted query for PMC: [Organism]/[All Fields] → [all]")
        
        self.logger.info("\n" + "="*70)
        self.logger.info(f"🔍 DIRECT {db_label.upper()} SEARCH")
        self.logger.info("="*70)
        self.logger.info(f"Query: {search_query}")
        self.logger.info(f"Database: {db_label}")
        self.logger.info(f"Max results: {max_results}")
        
        publications = []
        
        try:
            time.sleep(RATE_LIMIT)
            
            # Search the specified database
            handle = Entrez.esearch(
                db=db,
                term=search_query,
                retmax=max_results,
                sort='relevance'
            )
            result = Entrez.read(handle)
            handle.close()
            
            ids = result.get('IdList', [])
            total_count = int(result.get('Count', 0))
            
            self.logger.info(f"Found {total_count} total matches in {db_label}")
            self.logger.info(f"Retrieving {len(ids)} publications...")
            if progress_callback:
                progress_callback(
                    {
                        "stage": "searching",
                        "processed": 0,
                        "target": len(ids) if ids else max_results,
                        "message": f"Found {total_count} matches in {db_label}.",
                    }
                )
            
            if db == 'pmc':
                # PMC returns PMC IDs — convert to PMIDs via elink
                if progress_callback:
                    progress_callback(
                        {
                            "stage": "converting",
                            "processed": 0,
                            "target": len(ids) if ids else max_results,
                            "message": "Converting PMC IDs to PMID...",
                        }
                    )
                pmids = self._convert_pmc_to_pmids(ids, progress_callback=progress_callback)
                self.logger.info(f"Converted {len(pmids)} PMC articles to PubMed IDs")
            else:
                pmids = ids
            
            # Fetch metadata for each PMID
            for i, pmid in enumerate(pmids, 1):
                self.logger.info(f"  [{i}/{len(pmids)}] Fetching PMID:{pmid}")
                if progress_callback:
                    progress_callback(
                        {
                            "stage": "fetching",
                            "processed": i - 1,
                            "target": len(pmids) if pmids else max_results,
                            "message": f"Fetching PMID metadata ({i}/{len(pmids)})...",
                        }
                    )
                pub = extract_pubmed_metadata(pmid)
                if pub:
                    publications.append(pub)
                    self.stats["unique_pmids"].add(pmid)
                    self.stats["pubmed_hits"] += 1
                if progress_callback:
                    progress_callback(
                        {
                            "stage": "fetching",
                            "processed": i,
                            "target": len(pmids) if pmids else max_results,
                            "message": f"Fetched PMID metadata ({i}/{len(pmids)})...",
                        }
                    )
            
            self.stats["queries"] += 1
            
            self.logger.info(f"\n✓ Retrieved {len(publications)} publications")
            
        except Exception as e:
            self.logger.error(f"{db_label} search failed: {e}")
        
        return publications
    
    def _convert_pmc_to_pmids(
        self,
        pmc_ids: List[str],
        progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
    ) -> List[str]:
        """
        Convert PMC IDs to PubMed IDs using Entrez elink.
        Only gets the article's own PMID (pmc_pubmed), not referenced articles.
        
        Args:
            pmc_ids: List of PMC IDs (numeric, without 'PMC' prefix)
        
        Returns:
            List of PubMed IDs
        """
        if not pmc_ids:
            return []
        
        pmids = []
        # Process in smaller batches to avoid IncompleteRead errors
        batch_size = 20
        total_batches = (len(pmc_ids) + batch_size - 1) // batch_size
        
        for batch_num, i in enumerate(range(0, len(pmc_ids), batch_size), 1):
            batch = pmc_ids[i:i + batch_size]
            self.logger.info(f"  Converting PMC batch {batch_num}/{total_batches} ({len(batch)} IDs)...")
            if progress_callback:
                progress_callback(
                    {
                        "stage": "converting",
                        "processed": min(i, len(pmc_ids)),
                        "target": len(pmc_ids),
                        "message": f"Converting PMC IDs to PMID (batch {batch_num}/{total_batches})...",
                    }
                )
            
            # Retry logic for network reliability
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    time.sleep(RATE_LIMIT * (attempt + 1))  # Increasing delay on retries
                    handle = Entrez.elink(
                        dbfrom='pmc',
                        db='pubmed',
                        linkname='pmc_pubmed',  # Only the article's own PMID, NOT references
                        id=batch
                    )
                    results = Entrez.read(handle)
                    handle.close()
                    
                    for record in results:
                        link_sets = record.get('LinkSetDb', [])
                        for link_set in link_sets:
                            for link in link_set.get('Link', []):
                                pmid = link.get('Id', '')
                                if pmid and pmid not in pmids:
                                    pmids.append(pmid)
                    if progress_callback:
                        progress_callback(
                            {
                                "stage": "converting",
                                "processed": min(i + len(batch), len(pmc_ids)),
                                "target": len(pmc_ids),
                                "message": f"Converted PMC batch {batch_num}/{total_batches}.",
                            }
                        )
                    break  # Success, exit retry loop
                    
                except Exception as e:
                    if attempt < max_retries - 1:
                        self.logger.warning(f"⚠️ Batch {batch_num} attempt {attempt+1} failed: {e}. Retrying...")
                    else:
                        self.logger.warning(f"⚠️ Batch {batch_num} failed after {max_retries} attempts: {e}")
        
        return pmids
    
    def print_summary(self):
        """Print summary statistics."""
        self.logger.info("\n" + "="*70)
        self.logger.info("LINKOUT SUMMARY")
        self.logger.info("="*70)
        self.logger.info(f"Queries executed: {self.stats['queries']}")
        self.logger.info(f"BioSamples queried: {self.stats['biosamples_queried']}")
        self.logger.info(f"SRA accessions queried: {self.stats['accessions_queried']}")
        self.logger.info(f"PubMed hits: {self.stats['pubmed_hits']}")
        self.logger.info(f"Unique PMIDs: {len(self.stats['unique_pmids'])}")
        self.logger.info("="*70)
    
    def save_results(self, results: Dict, output_file: Path):
        """Save linkout results to JSON."""
        output_file = Path(output_file)
        
        data = {
            "metadata": {
                "fetched_at": datetime.now().isoformat(),
                "biosamples_queried": self.stats["biosamples_queried"],
                "accessions_queried": self.stats["accessions_queried"],
                "pubmed_hits": self.stats["pubmed_hits"],
                "unique_pmids": len(self.stats["unique_pmids"])
            },
            "linkout_results": results
        }
        
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        self.logger.info(f"\n✓ Results saved to: {output_file}")


# ============================================
# CONVENIENCE FUNCTIONS
# ============================================

def link_sra_to_pubmed(sra_data: Dict) -> Dict:
    """
    Link SRA experiments to PubMed publications.
    
    Args:
        sra_data: SRA JSON data with experiments
    
    Returns:
        Enriched data with PubMed links
    """
    linkout = LinkoutFetcher()
    
    # Extract unique BioSamples and accessions
    biosamples = set()
    accessions = set()
    
    for exp in sra_data.get('experiments', []):
        biosample = exp.get('biosample', '')
        if biosample:
            biosamples.add(biosample)
        
        exp_acc = exp.get('exp_accession', '')
        if exp_acc:
            accessions.add(exp_acc)
    
    biosamples = list(biosamples)[:20]  # Limit to first 20 to avoid rate limiting
    accessions = list(accessions)[:10]
    
    # Search
    pubmed_biosample_results = linkout.search_publications_for_biosamples(biosamples)
    pubmed_accession_results = linkout.search_publications_for_sra_accessions(accessions)
    
    # Combine results
    enriched_data = sra_data.copy()
    enriched_data['pubmed_links'] = {
        'via_biosamples': pubmed_biosample_results,
        'via_accessions': pubmed_accession_results,
        'unique_pmids': list(linkout.stats["unique_pmids"])
    }
    
    linkout.print_summary()
    
    return enriched_data


if __name__ == "__main__":
    # Test: Link SRA data to PubMed
    print("🔗 SRA to PubMed Linkout Test")
    print("="*70)
    
    linkout = LinkoutFetcher()
    
    # Test with sample BioProject
    print("\n1️⃣ Direct BioProject search:")
    bp_results = linkout.search_publications_for_bioproject('PRJNA1179470')
    
    if bp_results['publications']:
        print(f"\n Found {len(bp_results['publications'])} publications!")
        for pub in bp_results['publications'][:2]:
            print(f"\n  Title: {pub.get('title', '')[:60]}...")
            print(f"  Authors: {', '.join(pub.get('authors', []))}")
            print(f"  Year: {pub.get('year', '')}")
            print(f"  URL: {pub.get('url', '')}")
    else:
        print(" No publications found via BioProject ID")
    
    # Test with sample BioSamples
    print("\n\n2️⃣ BioSample search:")
    sample_biosamples = ['SAMN44494209', 'SAMN44494208', 'SAMN44494207']
    biosample_results = linkout.search_publications_for_biosamples(sample_biosamples)
    
    total_hits = sum(len(v) for v in biosample_results.values())
    print(f"\n Total hits: {total_hits}")
    
    linkout.print_summary()
