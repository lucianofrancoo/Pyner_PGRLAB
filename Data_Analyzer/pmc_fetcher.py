"""
PMC Full Text Fetcher
====================
Fetches full text from PubMed Central for papers with PMCID
"""

import logging
import re
from typing import Dict, Optional
from Bio import Entrez

logger = logging.getLogger(__name__)

# Configure Entrez (use same email as Fetcher_NCBI)
try:
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent.parent / "Fetcher_NCBI"))
    from config import NCBI_EMAIL, NCBI_API_KEY
    Entrez.email = NCBI_EMAIL
    if NCBI_API_KEY:
        Entrez.api_key = NCBI_API_KEY
except:
    Entrez.email = "user@example.com"
    logger.warning("Could not load NCBI credentials from Fetcher_NCBI/config.py")


class PMCFullTextFetcher:
    """Fetches full text from PubMed Central"""
    
    def __init__(self):
        self.max_text_length = 20000  # Characters to extract (avoid overwhelming LLM)
        self.max_methods_length = 8000  # More space for Methods section (has technique details)
    
    def fetch_full_text(self, pmcid: str) -> Optional[Dict[str, str]]:
        """
        Fetch full text sections from PMC
        
        Args:
            pmcid: PubMed Central ID (e.g., "PMC12513158" or "12513158")
            
        Returns:
            Dictionary with sections: {
                'abstract': str,
                'methods': str,
                'results': str,
                'full_text_preview': str
            }
            Returns None if fetch fails
        """
        # Clean PMCID
        pmcid = pmcid.strip().upper()
        if not pmcid.startswith('PMC'):
            pmcid = f'PMC{pmcid}'
        
        logger.info(f"Fetching full text from PMC: {pmcid}")
        
        try:
            # Fetch XML from PMC
            handle = Entrez.efetch(db="pmc", id=pmcid, rettype="xml", retmode="xml")
            xml_content = handle.read()
            handle.close()
                        # Decode bytes to string if necessary
            if isinstance(xml_content, bytes):
                xml_content = xml_content.decode('utf-8')
                        # Parse sections
            sections = self._parse_pmc_xml(xml_content)
            
            if sections:
                logger.info(f"✓ Successfully fetched full text for {pmcid}")
                return sections
            else:
                logger.warning(f"Could not parse sections from {pmcid}")
                return None
                
        except Exception as e:
            logger.warning(f"Failed to fetch full text for {pmcid}: {e}")
            return None
    
    def _parse_pmc_xml(self, xml_content: str) -> Optional[Dict[str, str]]:
        """
        Parse PMC XML and extract relevant sections
        
        Focuses on Methods and Results sections which contain most metadata
        """
        try:
            # Simple regex-based extraction (faster than full XML parsing)
            sections = {}
            
            # Extract abstract
            abstract_match = re.search(r'<abstract[^>]*>(.*?)</abstract>', xml_content, re.DOTALL | re.IGNORECASE)
            if abstract_match:
                abstract_text = self._clean_xml(abstract_match.group(1))
                sections['abstract'] = abstract_text[:2000]  # Limit abstract
            
            # Extract Methods section (contains organisms, tissues, conditions, techniques)
            # This section is important so we give it more space
            methods_patterns = [
                r'<sec[^>]*>.*?<title[^>]*>.*?(?:methods?|materials?|experimental|procedures?)[^<]*</title>(.*?)(?=</sec>|<sec)',
                r'<sec[^>]*sec-type=["\']methods["\'][^>]*>(.*?)(?=</sec>|<sec)',
                r'<sec[^>]*>.*?<title[^>]*>.*?(?:methods?|materials?)[^<]*</title>(.*?)(?=<sec|</sec>)'
            ]
            
            methods_text = ""
            for pattern in methods_patterns:
                match = re.search(pattern, xml_content, re.DOTALL | re.IGNORECASE)
                if match:
                    methods_text += self._clean_xml(match.group(1)) + " "
                    if len(methods_text) > self.max_methods_length:
                        break
            
            if methods_text:
                # Keep full Methods content - it has technique details
                sections['methods'] = methods_text[:self.max_methods_length]
            
            # Extract Results section (contains experimental details)
            results_patterns = [
                r'<sec[^>]*>.*?<title[^>]*>.*?(?:results?)[^<]*</title>(.*?)</sec>',
                r'<sec[^>]*sec-type=["\']results["\'][^>]*>(.*?)</sec>'
            ]
            
            results_text = ""
            for pattern in results_patterns:
                matches = re.finditer(pattern, xml_content, re.DOTALL | re.IGNORECASE)
                for match in matches:
                    results_text += self._clean_xml(match.group(1)) + " "
            
            if results_text:
                sections['results'] = results_text[:5000]  # Limit results
            
            # Create a preview: prioritize Methods (has technique details) then Results
            full_text_preview = ""
            if 'methods' in sections:
                # Give Methods more space since it has experimental technique details
                full_text_preview += sections['methods'][:9000] + " "
            if 'results' in sections:
                full_text_preview += sections['results'][:8000]
            
            if full_text_preview:
                sections['full_text_preview'] = full_text_preview[:self.max_text_length]
            
            return sections if sections else None
            
        except Exception as e:
            logger.error(f"Error parsing PMC XML: {e}")
            return None
    
    def _clean_xml(self, text: str) -> str:
        """Remove XML tags and clean up text"""
        # Remove all XML tags
        text = re.sub(r'<[^>]+>', '', text)
        
        # Decode common HTML entities
        text = text.replace('&lt;', '<')
        text = text.replace('&gt;', '>')
        text = text.replace('&amp;', '&')
        text = text.replace('&quot;', '"')
        text = text.replace('&#x2009;', ' ')
        text = text.replace('&#x00A0;', ' ')
        
        # Clean up whitespace
        text = re.sub(r'\s+', ' ', text)
        text = text.strip()
        
        return text


def test_pmc_fetcher():
    """Test PMC fetcher with a known PMCID"""
    fetcher = PMCFullTextFetcher()
    
    # Test with a known PMCID from your results
    test_pmcid = "PMC12513158"
    
    print(f"\nTesting PMC Full Text Fetcher with {test_pmcid}")
    print("=" * 80)
    
    sections = fetcher.fetch_full_text(test_pmcid)
    
    if sections:
        print("\n✓ Successfully fetched full text")
        print(f"\nSections found:")
        for section, content in sections.items():
            print(f"  - {section}: {len(content)} characters")
            if content:
                print(f"    Preview: {content[:200]}...")
        return True
    else:
        print("\n✗ Failed to fetch full text")
        return False


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test_pmc_fetcher()
