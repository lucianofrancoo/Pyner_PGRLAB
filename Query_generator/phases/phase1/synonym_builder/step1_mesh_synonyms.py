#!/usr/bin/env python3
"""
Step 1: Extract synonyms from MeSH ontology (desc2026.xml)

This script extracts official synonyms from the MeSH descriptor XML file.
These are curated by NLM and represent genuine synonyms, not just related terms.

Output: step1_mesh_synonyms.json
"""

import json
import logging
from pathlib import Path
from lxml import etree
from collections import defaultdict
from datetime import datetime

import config

# Setup logging
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def extract_mesh_synonyms(xml_path: Path) -> dict:
    """
    Extract synonyms from MeSH descriptor XML file.

    Structure of MeSH XML:
        <DescriptorRecord>
            <DescriptorUI>D018551</DescriptorUI>
            <DescriptorName>
                <String>Solanum lycopersicum</String>
            </DescriptorName>
            <ConceptList>
                <Concept PreferredConceptYN="Y">
                    <TermList>
                        <Term ConceptPreferredTermYN="Y">
                            <String>Solanum lycopersicum</String>  <!-- Preferred term -->
                        </Term>
                        <Term ConceptPreferredTermYN="N">
                            <String>Tomato</String>  <!-- SYNONYM -->
                        </Term>
                        <Term ConceptPreferredTermYN="N">
                            <String>Lycopersicon esculentum</String>  <!-- SYNONYM -->
                        </Term>
                    </TermList>
                </Concept>
            </ConceptList>
        </DescriptorRecord>

    Args:
        xml_path: Path to desc2026.xml

    Returns:
        dict: {
            "preferred_term_lowercase": {
                "mesh_ui": "D018551",
                "preferred_term": "Solanum lycopersicum",
                "synonyms": ["Tomato", "Tomatoes", "Lycopersicon esculentum"],
                "source": "mesh_ontology",
                "extracted_at": "2026-03-05T15:30:00"
            }
        }
    """
    logger.info(f"Starting MeSH synonym extraction from: {xml_path}")

    if not xml_path.exists():
        raise FileNotFoundError(f"MeSH file not found: {xml_path}")

    synonyms_dict = {}
    total_descriptors = 0
    descriptors_with_synonyms = 0
    total_synonym_count = 0

    # Parse XML efficiently with iterparse
    context = etree.iterparse(str(xml_path), events=('end',), tag='DescriptorRecord')

    for event, elem in context:
        total_descriptors += 1

        # Extract descriptor UI and name
        descriptor_ui_elem = elem.find('.//DescriptorUI')
        descriptor_name_elem = elem.find('.//DescriptorName/String')

        if descriptor_ui_elem is None or descriptor_name_elem is None:
            logger.warning(f"Skipping descriptor {total_descriptors}: missing UI or Name")
            elem.clear()
            continue

        descriptor_ui = descriptor_ui_elem.text
        preferred_term = descriptor_name_elem.text

        # Extract all terms from ConceptList
        all_terms = []
        preferred_terms = []

        for term in elem.xpath('.//ConceptList//Term'):
            term_string_elem = term.find('String')
            if term_string_elem is None:
                continue

            term_string = term_string_elem.text
            is_preferred = term.get('ConceptPreferredTermYN', 'N')

            if is_preferred == 'Y':
                preferred_terms.append(term_string)
            else:
                all_terms.append(term_string)

        # Extract synonyms (non-preferred terms)
        synonyms = []
        for term in all_terms:
            # Avoid duplicates and avoid adding the preferred term itself
            if term not in synonyms and term != preferred_term:
                synonyms.append(term)

        # Store if we have synonyms
        if synonyms:
            descriptors_with_synonyms += 1
            total_synonym_count += len(synonyms)

            # Use lowercase key for easier lookup
            key = preferred_term.lower()

            synonyms_dict[key] = {
                "mesh_ui": descriptor_ui,
                "preferred_term": preferred_term,
                "synonyms": synonyms,
                "synonym_count": len(synonyms),
                "source": "mesh_ontology",
                "extracted_at": datetime.now().isoformat()
            }

        # Progress logging
        if total_descriptors % config.LOG_PROGRESS_INTERVAL == 0:
            logger.info(f"Processed {total_descriptors:,} descriptors, "
                       f"found {descriptors_with_synonyms:,} with synonyms")

        # Clean up memory
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    logger.info("=" * 60)
    logger.info("MeSH Synonym Extraction Complete")
    logger.info("=" * 60)
    logger.info(f"Total descriptors processed:     {total_descriptors:,}")
    logger.info(f"Descriptors with synonyms:       {descriptors_with_synonyms:,}")
    logger.info(f"Total synonyms extracted:        {total_synonym_count:,}")
    logger.info(f"Average synonyms per descriptor: {total_synonym_count/descriptors_with_synonyms:.2f}")

    return synonyms_dict


def save_synonyms(synonyms_dict: dict, output_path: Path):
    """Save synonyms dictionary to JSON file."""
    logger.info(f"Saving synonyms to: {output_path}")

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(synonyms_dict, f, indent=2, ensure_ascii=False)

    # Calculate file size
    file_size_mb = output_path.stat().st_size / (1024 * 1024)
    logger.info(f"Saved {len(synonyms_dict):,} terms to {output_path.name} ({file_size_mb:.2f} MB)")


def validate_known_examples(synonyms_dict: dict):
    """Validate extraction with known examples."""
    logger.info("=" * 60)
    logger.info("Validating with known examples:")
    logger.info("=" * 60)

    test_cases = [
        "arabidopsis",
        "solanum lycopersicum",
        "rna-seq",
        "drought",
        "droughts",
        "oxidative stress"
    ]

    for term in test_cases:
        key = term.lower()
        if key in synonyms_dict:
            data = synonyms_dict[key]
            logger.info(f"\n✅ Found: {term}")
            logger.info(f"   Preferred term: {data['preferred_term']}")
            logger.info(f"   MeSH UI: {data['mesh_ui']}")
            logger.info(f"   Synonyms ({data['synonym_count']}):")
            for syn in data['synonyms'][:5]:  # Show first 5
                logger.info(f"      • {syn}")
            if data['synonym_count'] > 5:
                logger.info(f"      ... and {data['synonym_count'] - 5} more")
        else:
            logger.warning(f"❌ Not found: {term}")


def main():
    """Main execution."""
    logger.info("=" * 60)
    logger.info("STEP 1: MeSH Synonym Extraction")
    logger.info("=" * 60)
    logger.info(f"MeSH XML file: {config.MESH_DESC_FILE}")
    logger.info(f"Output file: {config.MESH_SYNONYMS_OUTPUT}")
    logger.info("")

    # Extract synonyms
    synonyms_dict = extract_mesh_synonyms(config.MESH_DESC_FILE)

    # Save to JSON
    save_synonyms(synonyms_dict, config.MESH_SYNONYMS_OUTPUT)

    # Validate with known examples
    validate_known_examples(synonyms_dict)

    logger.info("")
    logger.info("=" * 60)
    logger.info("✅ Step 1 Complete!")
    logger.info("=" * 60)
    logger.info(f"Next step: Implement step2_abstract_mining.py")


if __name__ == "__main__":
    main()
