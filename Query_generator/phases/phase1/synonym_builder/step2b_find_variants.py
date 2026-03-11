#!/usr/bin/env python3
"""
Step 2B: Find Text Variants from Collected Samples

Analyzes the sample texts collected in Step 2A to find orthographic variants
of each MeSH term using fuzzy matching and n-gram analysis.

Input: output/step2_mesh_consolidated.json (~1.2 GB)
Output: output/step2_text_variants.json (~50-100 MB)

This script identifies variants like:
- "RNA-Seq" → "RNAseq", "RNA seq", "rna sequencing"
- "Arabidopsis" → "A. thaliana", "arabidopsis thaliana"
- "Drought" → "water deficit", "water scarcity"

Created: 2026-03-06
"""

import json
import re
import logging
from pathlib import Path
from collections import Counter
from datetime import datetime
from difflib import SequenceMatcher
from multiprocessing import Pool, cpu_count
from functools import partial
import os

try:
    from Levenshtein import distance as levenshtein_distance
    HAS_LEVENSHTEIN = True
except ImportError:
    HAS_LEVENSHTEIN = False
    logging.warning("python-Levenshtein not installed. Using SequenceMatcher instead (slower).")

import config

# Setup logging
logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL),
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def calculate_similarity(str1: str, str2: str) -> float:
    """
    Calculate similarity between two strings (0.0 to 1.0).

    Uses Levenshtein distance if available, otherwise SequenceMatcher.
    """
    if not str1 or not str2:
        return 0.0

    if HAS_LEVENSHTEIN:
        # Levenshtein distance normalized to 0-1
        max_len = max(len(str1), len(str2))
        if max_len == 0:
            return 1.0
        return 1.0 - (levenshtein_distance(str1, str2) / max_len)
    else:
        # SequenceMatcher (slower but no dependencies)
        return SequenceMatcher(None, str1, str2).ratio()


def normalize_string(text: str) -> str:
    """
    Normalize string for comparison.

    Removes: hyphens, underscores, spaces
    Converts to lowercase
    """
    return text.lower().replace('-', '').replace('_', '').replace(' ', '')


def is_likely_variant(token: str, mesh_term: str, threshold: float = 0.75) -> bool:
    """
    Determine if a token is likely a variant of the MeSH term.

    Uses normalized string comparison with fuzzy matching.

    Args:
        token: Candidate token from text
        mesh_term: MeSH term to compare against
        threshold: Minimum similarity score (0.0 to 1.0, default: 0.75)

    Returns:
        True if likely a variant

    Note: Threshold of 0.75 focuses on high-precision matches:
        - High threshold reduces false positives (contextual words)
        - Genuine variants like "RNA-Seq" vs "RNAseq" still match (>0.85)
        - Multi-word variants are handled separately (e.g., "A. thaliana")
    """
    # Normalize both strings
    token_norm = normalize_string(token)
    mesh_norm = normalize_string(mesh_term)

    # Skip if too short
    if len(token_norm) < 3 or len(mesh_norm) < 3:
        return False

    # Skip if token and mesh_term are identical (already counted)
    if token_norm == mesh_norm:
        return True

    # For multi-word MeSH terms, check if token contains or is contained
    # This helps with variants like "Arabidopsis thaliana" vs "A. thaliana"
    mesh_words = mesh_norm.split()
    token_words = token_norm.split()

    # If MeSH is multi-word, require token to share at least one significant word
    if len(mesh_words) > 1 and len(token_words) > 0:
        # Check if token shares words with mesh term
        shared_words = set(mesh_words) & set(token_words)
        if not shared_words:
            # No shared words, unlikely to be a variant
            return False

    # Calculate similarity
    similarity = calculate_similarity(token_norm, mesh_norm)

    return similarity >= threshold


def tokenize_text(text: str, max_ngram: int = 4) -> list:
    """
    Extract tokens from text using n-grams, filtering common stop words.

    Extracts single words and multi-word phrases (up to max_ngram words).

    Args:
        text: Text to tokenize
        max_ngram: Maximum n-gram size (default: 4)

    Returns:
        List of tokens
    """
    # Common stop words to filter out from n-grams
    # These are words that often appear around terms but aren't part of synonyms
    STOP_WORDS = {
        'the', 'a', 'an', 'and', 'or', 'but', 'in', 'on', 'at', 'to', 'for',
        'of', 'with', 'by', 'from', 'as', 'is', 'was', 'are', 'were', 'been',
        'be', 'have', 'has', 'had', 'do', 'does', 'did', 'will', 'would',
        'should', 'could', 'may', 'might', 'can', 'this', 'that', 'these',
        'those', 'such', 'which', 'who', 'whom', 'what', 'where', 'when',
        'why', 'how', 'all', 'each', 'every', 'both', 'few', 'more', 'most',
        'other', 'some', 'such', 'than', 'too', 'very', 'we', 'us', 'our',
        'using', 'used', 'between', 'during', 'after', 'before', 'under',
        'over', 'through', 'into', 'upon'
    }

    # Split by separator
    parts = text.split('||')

    tokens = []

    for part in parts:
        # Clean text: keep letters, numbers, spaces, hyphens
        part = re.sub(r'[^\w\s\-]', ' ', part)

        # Split into words
        words = part.lower().split()

        # Extract n-grams
        for n in range(1, min(max_ngram + 1, len(words) + 1)):
            for i in range(len(words) - n + 1):
                ngram_words = words[i:i+n]

                # Skip n-grams that start or end with stop words
                if ngram_words[0] in STOP_WORDS or ngram_words[-1] in STOP_WORDS:
                    continue

                ngram = ' '.join(ngram_words)

                # Skip very short tokens
                if len(ngram) >= 3:
                    tokens.append(ngram)

    return tokens


def extract_organism_variants(text: str, genus: str) -> list:
    """
    Extract organism name variants using specific patterns.

    Patterns:
    - Full name: "Arabidopsis thaliana"
    - Abbreviated: "A. thaliana"
    - Without period: "A thaliana"

    Args:
        text: Text to search
        genus: Genus name (e.g., "Arabidopsis")

    Returns:
        List of variant strings found
    """
    # Common words that appear after genus but aren't species names
    NON_SPECIES_WORDS = {
        'and', 'or', 'the', 'a', 'an', 'in', 'on', 'at', 'to', 'for', 'of', 'with',
        'by', 'from', 'as', 'is', 'was', 'are', 'were', 'been', 'be', 'have',
        'has', 'had', 'genome', 'genomes', 'gene', 'genes', 'protein', 'proteins',
        'plants', 'plant', 'seed', 'seeds', 'leaf', 'leaves', 'root', 'roots',
        'flower', 'flowers', 'floral', 'seedling', 'seedlings', 'mutant', 'mutants',
        'cdna', 'mrna', 'dna', 'rna', 'containing', 'encoding', 'derived'
    }

    variants = []

    # Escape special regex characters in genus name
    genus_escaped = re.escape(genus)

    # Pattern 1: Full genus + species (lowercase species name, likely Latin)
    # Species names are typically lowercase in scientific text
    pattern1 = rf'\b{genus_escaped}\s+([a-z]\w+)\b'
    try:
        matches = re.findall(pattern1, text, re.IGNORECASE)
        for match in matches:
            # Filter out common non-species words
            if isinstance(match, str):
                species_word = match.lower()
            else:
                species_word = match[0].lower() if match else ""

            if species_word and species_word not in NON_SPECIES_WORDS:
                # Reconstruct full binomial name
                variants.append(f"{genus} {species_word}")
    except re.error:
        pass  # Skip if pattern is invalid

    # Pattern 2: Abbreviated genus + species (with period)
    # e.g., "A. thaliana"
    if len(genus) > 0:
        first_letter = re.escape(genus[0])
        pattern2 = rf'\b{first_letter}\.\s*([a-z]\w+)\b'
        try:
            matches = re.findall(pattern2, text, re.IGNORECASE)
            for match in matches:
                species_word = match.lower() if isinstance(match, str) else (match[0].lower() if match else "")
                if species_word and species_word not in NON_SPECIES_WORDS:
                    variants.append(f"{first_letter}. {species_word}")
        except re.error:
            pass

    return [v.strip() for v in variants]


def find_variants(mesh_term: str, sample_texts: list,
                  mesh_synonyms: list = None,
                  similarity_threshold: float = 0.75,
                  min_frequency: int = 2,
                  min_percentage: float = 0.02) -> dict:
    """
    Find orthographic variants of a MeSH term and its synonyms in sample texts.

    Args:
        mesh_term: The MeSH term to find variants for
        sample_texts: List of sample texts (up to 100)
        mesh_synonyms: List of official MeSH synonyms for this term
        similarity_threshold: Minimum similarity for fuzzy matching
        min_frequency: Minimum number of occurrences
        min_percentage: Minimum percentage of samples (0.0-1.0)

    Returns:
        Dictionary of variants with counts and percentages
    """
    if mesh_synonyms is None:
        mesh_synonyms = []

    all_tokens = []
    organism_tokens = []

    # Collect all terms to search for (main term + synonyms)
    terms_to_search = [mesh_term] + mesh_synonyms

    # Special handling for organism names (check if likely an organism)
    is_organism = False
    organism_genera = set()  # Track all genera from term and synonyms

    for term in terms_to_search:
        words = term.split()
        if len(words) >= 1 and words[0] and words[0][0].isupper():
            # Might be an organism name
            is_organism = True
            organism_genera.add(words[0])

    # Tokenize all sample texts
    for text in sample_texts:
        # Extract organism variants if applicable
        if is_organism:
            for genus in organism_genera:
                org_variants = extract_organism_variants(text, genus)
                organism_tokens.extend(org_variants)

        # Extract general tokens
        tokens = tokenize_text(text, max_ngram=config.MAX_NGRAM_SIZE)
        all_tokens.extend(tokens)

    # Count token frequencies
    token_counts = Counter(all_tokens)
    organism_counts = Counter(organism_tokens)

    # Find candidate variants
    variants = {}
    total_samples = len(sample_texts)

    # Process organism variants with more lenient matching
    if is_organism and organism_counts:
        for token, count in organism_counts.items():
            # For organism names, check if it's a valid binomial name
            token_lower = token.lower()
            genus_lower = genus.lower()

            # Accept if:
            # 1. Contains the full genus name (e.g., "Arabidopsis thaliana")
            # 2. OR starts with abbreviated genus (e.g., "A. thaliana")
            # AND looks like a binomial name (at least 2 words)
            is_valid_binomial = False
            if len(token.split()) >= 2:
                # Check full genus
                if genus_lower in token_lower:
                    is_valid_binomial = True
                # Check abbreviated (first letter + period)
                elif token_lower.startswith(genus_lower[0] + '.'):
                    is_valid_binomial = True

            if is_valid_binomial:
                percentage = count / total_samples
                if count >= min_frequency and percentage >= min_percentage:
                    variants[token.lower()] = {
                        "count": count,
                        "percentage": round(percentage, 4)
                    }

    # Process general tokens
    for token, count in token_counts.items():
        # Skip if already added as organism variant
        if token.lower() in variants:
            continue

        # Check if it's similar to the MeSH term OR any of its synonyms
        is_variant = False
        for search_term in terms_to_search:
            if is_likely_variant(token, search_term, threshold=similarity_threshold):
                is_variant = True
                break

        if is_variant:
            percentage = count / total_samples

            # Apply filters
            if count >= min_frequency and percentage >= min_percentage:
                variants[token] = {
                    "count": count,
                    "percentage": round(percentage, 4)
                }

    return variants


def calculate_confidence(percentage: float, count: int) -> str:
    """
    Calculate confidence level for a variant.

    Levels:
    - very_high: >50% of samples, >50 occurrences
    - high: >20% of samples, >20 occurrences
    - medium: >10% of samples, >10 occurrences
    - low: >5% of samples, >5 occurrences
    - very_low: anything else

    Args:
        percentage: Percentage of samples (0.0-1.0)
        count: Number of occurrences

    Returns:
        Confidence level string
    """
    if percentage >= 0.5 and count >= 50:
        return "very_high"
    elif percentage >= 0.2 and count >= 20:
        return "high"
    elif percentage >= 0.1 and count >= 10:
        return "medium"
    elif percentage >= 0.05 and count >= 5:
        return "low"
    else:
        return "very_low"


def process_single_term(args):
    """
    Process a single MeSH term to find variants.

    This function is designed to be called by multiprocessing workers.

    Args:
        args: Tuple of (key, term_data, similarity_threshold, min_confidence_idx)

    Returns:
        Tuple of (key, result_dict) or None if no variants found
    """
    key, term_data, similarity_threshold, min_confidence_idx = args

    confidence_order = ["very_low", "low", "medium", "high", "very_high"]

    mesh_term = term_data['preferred_term']
    sample_texts = term_data.get('sample_texts', [])
    article_count = term_data.get('article_count', 0)
    mesh_synonyms = term_data.get('mesh_synonyms', [])

    # Skip if no sample texts
    if not sample_texts:
        return None

    # Find variants (including variants of MeSH synonyms)
    variants = find_variants(
        mesh_term,
        sample_texts,
        mesh_synonyms=mesh_synonyms,
        similarity_threshold=similarity_threshold,
        min_frequency=config.ABSTRACT_MIN_FREQUENCY,
        min_percentage=config.ABSTRACT_MIN_PERCENTAGE
    )

    # Calculate confidence and filter
    filtered_variants = {}
    for variant, data in variants.items():
        confidence = calculate_confidence(data['percentage'], data['count'])
        confidence_idx = confidence_order.index(confidence)

        # Only include if meets minimum confidence
        if confidence_idx >= min_confidence_idx:
            filtered_variants[variant] = {
                **data,
                "confidence": confidence
            }

    # Return results if we found variants
    if filtered_variants:
        return (key, {
            "mesh_term": mesh_term,
            "article_count": article_count,
            "sample_count": len(sample_texts),
            "text_variants": filtered_variants
        })

    return None


def process_all_terms(consolidated_data: dict,
                      similarity_threshold: float = 0.75,
                      min_confidence: str = "low",
                      n_workers: int = None) -> dict:
    """
    Process all MeSH terms to find variants using parallel processing.

    Args:
        consolidated_data: Data from Step 2A
        similarity_threshold: Minimum similarity for fuzzy matching
        min_confidence: Minimum confidence level to include
        n_workers: Number of parallel workers (default: CPU count - 8)

    Returns:
        Dictionary of variants for each term
    """
    logger.info("=" * 60)
    logger.info("Starting PARALLEL variant detection for all MeSH terms")
    logger.info("=" * 60)
    logger.info(f"Total terms to process: {len(consolidated_data):,}")
    logger.info(f"Similarity threshold: {similarity_threshold}")
    logger.info(f"Minimum confidence: {min_confidence}")

    # Determine number of workers
    if n_workers is None:
        n_workers = max(1, cpu_count() - 8)  # Leave 8 cores free

    logger.info(f"CPU cores available: {cpu_count()}")
    logger.info(f"Parallel workers: {n_workers}")
    logger.info("")

    confidence_order = ["very_low", "low", "medium", "high", "very_high"]
    min_confidence_idx = confidence_order.index(min_confidence)

    # Prepare arguments for parallel processing
    args_list = [
        (key, term_data, similarity_threshold, min_confidence_idx)
        for key, term_data in consolidated_data.items()
    ]

    total_terms = len(args_list)
    results = {}
    processed = 0
    with_variants = 0
    total_variants = 0

    # Process in parallel with progress tracking
    logger.info(f"Processing {total_terms:,} terms in parallel...")

    with Pool(processes=n_workers) as pool:
        # Use imap_unordered for better progress tracking
        for result in pool.imap_unordered(process_single_term, args_list, chunksize=10):
            processed += 1

            if result is not None:
                key, term_result = result
                with_variants += 1
                total_variants += len(term_result['text_variants'])
                results[key] = term_result

            # Progress logging
            if processed % 1000 == 0:
                logger.info(f"Processed {processed:,}/{total_terms:,} terms "
                           f"({processed/total_terms*100:.1f}%) - "
                           f"Found {with_variants:,} terms with variants")

    logger.info("")
    logger.info("=" * 60)
    logger.info("Variant Detection Complete")
    logger.info("=" * 60)
    logger.info(f"Terms processed: {processed:,}")
    logger.info(f"Terms with variants: {with_variants:,} ({with_variants/processed*100:.1f}%)")
    logger.info(f"Total variants found: {total_variants:,}")
    if with_variants > 0:
        logger.info(f"Average variants/term: {total_variants/with_variants:.1f}")
    logger.info("")

    return results


def save_results(results: dict, output_path: Path):
    """Save variant results to JSON file."""
    logger.info(f"Saving results to: {output_path}")

    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    file_size_mb = output_path.stat().st_size / (1024 * 1024)
    logger.info(f"Saved {len(results):,} terms with variants to {output_path.name} ({file_size_mb:.2f} MB)")


def show_examples(results: dict, n: int = 10):
    """Show example variants for top N terms."""
    logger.info("=" * 60)
    logger.info(f"Example variants for top {n} terms:")
    logger.info("=" * 60)

    # Sort by number of variants
    sorted_results = sorted(
        results.items(),
        key=lambda x: len(x[1]['text_variants']),
        reverse=True
    )

    for i, (key, data) in enumerate(sorted_results[:n], 1):
        logger.info(f"\n{i}. {data['mesh_term']}")
        logger.info(f"   Articles: {data['article_count']:,}")
        logger.info(f"   Variants found: {len(data['text_variants'])}")

        # Show top 5 variants by count
        sorted_variants = sorted(
            data['text_variants'].items(),
            key=lambda x: x[1]['count'],
            reverse=True
        )

        for variant, vdata in sorted_variants[:5]:
            logger.info(f"      • {variant:40s} - {vdata['count']:3d} occurrences "
                       f"({vdata['percentage']*100:5.1f}%) [{vdata['confidence']}]")

        if len(sorted_variants) > 5:
            logger.info(f"      ... and {len(sorted_variants) - 5} more")


def main():
    """Main execution."""
    logger.info("=" * 60)
    logger.info("STEP 2B: TEXT VARIANT DETECTION (PARALLEL)")
    logger.info("=" * 60)
    logger.info(f"Started: {datetime.now()}")
    logger.info(f"Input file: {config.OUTPUT_DIR / 'step2_mesh_consolidated.json'}")
    logger.info(f"Output file: {config.OUTPUT_DIR / 'step2_text_variants.json'}")
    logger.info("")

    # Load consolidated data from Step 2A
    input_file = config.OUTPUT_DIR / "step2_mesh_consolidated.json"

    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        logger.error("Please run Step 2A first (batch_download_consolidated.py)")
        return

    logger.info("Loading consolidated data from Step 2A...")
    with open(input_file, 'r', encoding='utf-8') as f:
        consolidated_data = json.load(f)

    file_size_mb = input_file.stat().st_size / (1024 * 1024)
    logger.info(f"Loaded {len(consolidated_data):,} MeSH terms ({file_size_mb:.2f} MB)")
    logger.info("")

    # Process all terms to find variants (in parallel)
    results = process_all_terms(
        consolidated_data,
        similarity_threshold=0.70,  # Balanced threshold (tested with sample data)
        min_confidence="low",  # Include low confidence and above
        n_workers=None  # Auto-detect (CPU count - 8)
    )

    # Save results
    output_file = config.OUTPUT_DIR / "step2_text_variants.json"
    save_results(results, output_file)

    # Show examples
    show_examples(results, n=10)

    logger.info("")
    logger.info("=" * 60)
    logger.info("✅ Step 2B Complete!")
    logger.info("=" * 60)
    logger.info(f"Output: {output_file}")
    logger.info(f"Next step: Review results and proceed to Step 3 (BioWordVec)")


if __name__ == "__main__":
    main()
