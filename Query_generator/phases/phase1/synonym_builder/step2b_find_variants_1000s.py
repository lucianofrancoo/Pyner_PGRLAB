#!/usr/bin/env python3
"""
Step 2B modificado para usar step2_mesh_consolidated_parallel_1000s.json
"""
import sys
sys.path.insert(0, '.')

# Import everything from original step2b
from step2b_find_variants import *
import config

# Override the main function to use 1000-sample file
def main():
    """Main execution with 1000-sample data."""
    logger.info("=" * 60)
    logger.info("STEP 2B: TEXT VARIANT DETECTION (PARALLEL) - 1000 SAMPLES")
    logger.info("=" * 60)
    logger.info(f"Started: {datetime.now()}")
    
    input_file = config.OUTPUT_DIR / "step2_mesh_consolidated_parallel_1000s.json"
    output_file = config.OUTPUT_DIR / "step2_text_variants_1000s.json"
    
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output file: {output_file}")
    logger.info("")

    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        return

    logger.info("Loading consolidated data (1000 samples)...")
    with open(input_file, 'r', encoding='utf-8') as f:
        consolidated_data = json.load(f)

    file_size_mb = input_file.stat().st_size / (1024 * 1024)
    logger.info(f"Loaded {len(consolidated_data):,} MeSH terms ({file_size_mb:.2f} MB)")
    logger.info("")

    # Process all terms
    results = process_all_terms(
        consolidated_data,
        similarity_threshold=0.70,
        min_confidence="low",
        n_workers=None  # Auto-detect
    )

    # Save results
    save_results(results, output_file)
    show_examples(results, n=10)

    logger.info("")
    logger.info("=" * 60)
    logger.info("✅ Step 2B Complete!")
    logger.info("=" * 60)
    logger.info(f"Output: {output_file}")

if __name__ == "__main__":
    main()
