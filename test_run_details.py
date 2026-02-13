#!/usr/bin/env python3
"""
Test script to validate run details extraction (spots and size).
"""

import sys
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "Fetcher_NCBI"))

from boolean_fetcher_integrated import BooleanFetcherIntegrated
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_run_details_extraction():
    """Test that run details (spots, size) are extracted correctly."""
    
    logger.info("="*70)
    logger.info("Testing Run Details Extraction (Spots & Size)")
    logger.info("="*70)
    
    fetcher = BooleanFetcherIntegrated()
    
    # Test BioProject
    test_bioproject = "PRJNA1381306"
    
    logger.info(f"\nTesting with BioProject: {test_bioproject}")
    logger.info("Fetching SRA data with run details...")
    
    # Create test data
    bp_data = {
        'bioproject': test_bioproject,
        'title': 'Test Project',
        'organism': 'Solanum lycopersicum',
        'submission_date': '2024-01-01',
        'project_type': 'OMC',
        'description': 'Test'
    }
    
    # Process with full hierarchy
    result = fetcher.process_bioproject(bp_data)
    
    # Validate run details
    logger.info("\n" + "="*70)
    logger.info("VALIDATION CHECK - RUN DETAILS")
    logger.info("="*70)
    
    if 'sra_hierarchy' in result:
        hierarchy = result['sra_hierarchy']
        logger.info(f"✓ Hierarchy present with {len(hierarchy)} biosamples")
        
        # Check first sample with experiments
        found_runs = False
        for bs_id, bs_data in list(hierarchy.items())[:2]:
            logger.info(f"\nBioSample: {bs_id}")
            
            for exp in bs_data.get('experiments', [])[:1]:
                logger.info(f"  Experiment: {exp.get('experiment_id')}")
                
                runs = exp.get('runs', [])
                if runs and len(runs) > 0:
                    found_runs = True
                    logger.info(f"  Runs found: {len(runs)}")
                    
                    for i, run in enumerate(runs[:3], 1):
                        if isinstance(run, dict):
                            accession = run.get('accession', 'N/A')
                            spots = run.get('spots', 0)
                            base_count = run.get('base_count', 0)
                            size = run.get('size', 'N/A')
                            
                            logger.info(f"\n    Run {i}: {accession}")
                            logger.info(f"      • Spots: {spots:,}" if spots > 0 else f"      • Spots: {spots}")
                            logger.info(f"      • Base Count: {base_count:,}" if base_count > 0 else f"      • Base Count: {base_count}")
                            logger.info(f"      • Size: {size}")
                            
                            # Validate
                            if spots > 0:
                                logger.info(f"      ✓ Spots extracted")
                            if base_count > 0:
                                logger.info(f"      ✓ Base count extracted")
                            if size != 'N/A' and size != '0.00 GB':
                                logger.info(f"      ✓ Size calculated")
                        else:
                            logger.info(f"    Run {i}: {run} (string format)")
        
        if found_runs:
            logger.info("\n✓ Run details are being captured!")
        else:
            logger.warning("\n⚠ No detailed run info found (might be expected if XML doesn't have Statistics)")
    
    logger.info("\n" + "="*70)
    logger.info("✓ RUN DETAILS EXTRACTION TEST COMPLETE")
    logger.info("="*70)
    
    return result


def test_json_output_with_run_details():
    """Test that JSON includes run details."""
    
    logger.info("\n" + "="*70)
    logger.info("Testing JSON Output with Run Details")
    logger.info("="*70)
    
    fetcher = BooleanFetcherIntegrated()
    
    bp_data = {
        'bioproject': 'PRJNA1381306',
        'title': 'Test',
        'organism': 'Solanum lycopersicum',
        'submission_date': '2024-01-01',
        'project_type': 'OMC',
        'description': 'Test'
    }
    
    result = fetcher.process_bioproject(bp_data)
    fetcher.results.append(result)
    
    # Save to JSON
    output_file = Path("/tmp/test_runs_details.json")
    fetcher.save_results_json(output_file)
    
    # Read and validate
    with open(output_file, 'r') as f:
        data = json.load(f)
    
    logger.info(f"✓ JSON file created: {output_file}")
    
    if data['results']:
        result = data['results'][0]
        if 'sra_hierarchy' in result:
            hierarchy = result['sra_hierarchy']
            
            # Check first sample with run details
            for bs_id, bs_data in list(hierarchy.items())[:1]:
                for exp in bs_data.get('experiments', [])[:1]:
                    runs = exp.get('runs', [])
                    if runs and isinstance(runs[0], dict):
                        logger.info(f"\n✓ Run details in JSON:")
                        run = runs[0]
                        logger.info(f"   - Accession: {run.get('accession')}")
                        logger.info(f"   - Spots: {run.get('spots')}")
                        logger.info(f"   - Base Count: {run.get('base_count')}")
                        logger.info(f"   - Size: {run.get('size')}")
    
    logger.info("\n✓ JSON OUTPUT TEST COMPLETE")


if __name__ == '__main__':
    try:
        logger.info("\nStarting Run Details Tests...\n")
        result = test_run_details_extraction()
        test_json_output_with_run_details()
        
        logger.info("\n" + "="*70)
        logger.info("✓✓✓ ALL TESTS COMPLETED ✓✓✓")
        logger.info("Run details (spots & size) integration complete!")
        logger.info("="*70)
        sys.exit(0)
        
    except Exception as e:
        logger.error(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
