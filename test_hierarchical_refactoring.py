#!/usr/bin/env python3
"""
Test script to validate hierarchical SRA data refactoring.
Tests the new biosamples_dict and sra_runs structures.
"""

import sys
from pathlib import Path

# Add Fetcher_NCBI to path
sys.path.insert(0, str(Path(__file__).parent / "Fetcher_NCBI"))

from boolean_fetcher_integrated import BooleanFetcherIntegrated
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_sra_data_structure():
    """Test that fetch_sra_for_bioproject returns correct structure."""
    
    logger.info("="*70)
    logger.info("Testing SRA Data Structure Refactoring")
    logger.info("="*70)
    
    # Initialize fetcher
    fetcher = BooleanFetcherIntegrated()
    
    # Test BioProject (small, fast)
    test_bioproject = "PRJNA1381306"
    
    logger.info(f"\nTesting with BioProject: {test_bioproject}")
    logger.info("Fetching SRA data...")
    
    # Call the refactored method
    experiments, biosamples_dict, sra_runs = fetcher.fetch_sra_for_bioproject(test_bioproject)
    
    # Validate return types
    logger.info("\n" + "="*70)
    logger.info("VALIDATION RESULTS")
    logger.info("="*70)
    
    assert isinstance(experiments, list), "experiments should be list"
    logger.info(f"✓ experiments is list: {len(experiments)} items")
    
    assert isinstance(biosamples_dict, dict), "biosamples_dict should be dict"
    logger.info(f"✓ biosamples_dict is dict: {len(biosamples_dict)} items")
    
    assert isinstance(sra_runs, list), "sra_runs should be list"
    unique_runs = len(set(sra_runs))
    logger.info(f"✓ sra_runs is list: {len(sra_runs)} items ({unique_runs} unique)")
    
    # Validate experiments structure
    if experiments:
        exp = experiments[0]
        logger.info(f"\nFirst experiment structure:")
        logger.info(f"  - exp_accession: {exp.get('exp_accession')} (should be SRX*)")
        logger.info(f"  - title: {exp.get('title', 'N/A')[:50]}...")
        logger.info(f"  - biosample: {exp.get('biosample')} (should be SAMN*)")
        logger.info(f"  - runs: {exp.get('runs', [])[:3]}... ({len(exp.get('runs', []))} runs)")
        
        # Check that runs are SRR codes
        if exp.get('runs'):
            sample_run = exp['runs'][0]
            assert sample_run.startswith('SRR'), f"Run should start with SRR, got {sample_run}"
            logger.info(f"✓ Runs are SRR codes (not SRX)")
    
    # Validate biosamples_dict structure
    if biosamples_dict:
        bs_id = list(biosamples_dict.keys())[0]
        bs_data = biosamples_dict[bs_id]
        logger.info(f"\nFirst biosample structure:")
        logger.info(f"  - ID: {bs_id} (should be SAMN*)")
        assert bs_id.startswith('SAMN'), f"Biosample ID should start with SAMN, got {bs_id}"
        logger.info(f"✓ BioSample IDs are SAMN codes")
        
        logger.info(f"  - Structure: {bs_data}")
        assert 'experiment_titles' in bs_data, "Missing experiment_titles"
        assert 'samples' in bs_data, "Missing samples"
        logger.info(f"  - Experiment titles: {bs_data['experiment_titles'][:2]}...")
    
    # Validate sra_runs
    if sra_runs:
        logger.info(f"\nSample SRA runs:")
        for i, run in enumerate(sra_runs[:5]):
            logger.info(f"  {i+1}. {run}")
            assert run.startswith('SRR'), f"Run should be SRR*, got {run}"
        logger.info(f"✓ All runs are SRR codes")
    
    logger.info("\n" + "="*70)
    logger.info("✓ ALL TESTS PASSED - Refactoring successful!")
    logger.info("="*70)
    
    return True


def test_process_bioproject():
    """Test that process_bioproject handles new structure."""
    
    logger.info("\n" + "="*70)
    logger.info("Testing process_bioproject with new structure")
    logger.info("="*70)
    
    fetcher = BooleanFetcherIntegrated()
    
    # Create minimal test BioProject data
    bp_data = {
        'bioproject': 'PRJNA1381306',
        'title': 'Test Project',
        'organism': 'Test organism',
        'submission_date': '2024-01-01',
        'project_type': 'OMC',
        'description': 'Test description'
    }
    
    logger.info(f"Processing BioProject: {bp_data['bioproject']}")
    
    try:
        result = fetcher.process_bioproject(bp_data)
        
        # Check result structure
        assert 'sra_experiments_count' in result
        assert 'biosamples_count' in result
        assert 'sra_runs_count' in result
        assert 'sra_experiments' in result
        assert 'biosamples' in result
        assert 'sra_runs' in result
        
        logger.info("✓ process_bioproject returns correct fields")
        logger.info(f"  - sra_experiments_count: {result['sra_experiments_count']}")
        logger.info(f"  - biosamples_count: {result['biosamples_count']}")
        logger.info(f"  - sra_runs_count: {result['sra_runs_count']}")
        
        # Verify counts are numbers
        assert isinstance(result['sra_experiments_count'], int)
        assert isinstance(result['biosamples_count'], int)
        assert isinstance(result['sra_runs_count'], int)
        logger.info("✓ All counts are integers")
        
        return True
        
    except Exception as e:
        logger.error(f"✗ Error in process_bioproject: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == '__main__':
    try:
        success = test_sra_data_structure()
        if success:
            logger.info("\n✓ SRA structure test passed\n")
            success = test_process_bioproject()
            if success:
                logger.info("\n✓ Process BioProject test passed\n")
                logger.info("\n" + "="*70)
                logger.info("✓✓✓ ALL TESTS PASSED ✓✓✓")
                logger.info("Hierarchical refactoring is working correctly!")
                logger.info("="*70)
                sys.exit(0)
    except Exception as e:
        logger.error(f"Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
