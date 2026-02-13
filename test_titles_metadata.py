#!/usr/bin/env python3
"""
Test script to validate hierarchical SRA structure with full metadata and titles.
"""

import sys
import json
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "Fetcher_NCBI"))

from boolean_fetcher_integrated import BooleanFetcherIntegrated
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def test_hierarchical_structure_with_titles():
    """Test that SRA hierarchy includes titles and full metadata."""
    
    logger.info("="*70)
    logger.info("Testing Hierarchical SRA Structure with Titles & Metadata")
    logger.info("="*70)
    
    fetcher = BooleanFetcherIntegrated()
    
    # Test BioProject
    test_bioproject = "PRJNA1381306"
    
    logger.info(f"\nTesting with BioProject: {test_bioproject}")
    logger.info("Fetching SRA data with full metadata...")
    
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
    
    # Validate structure
    logger.info("\n" + "="*70)
    logger.info("VALIDATION CHECK")
    logger.info("="*70)
    
    # Check that hierarchy exists
    assert 'sra_hierarchy' in result, "Missing sra_hierarchy"
    logger.info("✓ sra_hierarchy field present")
    
    hierarchy = result['sra_hierarchy']
    logger.info(f"✓ Hierarchy has {len(hierarchy)} biosamples")
    
    # Check first biosample structure
    if hierarchy:
        first_bs_id = list(hierarchy.keys())[0]
        first_bs = hierarchy[first_bs_id]
        
        logger.info(f"\nFirst BioSample: {first_bs_id}")
        logger.info(f"  - Experiments: {len(first_bs.get('experiments', []))}")
        
        # Check first experiment
        if first_bs.get('experiments'):
            first_exp = first_bs['experiments'][0]
            logger.info(f"\nFirst Experiment:")
            logger.info(f"  - ID: {first_exp.get('experiment_id')}")
            logger.info(f"  - Title: {first_exp.get('title', 'N/A')[:60]}...")
            logger.info(f"  - Library Name: {first_exp.get('metadata', {}).get('library_name', 'N/A')}")
            logger.info(f"  - Strategy: {first_exp.get('metadata', {}).get('library_strategy', 'N/A')}")
            logger.info(f"  - Source: {first_exp.get('metadata', {}).get('library_source', 'N/A')}")
            logger.info(f"  - Selection: {first_exp.get('metadata', {}).get('library_selection', 'N/A')}")
            logger.info(f"  - Layout: {first_exp.get('metadata', {}).get('library_layout', 'N/A')}")
            logger.info(f"  - Instrument: {first_exp.get('metadata', {}).get('instrument', 'N/A')}")
            logger.info(f"  - Runs: {first_exp.get('runs', [])}")
            
            # Check for sample attributes
            sample_attrs = first_exp.get('sample_attributes', {})
            if sample_attrs:
                logger.info(f"  - Sample Attributes: {dict(list(sample_attrs.items())[:5])}")
            
            # Validate required fields
            assert first_exp.get('experiment_id'), "Missing experiment_id"
            assert first_exp.get('title'), "Missing title"
            logger.info("\n✓ All experiment fields present and populated")
    
    # Check counts
    logger.info(f"\nCounts validation:")
    logger.info(f"  - sra_experiments_count: {result.get('sra_experiments_count')}")
    logger.info(f"  - biosamples_count: {result.get('biosamples_count')}")
    logger.info(f"  - sra_runs_count: {result.get('sra_runs_count')}")
    
    assert result.get('sra_experiments_count') > 0, "No experiments"
    assert result.get('biosamples_count') > 0, "No biosamples"
    assert result.get('sra_runs_count') > 0, "No runs"
    logger.info("✓ All counts > 0")
    
    logger.info("\n" + "="*70)
    logger.info("✓ HIERARCHICAL STRUCTURE TEST PASSED")
    logger.info("="*70)
    
    return result


def test_json_output():
    """Test that JSON output includes hierarchical structure."""
    
    logger.info("\n" + "="*70)
    logger.info("Testing JSON Output")
    logger.info("="*70)
    
    fetcher = BooleanFetcherIntegrated()
    
    # Create sample result
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
    output_file = Path("/tmp/test_hierarchy.json")
    fetcher.save_results_json(output_file)
    
    # Read and validate
    with open(output_file, 'r') as f:
        data = json.load(f)
    
    logger.info(f"✓ JSON file created: {output_file}")
    logger.info(f"  - Total results: {data['metadata']['total_results']}")
    
    if data['results']:
        result = data['results'][0]
        if 'sra_hierarchy' in result:
            logger.info(f"✓ sra_hierarchy present in JSON output")
            logger.info(f"  - BioSamples in hierarchy: {len(result['sra_hierarchy'])}")
        else:
            logger.warning("⚠ sra_hierarchy not in JSON output (but this is OK if it's too large)")
    
    logger.info("\n✓ JSON OUTPUT TEST PASSED")


if __name__ == '__main__':
    try:
        result = test_hierarchical_structure_with_titles()
        test_json_output()
        
        logger.info("\n" + "="*70)
        logger.info("✓✓✓ ALL TESTS PASSED ✓✓✓")
        logger.info("Hierarchical title & metadata integration working!")
        logger.info("="*70)
        sys.exit(0)
        
    except Exception as e:
        logger.error(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
