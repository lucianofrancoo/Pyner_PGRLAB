#!/usr/bin/env python3
"""
Test script for Paper Visualization
=====================================
Generates an interactive HTML network graph from fetcher/analyzer output.

Usage:
    # Auto-detects classified TSV if it exists:
    python3 test_visualization.py output/pubmed_fetch_20260303_145749.json

    # Force classified TSV:
    python3 test_visualization.py output/pubmed_fetch_20260303_145749.json \\
        --classified output/classified_papers_20260303_145053.tsv

    # Custom output:
    python3 test_visualization.py output/pubmed_fetch_20260303_145749.json -o my_viz.html
"""

import sys
import os
import argparse
import logging
from pathlib import Path

ROOT_DIR = Path(__file__).parent
sys.path.insert(0, str(ROOT_DIR))

from Data_visualization.paper_visualizer import PaperVisualizer

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(description='Generate interactive HTML paper network')
    parser.add_argument('json_path', help='Path to pubmed_fetch_*.json')
    parser.add_argument('--classified', '-c', help='Path to classified_papers_*.tsv (auto-detected if omitted)', default=None)
    parser.add_argument('--output', '-o', help='Output HTML path', default=None)
    
    args = parser.parse_args()
    
    if not os.path.exists(args.json_path):
        logger.error(f"File not found: {args.json_path}")
        sys.exit(1)
    
    # Output path
    if args.output:
        output_path = args.output
    else:
        import re
        basename = os.path.basename(args.json_path)
        match = re.search(r'pubmed_fetch_(\d+_\d+)', basename)
        timestamp = match.group(1) if match else 'output'
        output_path = os.path.join(os.path.dirname(args.json_path) or 'output', f'paper_network_{timestamp}.html')
    
    # Generate
    print("\n" + "=" * 60)
    print("🔬 PYNER — Paper Network Visualization")
    print("=" * 60)
    
    viz = PaperVisualizer()
    n_papers = viz.load_data(args.json_path)
    print(f"  📄 Loaded {n_papers} papers from JSON")
    
    # Auto-detect or use explicit classified path
    classified_path = args.classified
    if not classified_path:
        classified_path = viz.auto_detect_classified(args.json_path)
    
    if classified_path and os.path.exists(classified_path):
        n_classified = viz.load_classified_data(classified_path)
        print(f"  📊 Auto-loaded {n_classified} classified results")
        print(f"  ✓ Mode: Enriched (hierarchical by relevance)")
    else:
        print(f"  ✓ Mode: Basic (no classified data found)")
    
    final_path = viz.generate_html(output_path)
    print(f"\n  ✓ Generated: {final_path}")
    print(f"  📦 Size: {os.path.getsize(final_path) / 1024:.1f} KB")
    print("=" * 60)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
