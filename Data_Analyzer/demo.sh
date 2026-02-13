#!/bin/bash
#
# Quick Demo - Data Analyzer
# ===========================
# Analyzes the example PubMed results
#

set -e

cd "$(dirname "$0")"

echo "=================================="
echo "DATA ANALYZER - Quick Demo"
echo "=================================="
echo ""

# Find most recent pubmed_results JSON
LATEST_JSON=$(ls -t ../pubmed_results_*.json 2>/dev/null | head -1)

if [ -z "$LATEST_JSON" ]; then
    echo "ERROR: No PubMed results found"
    echo "Run the fetcher first: bash test_fetcher_integrator.sh"
    exit 1
fi

echo "📄 Input: $LATEST_JSON"
echo ""

# Run analyzer
python3 paper_analyzer.py "$LATEST_JSON"

# Show output
LATEST_CSV=$(ls -t output/classified_papers_*.csv 2>/dev/null | head -1)

if [ -f "$LATEST_CSV" ]; then
    echo ""
    echo "📊 Results preview:"
    echo "-----------------------------------"
    head -3 "$LATEST_CSV" | cut -c1-120
    echo "..."
    echo ""
    echo "✓ Full results in: $LATEST_CSV"
fi
