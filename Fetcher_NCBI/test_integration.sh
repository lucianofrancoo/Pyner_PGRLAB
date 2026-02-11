#!/bin/bash
#
# Integration Test: Query Generator → Fetcher_NCBI
# =================================================
# Tests the complete pipeline from natural language to NCBI data retrieval
#

set -e

echo "=========================================="
echo "INTEGRATION TEST: Query Generator → Fetcher"
echo "=========================================="

# Configuration
TEST_PHRASE="arabidopsis drought stress rna-seq"
MAX_RESULTS=10
OUTPUT_FILE="data/integration_test_results.json"

# Step 1: Generate query
echo ""
echo "Step 1: Generating NCBI query from natural language..."
echo "  Input: \"$TEST_PHRASE\""
echo ""

cd /home/lahumada/disco1/Pyner_PGRLAB/Query_generator/phases/phase3

# Run query generator in quick mode and extract final query
QUERY=$(echo "$TEST_PHRASE" | python api/main.py --quick 2>&1 | grep -A1 "✅ NCBI Query:" | tail -1 | xargs)

if [ -z "$QUERY" ]; then
    echo "❌ Failed to generate query"
    exit 1
fi

echo "✅ Generated query (truncated):"
echo "  ${QUERY:0:150}..."
echo ""

# Step 2: Fetch data from NCBI
echo "Step 2: Fetching data from NCBI SRA..."
echo "  Max Results: $MAX_RESULTS"
echo ""

cd /home/lahumada/disco1/Pyner_PGRLAB/Fetcher_NCBI

python main.py -q "$QUERY" -m "$MAX_RESULTS" -o "$OUTPUT_FILE"

# Step 3: Verify results
echo ""
echo "Step 3: Verifying results..."
echo ""

if [ -f "$OUTPUT_FILE" ]; then
    RESULT_COUNT=$(cat "$OUTPUT_FILE" | python -c "import sys, json; data=json.load(sys.stdin); print(data['metadata']['total_results'])")
    echo "✅ Integration test PASSED"
    echo "  Results saved: $OUTPUT_FILE"
    echo "  Records fetched: $RESULT_COUNT"
    echo ""
    echo "Sample result:"
    cat "$OUTPUT_FILE" | python -c "import sys, json; data=json.load(sys.stdin); r=data['results'][0] if data['results'] else {}; print(f\"  BioProject: {r.get('bioproject', 'N/A')}\"); print(f\"  Organism: {r.get('organism', 'N/A')}\"); print(f\"  Strategy: {r.get('library_strategy', 'N/A')}\"); print(f\"  Title: {r.get('title', 'N/A')[:60]}...\")"
else
    echo "❌ Integration test FAILED"
    echo "  Output file not created"
    exit 1
fi

echo ""
echo "=========================================="
echo "INTEGRATION TEST COMPLETE"
echo "=========================================="
