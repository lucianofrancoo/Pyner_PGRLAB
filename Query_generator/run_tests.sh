#!/usr/bin/env bash
set -e
echo "Running Query_generator quick tests"
python3 Query_generator/quick_test.py
echo
echo "Running example test_queries"
python3 Query_generator/test_queries.py "virus humanos"

echo
echo "Done. Use --interactive for interactive mode:"
echo "  python3 Query_generator/test_queries.py --interactive"
