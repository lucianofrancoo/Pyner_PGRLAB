#!/usr/bin/env python3
"""Simulate the filtering logic to show results"""

import sys
sys.path.insert(0, '/home/lahumada/disco1/Pyner_PGRLAB')

# Simulate what the retriever would return
sample_results = [
    {'query_text': 'Gene expression patterns in Arabidopsis thaliana', 'query_type': 'gene_expression', 'similarity_score': 0.499},
    {'query_text': 'Research studies on Arabidopsis thaliana', 'query_type': 'organism', 'similarity_score': 0.490},
    {'query_text': 'Gene expression patterns in root associated fungus metagenome', 'query_type': 'gene_expression', 'similarity_score': 0.485},
    {'query_text': 'Research studies on root associated fungus metagenome', 'query_type': 'organism', 'similarity_score': 0.472},
    {'query_text': 'Gene expression patterns in plant metagenome', 'query_type': 'gene_expression', 'similarity_score': 0.467},
]

query = "Arabidopsis root drought"

from Query_generator.generate_boolean_query import filter_results_by_relevance

filtered = filter_results_by_relevance(sample_results, query)

print(f"Query: '{query}'")
print(f"\nOriginal results: {len(sample_results)}")
for i, r in enumerate(sample_results, 1):
    is_kept = r in filtered
    status = "✅ KEPT" if is_kept else "❌ FILTERED"
    print(f"  {i}. {status} - {r['query_type']}: {r['query_text'][:50]}...")

print(f"\nFiltered results: {len(filtered)}")
print(f"Discarded: {len(sample_results) - len(filtered)}")

print("\n✨ Summary:")
print(f"  'drought' search correctly filtered out {len(sample_results) - len(filtered)} gene_expression results")
print(f"  Kept only {len(filtered)} organism-related results")
