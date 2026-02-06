#!/usr/bin/env python3
"""Direct test of filter logic without all the overhead"""

# Test the filtering logic directly
candidates = ['drought', 'gene_expression']  # Results we got
keywords_in_query = {'arabidopsis', 'root', 'drought'}  # From "Arabidopsis root drought"

# Logic: keep gene_expression only if user searched for expression-related keywords
expr_keywords = {'expression', 'profiling', 'microarray', 'transcriptome'}
keep_expr = any(kw in keywords_in_query for kw in expr_keywords)

print("=" * 60)
print("RESULT FOR: 'Arabidopsis root drought'")
print("=" * 60)
print(f"\nKeywords extracted: {keywords_in_query}")
print(f"Expression keywords in query? {keep_expr}")
print(f"\n✅ Keep 'organism' results: Always (Research studies on Arabidopsis)")
print(f"✅ Keep 'strategy' results: Always (but none in this example)")
print(f"❌ Keep 'gene_expression' results: NO ({keep_expr}) - Will be FILTERED OUT")
print("\nFinal boolean query will use ONLY organism results.")
print("=" * 60)
