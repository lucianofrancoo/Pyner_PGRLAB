#!/usr/bin/env python3
import sys
sys.path.insert(0, '/home/lahumada/disco1/Pyner_PGRLAB')

from Query_generator.generate_boolean_query import extract_query_keywords, get_relevant_query_types

# Test 1: "Arabidopsis root drought"
query1 = "Arabidopsis root drought"
keywords1 = extract_query_keywords(query1)
relevant1 = get_relevant_query_types(query1)

print(f"Query: '{query1}'")
print(f"Keywords: {keywords1}")
print(f"Relevant query types: {relevant1}")
print()

# Test 2: "Mouse RNAseq expression"
query2 = "Mouse RNAseq expression"
keywords2 = extract_query_keywords(query2)
relevant2 = get_relevant_query_types(query2)

print(f"Query: '{query2}'")
print(f"Keywords: {keywords2}")
print(f"Relevant query types: {relevant2}")
