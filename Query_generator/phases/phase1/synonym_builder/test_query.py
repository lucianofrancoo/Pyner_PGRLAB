#!/usr/bin/env python3
"""
Quick test script to find semantic synonyms for specific terms
Usage: python3 test_query.py "drought" "arabidopsis" "rna-seq"
"""
import json
import sys
from pathlib import Path
from sentence_transformers import SentenceTransformer, util
import torch

# Paths
SCRIPT_DIR = Path(__file__).parent
MESH_FILE = SCRIPT_DIR / "output" / "step1_mesh_synonyms.json"
CACHE_DIR = SCRIPT_DIR / "pubmedbert_cache"
EMBEDDINGS_CACHE = CACHE_DIR / "mesh_embeddings.pt"
TERMS_CACHE = CACHE_DIR / "mesh_terms.json"

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 test_query.py <term1> [term2] [term3] ...")
        print("\nExample:")
        print('  python3 test_query.py "drought" "arabidopsis" "rna-seq"')
        print('  python3 test_query.py "water stress" "gene expression"')
        sys.exit(1)

    query_terms = sys.argv[1:]

    print("="*80)
    print("PubMedBERT Semantic Synonym Finder - Quick Test")
    print("="*80)
    print()

    # Load model
    print("Loading PubMedBERT model...")
    model = SentenceTransformer('pritamdeka/S-PubMedBert-MS-MARCO')
    print("✓ Model loaded")
    print()

    # Load MeSH data
    print("Loading MeSH terms...")
    with open(MESH_FILE) as f:
        mesh_data = json.load(f)
    preferred_terms = [entry["preferred_term"] for entry in mesh_data.values()]
    print(f"✓ {len(preferred_terms):,} MeSH terms loaded")
    print()

    # Load or create embeddings
    if EMBEDDINGS_CACHE.exists() and TERMS_CACHE.exists():
        print("Loading cached embeddings...")
        embeddings = torch.load(EMBEDDINGS_CACHE)
        with open(TERMS_CACHE) as f:
            cached_terms = json.load(f)

        if cached_terms == preferred_terms:
            print(f"✓ Loaded from cache: {embeddings.shape}")
        else:
            print("Cache mismatch, re-encoding...")
            embeddings = model.encode(preferred_terms, show_progress_bar=True, convert_to_tensor=True)
    else:
        print("Creating embeddings (first time, ~1-2 min)...")
        embeddings = model.encode(preferred_terms, show_progress_bar=True, convert_to_tensor=True)

        # Save cache
        CACHE_DIR.mkdir(exist_ok=True)
        torch.save(embeddings, EMBEDDINGS_CACHE)
        with open(TERMS_CACHE, 'w') as f:
            json.dump(preferred_terms, f)
        print(f"✓ Cached at {CACHE_DIR}/")

    print()
    print("="*80)
    print("RESULTS")
    print("="*80)
    print()

    # Process each query
    for query in query_terms:
        print(f"🔍 Query: '{query}'")
        print()

        # Encode query
        query_embedding = model.encode(query, convert_to_tensor=True)

        # Compute similarities
        similarities = util.cos_sim(query_embedding, embeddings)[0]

        # Get top 10 matches
        top_indices = similarities.argsort(descending=True)[:10]

        # Show results
        for rank, idx in enumerate(top_indices, 1):
            term = preferred_terms[idx.item()]
            score = similarities[idx].item()
            bar = "█" * int(score * 20)

            print(f"   {rank:2d}. {term:45s} {score:.4f}  {bar}")

        print()
        print("-"*80)
        print()

if __name__ == "__main__":
    main()
