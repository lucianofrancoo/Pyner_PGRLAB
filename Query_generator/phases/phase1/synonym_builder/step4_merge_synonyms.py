#!/usr/bin/env python3
"""
Step 4: Merge Synonym Dictionaries

Combines synonyms from 3 sources into a unified dictionary:
  1. Step 1: MeSH official synonyms (curated by NLM)
  2. Step 2B: Text variants from abstracts (real-world usage)
  3. Step 3: PubMedBERT semantic synonyms (contextual similarity)

Output: final_synonym_dictionary.json
  - Optimized for O(1) lookup in production
  - Includes confidence scores for ranking
  - Ready to use in Query Generator
"""
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# Paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / "output"

# Input files
MESH_SYNONYMS = OUTPUT_DIR / "step1_mesh_synonyms.json"
TEXT_VARIANTS = OUTPUT_DIR / "step2_text_variants_1000s.json"
PUBMEDBERT_SYNONYMS = OUTPUT_DIR / "step3_pubmedbert_synonyms.json"

# Output file
FINAL_DICTIONARY = OUTPUT_DIR / "final_synonym_dictionary.json"

def load_json(filepath):
    """Load JSON file"""
    print(f"Loading {filepath.name}...")
    with open(filepath) as f:
        data = json.load(f)
    return data

def merge_synonyms():
    """
    Merge all synonym sources into unified dictionary

    Structure:
    {
      "drought": {
        "preferred_term": "Droughts",
        "mesh_ui": "D055864",
        "synonyms": [
          {"term": "Drought", "source": "mesh", "confidence": "official"},
          {"term": "water deficit", "source": "text", "confidence": 0.85, "frequency": 45},
          {"term": "Desiccation", "source": "pubmedbert", "confidence": 0.9377}
        ],
        "total_synonyms": 15
      }
    }
    """
    print("="*80)
    print("STEP 4: MERGING SYNONYM DICTIONARIES")
    print("="*80)
    print()

    # Load all sources
    mesh_data = load_json(MESH_SYNONYMS)
    text_variants = load_json(TEXT_VARIANTS)

    # Step 3 might not exist yet
    if PUBMEDBERT_SYNONYMS.exists():
        pubmedbert_data = load_json(PUBMEDBERT_SYNONYMS)
        # Extract just the synonyms dict (skip metadata)
        if "synonyms" in pubmedbert_data:
            pubmedbert_synonyms = pubmedbert_data["synonyms"]
        else:
            pubmedbert_synonyms = pubmedbert_data
    else:
        print("⚠ Step 3 (PubMedBERT) not found, skipping...")
        pubmedbert_synonyms = {}

    print()
    print(f"✓ MeSH terms: {len(mesh_data):,}")
    print(f"✓ Text variant terms: {len(text_variants):,}")
    print(f"✓ PubMedBERT terms: {len(pubmedbert_synonyms):,}")
    print()

    # Build unified dictionary
    print("Merging synonyms...")
    unified = {}
    stats = {
        "total_terms": 0,
        "mesh_only": 0,
        "mesh_plus_text": 0,
        "mesh_plus_pubmedbert": 0,
        "all_three_sources": 0,
        "total_synonyms": 0
    }

    for term_key, entry in mesh_data.items():
        preferred_term = entry["preferred_term"]
        mesh_ui = entry["mesh_ui"]

        # Start with MeSH official synonyms
        synonyms = []

        # Add MeSH synonyms (highest confidence)
        for syn in entry["synonyms"]:
            synonyms.append({
                "term": syn,
                "source": "mesh",
                "confidence": "official"
            })

        # Track which sources we use
        has_text = False
        has_pubmedbert = False

        # Add text variants
        if term_key in text_variants:
            has_text = True
            for variant in text_variants[term_key].get("variants", []):
                # Skip if already in MeSH synonyms (case-insensitive)
                mesh_syns_lower = set([s.lower() for s in entry["synonyms"]])
                if variant["variant"].lower() not in mesh_syns_lower and \
                   variant["variant"].lower() != preferred_term.lower():
                    synonyms.append({
                        "term": variant["variant"],
                        "source": "text",
                        "confidence": round(variant["confidence"], 4),
                        "frequency": variant["frequency"]
                    })

        # Add PubMedBERT semantic synonyms
        if term_key in pubmedbert_synonyms:
            has_pubmedbert = True
            for syn in pubmedbert_synonyms[term_key].get("pubmedbert_synonyms", []):
                # Skip if already added (case-insensitive check)
                existing_terms_lower = set([s["term"].lower() for s in synonyms])
                if syn["term"].lower() not in existing_terms_lower and \
                   syn["term"].lower() != preferred_term.lower():
                    synonyms.append({
                        "term": syn["term"],
                        "source": "pubmedbert",
                        "confidence": syn["confidence"]
                    })

        # Store unified entry
        unified[term_key] = {
            "preferred_term": preferred_term,
            "mesh_ui": mesh_ui,
            "synonyms": synonyms,
            "total_synonyms": len(synonyms)
        }

        # Update statistics
        stats["total_terms"] += 1
        stats["total_synonyms"] += len(synonyms)

        if has_text and has_pubmedbert:
            stats["all_three_sources"] += 1
        elif has_text:
            stats["mesh_plus_text"] += 1
        elif has_pubmedbert:
            stats["mesh_plus_pubmedbert"] += 1
        else:
            stats["mesh_only"] += 1

    print(f"✓ Merged {stats['total_terms']:,} terms")
    print()

    # Print statistics
    print("="*80)
    print("STATISTICS")
    print("="*80)
    print(f"Total terms: {stats['total_terms']:,}")
    print(f"Total synonyms: {stats['total_synonyms']:,}")
    print(f"Avg synonyms/term: {stats['total_synonyms']/stats['total_terms']:.2f}")
    print()
    print("Source distribution:")
    print(f"  MeSH only: {stats['mesh_only']:,} ({stats['mesh_only']/stats['total_terms']*100:.1f}%)")
    print(f"  MeSH + Text: {stats['mesh_plus_text']:,} ({stats['mesh_plus_text']/stats['total_terms']*100:.1f}%)")
    print(f"  MeSH + PubMedBERT: {stats['mesh_plus_pubmedbert']:,} ({stats['mesh_plus_pubmedbert']/stats['total_terms']*100:.1f}%)")
    print(f"  All 3 sources: {stats['all_three_sources']:,} ({stats['all_three_sources']/stats['total_terms']*100:.1f}%)")
    print()

    return unified, stats

def save_dictionary(unified, stats):
    """Save final dictionary to JSON"""
    output_data = {
        "metadata": {
            "created_at": datetime.now().isoformat(),
            "version": "1.0",
            "description": "Unified synonym dictionary from MeSH, text variants, and PubMedBERT",
            "sources": {
                "mesh": str(MESH_SYNONYMS.name),
                "text_variants": str(TEXT_VARIANTS.name),
                "pubmedbert": str(PUBMEDBERT_SYNONYMS.name) if PUBMEDBERT_SYNONYMS.exists() else None
            },
            "statistics": stats
        },
        "dictionary": unified
    }

    print(f"Saving final dictionary to {FINAL_DICTIONARY.name}...")
    with open(FINAL_DICTIONARY, 'w') as f:
        json.dump(output_data, f, indent=2)

    file_size_mb = FINAL_DICTIONARY.stat().st_size / (1024 * 1024)
    print(f"✓ Saved: {FINAL_DICTIONARY} ({file_size_mb:.1f} MB)")

def show_examples(unified, n=5):
    """Show example entries"""
    print()
    print("="*80)
    print(f"EXAMPLE ENTRIES (first {n} terms with multiple sources)")
    print("="*80)
    print()

    count = 0
    for term_key, entry in unified.items():
        # Show entries with multiple sources
        sources = set([s["source"] for s in entry["synonyms"]])
        if len(sources) > 1:
            print(f"📌 {entry['preferred_term']} ({entry['mesh_ui']})")
            print(f"   Total synonyms: {entry['total_synonyms']}")

            # Group by source
            by_source = defaultdict(list)
            for syn in entry["synonyms"]:
                by_source[syn["source"]].append(syn)

            # Show each source
            for source in ["mesh", "text", "pubmedbert"]:
                if source in by_source:
                    print(f"   {source.upper()}:")
                    for syn in by_source[source][:3]:  # Max 3 per source
                        conf = syn.get("confidence", "official")
                        freq = f", freq={syn['frequency']}" if "frequency" in syn else ""
                        print(f"      • {syn['term']} (confidence={conf}{freq})")
            print()

            count += 1
            if count >= n:
                break

def main():
    print("="*80)
    print("STEP 4: MERGE SYNONYM DICTIONARIES")
    print("="*80)
    print()
    print("This script combines:")
    print("  1. MeSH official synonyms (NLM curated)")
    print("  2. Text variants from 26M PubMed abstracts")
    print("  3. PubMedBERT semantic synonyms (optional)")
    print()
    print("Output: final_synonym_dictionary.json")
    print("  → Ready for O(1) lookup in Query Generator")
    print()
    print("-"*80)
    print()

    # Merge all sources
    unified, stats = merge_synonyms()

    # Save final dictionary
    save_dictionary(unified, stats)

    # Show examples
    show_examples(unified, n=5)

    print()
    print("="*80)
    print("DONE!")
    print("="*80)
    print()
    print("Next steps:")
    print("  1. Review: output/final_synonym_dictionary.json")
    print("  2. Integrate with Query Generator (phase3/api/)")
    print("  3. Test query expansion with real user queries")
    print()
    print("Example usage in code:")
    print("""
    import json
    with open('final_synonym_dictionary.json') as f:
        data = json.load(f)
        synonyms = data['dictionary']

    # Lookup synonyms
    term = "drought"
    if term in synonyms:
        print(synonyms[term]['synonyms'])
    """)

if __name__ == "__main__":
    main()
