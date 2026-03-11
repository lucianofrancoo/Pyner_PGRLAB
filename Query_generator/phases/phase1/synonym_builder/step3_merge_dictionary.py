#!/usr/bin/env python3
"""
Step 3: Merge MeSH + Text Variants into Final Dictionary

Combines:
  1. Step 1: MeSH official synonyms (curated by NLM)
  2. Step 2B: Text variants from 26M PubMed abstracts

Output: final_synonym_dictionary.json
  - Optimized for O(1) lookup in production
  - Used by query_expander.py for exact matches
  - PubMedBERT used only for on-the-fly semantic fallback

Note: This replaces the old step4 that included PubMedBERT pre-computed synonyms.
      Now PubMedBERT is used only as fallback for unknown terms.
"""
import json
from pathlib import Path
from datetime import datetime

# Paths
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / "output"

# Input files
MESH_SYNONYMS = OUTPUT_DIR / "step1_mesh_synonyms.json"
TEXT_VARIANTS = OUTPUT_DIR / "step2_text_variants_1000s.json"

# Output file
FINAL_DICTIONARY = OUTPUT_DIR / "final_synonym_dictionary.json"

def load_json(filepath):
    """Load JSON file"""
    print(f"Loading {filepath.name}...")
    with open(filepath) as f:
        data = json.load(f)
    print(f"✓ Loaded {len(data):,} entries")
    return data

def merge_mesh_and_text():
    """
    Merge MeSH synonyms with text variants

    Structure:
    {
      "drought": {
        "preferred_term": "Droughts",
        "mesh_ui": "D055864",
        "synonyms": [
          {"term": "Drought", "source": "mesh", "confidence": "official"},
          {"term": "water deficit", "source": "text", "confidence": 0.85, "frequency": 45}
        ]
      }
    }
    """
    print("\n" + "="*80)
    print("MERGING MeSH + TEXT VARIANTS")
    print("="*80)
    print()

    # Load sources
    mesh_data = load_json(MESH_SYNONYMS)
    text_variants = load_json(TEXT_VARIANTS)
    print()

    # Build unified dictionary
    print("Merging synonyms...")
    unified = {}
    stats = {
        "total_terms": len(mesh_data),
        "mesh_only": 0,
        "mesh_plus_text": 0,
        "total_synonyms": 0,
        "total_mesh_synonyms": 0,
        "total_text_variants": 0
    }

    for term_key, entry in mesh_data.items():
        preferred_term = entry["preferred_term"]
        mesh_ui = entry["mesh_ui"]

        # Start with MeSH official synonyms
        synonyms = []

        # Add MeSH synonyms
        for syn in entry["synonyms"]:
            synonyms.append({
                "term": syn,
                "source": "mesh",
                "confidence": "official"
            })

        mesh_count = len(synonyms)
        stats["total_mesh_synonyms"] += mesh_count

        # Add text variants (if exist)
        text_count = 0
        if term_key in text_variants:
            mesh_syns_lower = set([s.lower() for s in entry["synonyms"]])

            # text_variants structure: {"variant_term": {"count": N, "percentage": X, "confidence": "high"}}
            for variant_term, variant_data in text_variants[term_key].get("text_variants", {}).items():
                variant_lower = variant_term.lower()

                # Skip if already in MeSH synonyms
                if variant_lower in mesh_syns_lower or variant_lower == preferred_term.lower():
                    continue

                # Convert confidence level to numeric
                confidence_map = {
                    "very_high": 0.95,
                    "high": 0.85,
                    "medium": 0.75,
                    "low": 0.65
                }
                confidence_str = variant_data.get("confidence", "medium")
                confidence_numeric = confidence_map.get(confidence_str, 0.75)

                synonyms.append({
                    "term": variant_term,
                    "source": "text",
                    "confidence": confidence_numeric,
                    "frequency": variant_data.get("count", 0),
                    "percentage": variant_data.get("percentage", 0)
                })
                text_count += 1

        stats["total_text_variants"] += text_count

        # Store unified entry
        unified[term_key] = {
            "preferred_term": preferred_term,
            "mesh_ui": mesh_ui,
            "synonyms": synonyms,
            "synonym_count": {
                "total": len(synonyms),
                "mesh": mesh_count,
                "text": text_count
            }
        }

        # Update statistics
        stats["total_synonyms"] += len(synonyms)

        if text_count > 0:
            stats["mesh_plus_text"] += 1
        else:
            stats["mesh_only"] += 1

    print(f"✓ Merged {stats['total_terms']:,} terms")
    print()

    return unified, stats

def save_dictionary(unified, stats):
    """Save final dictionary to JSON"""
    output_data = {
        "metadata": {
            "created_at": datetime.now().isoformat(),
            "version": "2.0",
            "description": "Unified synonym dictionary from MeSH official + text variants",
            "strategy": "Pre-computed dictionary for exact matches, PubMedBERT fallback for unknown terms",
            "sources": {
                "mesh": str(MESH_SYNONYMS.name),
                "text_variants": str(TEXT_VARIANTS.name)
            },
            "statistics": stats,
            "usage": {
                "exact_match": "O(1) lookup in this dictionary",
                "semantic_fallback": "Use query_expander.py with PubMedBERT for unknown terms"
            }
        },
        "dictionary": unified
    }

    print("Saving final dictionary...")
    with open(FINAL_DICTIONARY, 'w') as f:
        json.dump(output_data, f, indent=2)

    file_size_mb = FINAL_DICTIONARY.stat().st_size / (1024 * 1024)
    print(f"✓ Saved: {FINAL_DICTIONARY.name} ({file_size_mb:.1f} MB)")
    print(f"  Location: {FINAL_DICTIONARY}")

def print_statistics(stats):
    """Print detailed statistics"""
    print("\n" + "="*80)
    print("STATISTICS")
    print("="*80)
    print(f"Total MeSH terms: {stats['total_terms']:,}")
    print(f"Total synonyms: {stats['total_synonyms']:,}")
    print(f"  → MeSH official: {stats['total_mesh_synonyms']:,}")
    print(f"  → Text variants: {stats['total_text_variants']:,}")
    print()
    print(f"Avg synonyms per term: {stats['total_synonyms']/stats['total_terms']:.2f}")
    print()
    print("Enrichment:")
    print(f"  MeSH only: {stats['mesh_only']:,} ({stats['mesh_only']/stats['total_terms']*100:.1f}%)")
    print(f"  MeSH + Text variants: {stats['mesh_plus_text']:,} ({stats['mesh_plus_text']/stats['total_terms']*100:.1f}%)")
    print()

def show_examples(unified, n=10):
    """Show example entries"""
    print("="*80)
    print(f"EXAMPLE ENTRIES (first {n} terms with text variants)")
    print("="*80)
    print()

    count = 0
    for term_key, entry in unified.items():
        if entry["synonym_count"]["text"] > 0:
            print(f"📌 {entry['preferred_term']} ({entry['mesh_ui']})")
            print(f"   Total: {entry['synonym_count']['total']} synonyms " +
                  f"({entry['synonym_count']['mesh']} MeSH + {entry['synonym_count']['text']} text)")

            # Show MeSH
            mesh_syns = [s for s in entry["synonyms"] if s["source"] == "mesh"]
            if mesh_syns:
                print(f"   MeSH: {', '.join([s['term'] for s in mesh_syns[:3]])}")

            # Show text
            text_syns = [s for s in entry["synonyms"] if s["source"] == "text"]
            if text_syns:
                print(f"   Text variants:")
                for syn in text_syns[:5]:
                    print(f"      • {syn['term']} (conf={syn['confidence']}, freq={syn['frequency']})")

            print()
            count += 1
            if count >= n:
                break

def main():
    print("="*80)
    print("STEP 3: MERGE MeSH + TEXT VARIANTS")
    print("="*80)
    print()
    print("This creates a pre-computed dictionary for fast lookups.")
    print("PubMedBERT is used separately for semantic fallback on unknown terms.")
    print()

    # Merge
    unified, stats = merge_mesh_and_text()

    # Save
    save_dictionary(unified, stats)

    # Statistics
    print_statistics(stats)

    # Examples
    show_examples(unified, n=10)

    print()
    print("="*80)
    print("DONE!")
    print("="*80)
    print()
    print("Next steps:")
    print("  1. Test with: python3 query_expander.py 'drought' 'arabidopsis'")
    print("  2. Test fallback: python3 query_expander.py 'water scarcity'")
    print("  3. Integrate with Phase 3 Query Generator")
    print()

if __name__ == "__main__":
    main()
