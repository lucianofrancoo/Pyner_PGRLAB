#!/usr/bin/env python3
"""Quick test of refined strategy synonyms"""

import sys
sys.path.insert(0, '.')

from Query_generator.generate_boolean_query import (
    RNASEQ_SYNS, GENE_EXPR_SYNS, CHIP_SYNS, WGS_SYNS, QPCR_SYNS, METAGENOMICS_SYNS
)

print("=== Strategy-specific Synonyms ===\n")
print(f"RNASEQ_SYNS:\n  {RNASEQ_SYNS}\n")
print(f"GENE_EXPR_SYNS:\n  {GENE_EXPR_SYNS}\n")
print(f"CHIP_SYNS:\n  {CHIP_SYNS}\n")
print(f"WGS_SYNS:\n  {WGS_SYNS}\n")
print(f"QPCR_SYNS:\n  {QPCR_SYNS}\n")
print(f"METAGENOMICS_SYNS:\n  {METAGENOMICS_SYNS}\n")

print("✅ All synonyms loaded successfully!")
print("\nKey improvements:")
print("  - RNA-Seq now limited to specific synonyms (not 'gene expression')")
print("  - Gene expression uses microarray, profiling (separate from sequencing)")
print("  - Each strategy has dedicated synonyms (WGS, ChIP-seq, qPCR, etc.)")
print("  - No overlap between strategies")
