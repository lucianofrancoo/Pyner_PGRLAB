"""
Stage 3 - Final KB processor
Solo: Organismos + Estrategias asociadas (counts)
Sin: Listas de bioprojects, experiments detail, sources, selections
"""

import json
import xml.etree.ElementTree as ET
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import sys
import time

sys.path.insert(0, str(Path(__file__).parent.parent))
from config import NCBI_SRA_PATH, OUTPUT_DIR

print("="*80, flush=True)
print("🔬 PYNER PHASE 1 - STAGE 3 FINAL (SIMPLIFIED)", flush=True)
print("="*80, flush=True)
print(f"📁 Source: {NCBI_SRA_PATH}", flush=True)
print(f"📁 Output: {OUTPUT_DIR}", flush=True)
print("="*80, flush=True)

# Solo estadísticas necesarias
stats = {
    "organisms": defaultdict(int),  # organism_name: count
    "organism_strategies": defaultdict(lambda: defaultdict(int)),  # organism: {strategy: count}
    "unique_strategies": set(),
    "total_samples": 0,
    "errors": 0
}

def parse_xml_safe(xml_file):
    """Parse XML safely"""
    try:
        tree = ET.parse(xml_file)
        return tree.getroot()
    except:
        return None

def process_bioproject(bioproject_dir):
    """Process one bioproject - only organisms and strategies"""
    bioproject_id = bioproject_dir.name
    
    # Collect strategies from this bioproject
    bioproject_strategies = set()
    exp_file = bioproject_dir / f"{bioproject_id}.experiment.xml"
    if exp_file.exists():
        root = parse_xml_safe(exp_file)
        if root is not None:
            for exp in root.findall('.//EXPERIMENT'):
                # Strategy
                strategy = exp.find('.//LIBRARY_STRATEGY')
                if strategy is not None and strategy.text:
                    bioproject_strategies.add(strategy.text)
                    stats["unique_strategies"].add(strategy.text)
    
    # Collect organisms and link to strategies
    sample_file = bioproject_dir / f"{bioproject_id}.sample.xml"
    if sample_file.exists():
        root = parse_xml_safe(sample_file)
        if root is not None:
            for sample in root.findall('.//SAMPLE'):
                stats["total_samples"] += 1
                
                # Organism
                organism = sample.find('.//SCIENTIFIC_NAME')
                if organism is not None and organism.text:
                    org_name = organism.text
                    stats["organisms"][org_name] += 1
                    
                    # Link organism to strategies from this bioproject
                    for strategy in bioproject_strategies:
                        stats["organism_strategies"][org_name][strategy] += 1

# Main loop
print("\n🚀 Starting processing...", flush=True)
print(f"📁 Scanning {NCBI_SRA_PATH}...", flush=True)
start_time = time.time()

i = 0
for bioproject_dir in NCBI_SRA_PATH.iterdir():
    if not bioproject_dir.is_dir():
        continue
    
    i += 1
    try:
        process_bioproject(bioproject_dir)
        
        # Progress every 10K
        if i % 10000 == 0:
            elapsed = time.time() - start_time
            rate = i / elapsed
            
            print(f"✅ Progress: {i:,} files | "
                  f"Rate: {rate:.0f} files/sec | "
                  f"Elapsed: {elapsed/3600:.1f}h", flush=True)
            print(f"   Organisms: {len(stats['organisms'])} | "
                  f"Strategies: {len(stats['unique_strategies'])}", flush=True)
    
    except Exception as e:
        stats["errors"] += 1

total = i
elapsed = time.time() - start_time
print(f"\n✅ Processing complete: {total:,} files in {elapsed/3600:.2f} hours", flush=True)
print(f"   Rate: {total/elapsed:.1f} files/sec", flush=True)

# Save results
print("\n💾 Saving results...", flush=True)

output_data = {
    "stage": 3,
    "timestamp": datetime.now().isoformat(),
    "files_processed": total,
    "processing_time_hours": elapsed / 3600,
    "statistics": {
        "total_samples": stats["total_samples"],
        "unique_organisms": len(stats["organisms"]),
        "unique_strategies": len(stats["unique_strategies"]),
        "errors": stats["errors"]
    },
    "organisms": dict(stats["organisms"]),  # organism: count
    "strategies": sorted(list(stats["unique_strategies"])),
    "organism_strategies": {
        org: dict(strategies)
        for org, strategies in stats["organism_strategies"].items()
    }
}

output_file = OUTPUT_DIR / "stage3_kb_full.json"
with open(output_file, 'w') as f:
    json.dump(output_data, f, indent=2)

print(f"✅ Saved to: {output_file}", flush=True)
print("\n📊 FINAL STATS:", flush=True)
print(f"   Organisms: {len(stats['organisms']):,}", flush=True)
print(f"   Strategies: {len(stats['unique_strategies'])}", flush=True)
print(f"   Samples: {stats['total_samples']:,}", flush=True)
print(f"   Errors: {stats['errors']:,}", flush=True)
print("\n" + "="*80, flush=True)
