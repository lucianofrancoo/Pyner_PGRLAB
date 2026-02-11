"""
Simple Stage 3 processor - sin deadlocks
Procesa todos los archivos de forma secuencial pero eficiente
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
print("🔬 PYNER PHASE 1 - STAGE 3 SIMPLE", flush=True)
print("="*80, flush=True)
print(f"📁 Source: {NCBI_SRA_PATH}", flush=True)
print(f"📁 Output: {OUTPUT_DIR}", flush=True)
print("="*80, flush=True)

# Estadísticas
stats = {
    "organisms": defaultdict(lambda: {"count": 0, "studies": set()}),
    "strategies": defaultdict(int),
    "sources": defaultdict(int),
    "selections": defaultdict(int),
    "organism_strategies": defaultdict(lambda: defaultdict(int)),  # organism -> strategy -> count
    "total_experiments": 0,
    "total_samples": 0,
    "total_runs": 0,
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
    """Process one bioproject and cross-reference organisms with strategies"""
    bioproject_id = bioproject_dir.name
    
    # First pass: collect strategies from this bioproject
    bioproject_strategies = set()
    exp_file = bioproject_dir / f"{bioproject_id}.experiment.xml"
    if exp_file.exists():
        root = parse_xml_safe(exp_file)
        if root is not None:
            for exp in root.findall('.//EXPERIMENT'):
                stats["total_experiments"] += 1
                
                # Strategy
                strategy = exp.find('.//LIBRARY_STRATEGY')
                if strategy is not None and strategy.text:
                    stats["strategies"][strategy.text] += 1
                    bioproject_strategies.add(strategy.text)
                
                # Source
                source = exp.find('.//LIBRARY_SOURCE')
                if source is not None and source.text:
                    stats["sources"][source.text] += 1
                
                # Selection
                selection = exp.find('.//LIBRARY_SELECTION')
                if selection is not None and selection.text:
                    stats["selections"][selection.text] += 1
    
    # Second pass: collect organisms and cross-reference with strategies
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
                    stats["organisms"][org_name]["count"] += 1
                    stats["organisms"][org_name]["studies"].add(bioproject_id)
                    
                    # Cross-reference: link this organism with all strategies in this bioproject
                    for strategy in bioproject_strategies:
                        stats["organism_strategies"][org_name][strategy] += 1
    
    # Run XML
    run_file = bioproject_dir / f"{bioproject_id}.run.xml"
    if run_file.exists():
        root = parse_xml_safe(run_file)
        if root is not None:
            runs = root.findall('.//RUN')
            stats["total_runs"] += len(runs)

# Main loop
print("\n🚀 Starting processing...", flush=True)
print(f"📁 Scanning {NCBI_SRA_PATH}...", flush=True)
start_time = time.time()

# Use generator to avoid loading all 6M dirs in memory
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
                  f"Strategies: {len(stats['strategies'])} | "
                  f"Experiments: {stats['total_experiments']:,}", flush=True)
    
    except Exception as e:
        stats["errors"] += 1
        if i % 10000 == 0:  # Print errors sample
            print(f"⚠️  Error in {bioproject_dir.name}: {e}", flush=True)

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
        "total_experiments": stats["total_experiments"],
        "total_samples": stats["total_samples"],
        "total_runs": stats["total_runs"],
        "unique_organisms": len(stats["organisms"]),
        "unique_strategies": len(stats["strategies"]),
        "unique_sources": len(stats["sources"]),
        "unique_selections": len(stats["selections"]),
        "errors": stats["errors"]
    },
    "organisms": {
        k: {"count": v["count"], "studies": sorted(list(v["studies"]))}
        for k, v in stats["organisms"].items()
    },
    "strategies": dict(stats["strategies"]),
    "sources": dict(stats["sources"]),
    "selections": dict(stats["selections"]),
    "organism_strategies": {
        org: dict(strategies)
        for org, strategies in stats["organism_strategies"].items()
    }
}

output_file = OUTPUT_DIR / "stage3_indices_full.json"
with open(output_file, 'w') as f:
    json.dump(output_data, f, indent=2)

print(f"✅ Saved to: {output_file}", flush=True)
print("\n📊 FINAL STATS:", flush=True)
print(f"   Organisms: {len(stats['organisms'])}", flush=True)
print(f"   Strategies: {len(stats['strategies'])}", flush=True)
print(f"   Experiments: {stats['total_experiments']:,}", flush=True)
print(f"   Samples: {stats['total_samples']:,}", flush=True)
print(f"   Runs: {stats['total_runs']:,}", flush=True)
print("\n" + "="*80, flush=True)
