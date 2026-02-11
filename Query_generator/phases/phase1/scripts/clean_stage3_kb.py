"""
Clean Stage 3 KB
- Collapse organism names to genus + species
- Normalize whitespace
- Preserve strategies counts per organism
- Tag hybrids
"""

import json
import re
from pathlib import Path
from collections import defaultdict
from datetime import datetime

INPUT_FILE = Path("/home/lahumada/disco1/Pyner_PGRLAB/Query_generator/phases/phase1/output/stage3_kb_full.json")
OUTPUT_FILE = Path("/home/lahumada/disco1/Pyner_PGRLAB/Query_generator/phases/phase1/output/stage3_kb_reduced.json")


def normalize_name(name: str) -> str:
    return re.sub(r"\s+", " ", name.strip())


def collapse_to_scientific(name: str) -> str:
    name = normalize_name(name)
    parts = name.split(" ")
    if len(parts) < 2:
        return name

    # Handle "Candidatus Genus species"
    if parts[0].lower() == "candidatus" and len(parts) >= 3:
        return " ".join(parts[:3])

    # Collapse to Genus species
    return " ".join(parts[:2])


def is_hybrid(name: str) -> bool:
    return " x " in name


def main() -> None:
    with INPUT_FILE.open() as f:
        data = json.load(f)

    organisms = data.get("organisms", {})
    organism_strategies = data.get("organism_strategies", {})
    strategies = data.get("strategies", [])

    collapsed_counts = defaultdict(int)
    collapsed_strategies = defaultdict(lambda: defaultdict(int))
    organism_tags = defaultdict(set)

    for org_name, count in organisms.items():
        clean_name = collapse_to_scientific(org_name)
        collapsed_counts[clean_name] += count

        if is_hybrid(org_name):
            organism_tags[clean_name].add("hybrid")

        # Merge strategy counts
        org_strats = organism_strategies.get(org_name, {})
        for strat, strat_count in org_strats.items():
            collapsed_strategies[clean_name][strat] += strat_count

    output = {
        "stage": data.get("stage", 3),
        "timestamp": datetime.now().isoformat(),
        "source_file": str(INPUT_FILE),
        "statistics": {
            "unique_organisms_before": len(organisms),
            "unique_organisms_after": len(collapsed_counts),
            "unique_strategies": len(strategies) if isinstance(strategies, list) else len(strategies.keys()),
        },
        "organisms": dict(collapsed_counts),
        "strategies": strategies,
        "organism_strategies": {
            org: dict(strats) for org, strats in collapsed_strategies.items()
        },
        "organism_tags": {
            org: sorted(list(tags)) for org, tags in organism_tags.items()
        },
    }

    with OUTPUT_FILE.open("w") as f:
        json.dump(output, f, indent=2)

    print(f"Saved cleaned KB to: {OUTPUT_FILE}")
    print(f"Unique organisms before: {len(organisms):,}")
    print(f"Unique organisms after:  {len(collapsed_counts):,}")


if __name__ == "__main__":
    main()
