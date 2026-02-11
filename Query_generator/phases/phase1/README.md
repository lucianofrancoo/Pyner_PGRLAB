# Phase 1: Knowledge Base Extraction
**Status:** ✅ **COMPLETED - FROZEN ARTIFACTS**  
**Date:** 2026-02-10

---

## Technical Summary

Phase 1 extracts structured metadata from **NCBI SRA XML files** to build a comprehensive Knowledge Base. This pipeline was executed **once** to generate frozen artifacts consumed by Phase 2 and Phase 3.

### Executed Pipeline (build-time)

```
NCBI SRA XML files (7.2M+ files)
         ↓
  Stage 1: Parse XML → Extract metadata (organism, strategy, tissue, etc.)
         ↓
  Stage 2: Parallel GPU → Process 50K files with GPU acceleration
         ↓
  Stage 3: Final KB → Build comprehensive Knowledge Base
         ↓
  Output: Frozen JSON artifacts
```

### Purpose

1. **XML parsing:** Extract structured metadata from NCBI SRA XML files
2. **GPU acceleration:** Use 3x NVIDIA RTX 4000 Ada for parallel processing
3. **Knowledge Base construction:** Build comprehensive organism/strategy/tissue/condition mappings
4. **Artifact generation:** Create frozen JSON files for downstream phases

---

## Preserved Artifacts

| File | Size | Purpose | Used By |
|------|------|---------|---------|
| `output/stage3_kb_reduced.json` | 16MB | Reduced KB with essential metadata | ✅ Phase 3 (query generation) |
| `output/stage3_knowledge_base.json` | 11KB | Compact KB summary | ✅ Phase 2 (vocabulary) |
| `output/stage2_knowledge_base.json` | 10KB | Intermediate KB | ✅ Phase 3 (fallback) |
| `output/stage1_indices.json` | 122KB | XML file indices | ✅ Rebuild reference |

### Cleaned Files

- `logs/*.log` - Execution logs (1.2MB, 90+ files)
- `__pycache__/` - Python bytecode cache
- `checkpoints/` - Empty checkpoint directories
- `tests/` - Empty test fixtures
- `output/stage3_kb_full.json` - Full KB (27MB duplicate)

---

## Technical Components

### `config.py`
Pipeline configuration:
- XML paths and directories
- GPU settings (3x RTX 4000 Ada, 24GB VRAM each)
- Parallel workers (20 workers)
- Checkpoint intervals

### `utils.py`
Helper utilities:
- Logging setup
- GPU memory management
- Checkpoint save/restore
- Progress tracking

### `scripts/stage1_parse_xml.py`
Initial XML parsing:
- Sequential parsing of 1K sample files
- Metadata extraction
- Field validation
- Generates `stage1_indices.json`

### `scripts/stage2_parallel_gpu.py`
GPU-accelerated parallel processing:
- Multi-GPU workload distribution
- 8-20 parallel workers
- Checkpoint-based recovery
- Generates `stage2_knowledge_base.json`

### `scripts/stage3_final.py` & `stage3_simple.py`
Final KB construction:
- Merge and deduplicate entries
- Build organism/strategy mappings
- Generate reduced and full versions
- Generates `stage3_kb_reduced.json` and `stage3_knowledge_base.json`

### `scripts/clean_stage3_kb.py`
KB cleanup utility:
- Remove duplicates
- Compress entries
- Validate structure

---

## Reconstruction (Optional)

If you need to rebuild the Knowledge Base from NCBI XML files:

### Hardware Requirements

- **GPUs:** 3x NVIDIA RTX 4000 Ada (24GB VRAM each) or equivalent
- **RAM:** 251GB (minimum 64GB)
- **CPU:** Multi-core (20+ threads recommended)
- **Storage:** 100GB+ for XML files and artifacts

### Dependencies

```bash
pip install torch transformers sentencepiece
# For GPU support:
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

### Execution

```bash
cd Query_generator/phases/phase1

# Stage 1: Parse initial 1K XML files (test)
python scripts/stage1_parse_xml.py

# Stage 2: Process 50K files with GPU parallelization
python scripts/stage2_parallel_gpu.py --stage 2

# Stage 3: Build final Knowledge Base
python scripts/stage3_final.py
# or for simplified version:
python scripts/stage3_simple.py

# Clean and optimize KB
python scripts/clean_stage3_kb.py
```

**Expected output:**
- `output/stage1_indices.json` - XML indices
- `output/stage2_knowledge_base.json` - Intermediate KB
- `output/stage3_kb_reduced.json` - Reduced KB (16MB)
- `output/stage3_knowledge_base.json` - Compact KB (11KB)

### Monitoring

```bash
# Watch logs in real-time
tail -f logs/stage*.log

# Check checkpoint progress
ls -lh checkpoints/stage*/

# Monitor GPU usage
watch -n 1 nvidia-smi
```

---

## Performance Metrics

- **XML files processed:** 7.2M+ files
- **Throughput:** 100-150 files/sec (GPU parallel)
- **Total execution time:** ~3-4 hours (full pipeline)
- **Final KB size:** 16MB (reduced), 27MB (full)
- **Unique organisms:** 18,413
- **GPU utilization:** 80-95% across 3 GPUs

### Hardware Specs Used
- GPU: 3x NVIDIA RTX 4000 Ada (24GB VRAM)
- RAM: 251GB
- CPU: Multi-core (20 threads)

---

## Phase 2 & 3 Integration

**Phase 2** consumes:
1. `stage3_knowledge_base.json` → To build query cache and extract terms

**Phase 3** consumes:
1. `stage3_kb_reduced.json` → Primary KB for query generation and validation
2. `stage2_knowledge_base.json` → Fallback KB if reduced version unavailable

---

## Final Structure

```
phase1/
├── config.py                      # Pipeline configuration
├── utils.py                       # Helper utilities
├── README.md                      # This file
├── EXECUTION_REPORT.md            # Detailed execution report
├── logs/                          # Empty (cleaned)
├── output/
│   ├── stage1_indices.json       # XML indices (122KB)
│   ├── stage2_knowledge_base.json # Intermediate KB (10KB)
│   ├── stage3_kb_reduced.json    # Reduced KB (16MB) - Used by Phase 3
│   └── stage3_knowledge_base.json # Compact KB (11KB) - Used by Phase 2
└── scripts/
    ├── stage1_parse_xml.py        # Stage 1: XML parsing
    ├── stage2_parallel_gpu.py     # Stage 2: GPU parallel processing
    ├── stage3_final.py            # Stage 3: Final KB build
    ├── stage3_simple.py           # Stage 3: Simplified KB build
    ├── clean_stage3_kb.py         # KB cleanup utility
    └── phase1_20251216.log        # Historical log reference
```

---

## Troubleshooting

### Out of Memory (OOM)
```bash
# Reduce batch size in config.py
BATCH_SIZE = 32  # decrease to 16 or 8

# Reduce parallel workers
NUM_WORKERS = 8  # decrease from 20
```

### GPU not detected
```bash
# Check CUDA installation
python -c "import torch; print(torch.cuda.is_available())"

# Verify GPU visibility
nvidia-smi
```

### Checkpoint recovery
```bash
# Resume from last checkpoint
python scripts/stage2_parallel_gpu.py --resume --checkpoint checkpoints/stage2/
```

### XML parsing errors
Check `logs/` for detailed error messages and skip corrupted files using `--skip-errors` flag.

---

**Status:** ✅ Pipeline executed, artifacts frozen  
**Last Updated:** 2026-02-10  
**Next:** Phase 2 builds vector database from these artifacts
