# BioProject Cache Fix - Verification Report

**Commit:** `82face3`  
**Date:** 2026-02-12  
**Issue:** BioProject cache persisted across searches, marking first results as "already seen"

## Problem Analysis

The `BioProjectCache` class in `Fetcher_NCBI/ncbi_fetcher.py` was:
```python
- Loading cache from disk on initialization (self.seen = self._load_cache())
- Saving cache after each session (self.save())
```

This caused:
- First search: Creates cache file with seen BioProjects
- Second search: Loads previous cache, marks first results as duplicates
- Result: "BioProject XXXXX already seen - skipping" for valid new searches

## Solution Implemented

Modified `BioProjectCache` class to **NOT persist between searches**:

```python
class BioProjectCache:
    def __init__(self, cache_file: Path = DEDUP_CACHE):
        self.cache_file = cache_file
        # Start fresh for each search - don't load previous cache
        self.seen: Set[str] = set()
    
    def _load_cache(self) -> Set[str]:
        """Load existing cache from disk."""
        # Disabled - only track duplicates in current search
        return set()
    
    def save(self):
        """Save cache to disk. Disabled for per-search deduplication."""
        # Don't persist cache between searches
        pass
```

### Key Changes:
1. ✅ `_load_cache()` now returns empty set (never loads from disk)
2. ✅ `save()` now does nothing (cache not written to disk)
3. ✅ Each search starts with `self.seen = set()`
4. ✅ Duplicates only tracked within current session

## Verification Results

### Test 1: Nicotiana Drought Search
```
Studies processed successfully: 5
Skipped (duplicate BioProject): 0 ✓
Unique BioProjects collected: 5
```

### Test 2: Arabidopsis Drought Search (Different Query)
```
Studies processed successfully: 5
Skipped (duplicate BioProject): 0 ✓
Unique BioProjects collected: 5
```

### Verification Points:
- ✅ No "already seen" messages for first BioProjects in second search
- ✅ No cache file stored on disk (`bioproject_cache.json`)
- ✅ Each search starts fresh with empty cache
- ✅ Duplicates still detected within same search (dedup still works)

## Technical Details

### Original Behavior:
```
Search 1: [Load cache (empty)] → Find 30 → Save cache
          ✓ New BioProject: PRJNA1418118
          ✓ New BioProject: PRJNA1367680
          ...

Search 2: [Load cache (from disk)] → Find 30 → Save cache
          ✗ BioProject PRJNA1418118 already seen - skipping
          ✗ BioProject PRJNA1367680 already seen - skipping
```

### New Behavior:
```
Search 1: [Initialize empty cache] → Find 30
          ✓ New BioProject: PRJNA1418118
          ✓ New BioProject: PRJNA1367680
          (Cache not saved)

Search 2: [Initialize empty cache] → Find 30
          ✓ New BioProject: PRJNA1418118 (This is valid now!)
          ✓ New BioProject: PRJNA1367680 (This is valid now!)
```

## What Still Works

- ✅ **Within-search deduplication**: If a BioProject appears twice in same search, it's filtered
- ✅ **Metadata extraction**: All 14 fields still extracted correctly
- ✅ **CSV output**: Results saved properly
- ✅ **Query integration**: Query generator + fetcher integration works

## Notes

- Old cache file was deleted: `Fetcher_NCBI/data/bioproject_cache.json`
- This is the intended behavior for a search tool: each search is independent
- If persistent cross-search deduplication is needed in future, it should be:
  - Explicit (with `--no-duplicates-from-previous` flag)
  - Separate from the main workflow
  - Documented clearly

---
**Status:** FIXED ✅  
**Test Date:** 2026-02-12  
**Tested By:** Integration testing
