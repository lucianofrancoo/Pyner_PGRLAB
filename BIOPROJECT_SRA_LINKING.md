# BioProject to SRA/BioSample Linking - Technical Analysis

**Date:** 2026-02-12  
**Conclusion:** Limited direct API access, but solutions exist

## What I Tested

### ✅ What Works
- ✅ Access BioProject metadata (title, description, accession, registration date)
- ✅ Get full XML of BioProject including metadata
- ✅ Extract BioProject ID and accession code

### ❌ What Doesn't Work (NCBI API Limitation)
- ❌ Direct search `PRJNA######[BioProject]` in SRA database
- ❌ Direct search `PRJNA######[BioProject]` in BioSample database
- ❌ BioProject XML does not contain embedded SRA/BioSample IDs
- ❌ No reverse linkage available through API

## Architecture Understanding

```
NCBI Database Structure:
├── SRA Database
│   ├── SRA Experiment Records (e.g., SRX###)
│   │   ├── References → BioProject (one-way link)
│   │   ├── References → BioSample
│   │   └── Contains run data
│   
├── BioSample Database
│   ├── BioSample Records (SAMN###)
│   │   ├── References → BioProject (one-way link)
│   │   └── Contains sample metadata
│   
└── BioProject Database
    ├── BioProject Records (PRJNA###)
    │   ├── NO direct links OUT to SRA/BioSample
    │   └── Only receives references FROM other databases
```

## Available Solutions

### Option 1: Use BioProject Web Programmatic Approach
Scrape NCBI BioProject web page HTML which DOES contain the counts:
```
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1418118
→ Page contains: "Studies: X", "Sample Metadata: Y", etc.
```
**Pros:** Works, has the data  
**Cons:** Web scraping (fragile, slower)

### Option 2: Use NCBI SRA Toolkit Locally
Download $NCBI_SRA/sra-stat tool to query local cache:
```bash
sra-stat --xml PRJNA1418118
```
**Pros:** Fast, reliable  
**Cons:** Requires local tools, setup

### Option 3: Accept Current Limitations
Keep BioProject as metadata layer only:
```
CSV Output:
- bioproject → ✅ Direct from BioProject API
- title → ✅ Direct from BioProject API
- submission_date → ✅ From Registration_Date  
- description → ✅ Direct from BioProject API

Cannot include:
- sra_experiments_count → Not in API
- biosamples_count → Not in API
```
**Pros:** Simple, reliable, no extra overhead  
**Cons:** Limited information

### Option 4: Manual Query in SRA After BioProject Retrieval
User can create separate workflow:
1. Get BioProjects from your query (✅ we do this)
2. For each BioProject, user manually searches SRA:
   ```
   https://www.ncbi.nlm.nih.gov/sra?term=PRJNA1418118
   ```

## Recommendation

**Current Implementation (Option 3):**
- ✅ Keep extraction simple and reliable
- ✅ Provide BioProject-level metadata only
- ✅ Focus on what NCBI API actually provides
- ✅ Users can search SRA separately if needed

**If more data needed (Option 1, future enhancement):**
- Could add web scraping to extract study/sample counts
- Would require BeautifulSoup or Selenium
- Adds complexity and dependency

## Conclusion

NCBI's API design has BioProject as a **metadata container**, not a query interface for contained data. The linking goes SRA→BioProject, not the reverse. This is by design - allows orphaned projects and maintains referential integrity.

For now, stick with current simple CSV format. If users need SRA counts, they have options:
1. Use NCBI web search separately
2. Implement web scraping (future enhancement)
3. Use SRA Toolkit locally
