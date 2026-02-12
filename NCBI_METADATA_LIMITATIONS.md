# NCBI BioProject Metadata Limitations

**Date:** 2026-02-12  
**Issue:** Several metadata fields are NOT available from NCBI BioProject API

## Analysis of Available Fields

### Test Query
```
arabidopsis[Organism] AND drought
```

### What NCBI BioProject `esummary` Response Contains

✅ **Available Fields** (populated when data exists):
```
- Project_Acc          → BioProject accession (PRJNA######)
- Project_Title        → Study title
- Project_Description  → Detailed description
- Project_Data_Type    → Type of data (e.g., "Raw sequence reads")
- Organism_Name        → Organism name (frequently empty in response)
- Organism_Label       → Organism label (frequently empty in response)
- Organism_Strain      → Strain information (frequently empty in response)
- Registration_Date    → When registered in NCBI (e.g., "2026/02/03")
- Project_MethodType   → Sequencing/other method
- Submitter_Organization → Who submitted the project
```

❌ **NOT Available in BioProject esummary**:
```
- Submission_Date      → No specific submission date available (only Registration_Date)
- Public_Date          → Not exposed by API
- Total_Studies        → Not in BioProject metadata (would need SRA database)
- Total_Runs           → Not in BioProject metadata (would need SRA database)
- Publications         → Not in BioProject metadata
- Experimental_Design  → Not in BioProject metadata
- Cultivar             → Not in BioProject metadata
```

## Root Cause

NCBI BioProject `esummary` API returns **limited summary-level metadata**. 

The missing fields are:
1. **Either SRA-level details** (Total_Studies, Total_Runs) - available through SRA database, not BioProject
2. **Or enrichment details** (Publications, Experimental_Design) - not exposed by NCBI API
3. **Or redundant dates** (Submission_Date, Public_Date) - only Registration_Date provided

## Impact on CSV Output

Current columns and their status:

| Column | Source | Status | Value |
|--------|--------|--------|-------|
| `sra_id` | Study ID | ✅ Populated | From search results |
| `bioproject` | Project_Acc | ✅ Populated | PRJNA### |
| `title` | Project_Title | ✅ Populated | Full title |
| `organisms` | Organism_Name/Label | ⚠️  Often empty | When available in metadata |
| `strain` | Organism_Strain | ⚠️ Often empty | When available in metadata |
| `cultivar` | N/A | ❌ Empty | Not in NCBI BioProject |
| `submission_date` | Registration_Date | ⚠️ Limited | Available but named differently |
| `public_date` | N/A | ❌ Empty | Not exposed by API |
| `total_studies` | N/A | ❌ Empty | Not in BioProject (SRA only) |
| `total_runs` | N/A | ❌ Empty | Not in BioProject (SRA only) |
| `publications` | N/A | ❌ Empty | Not in NCBI BioProject API |
| `experimental_design` | N/A | ❌ Empty | Not in NCBI API response |
| `description` | Project_Description | ✅ Populated | Full description |
| `fetched_at` | System | ✅ Populated | Fetch timestamp |

## Solutions

### Current Implementation (Simple)
- Use only the fields NCBI actually provides
- Accept that organism/strain/dates are often empty in BioProject metadata
- CSV columns remain but are empty unless NCBI provided the data

### Alternative 1: Enrich with SRA Query
- Keep BioProject query for filtering
- Make additional SRA queries to get:
  - Total_Studies and Total_Runs
  - Potentially publication info
- **Trade-off:** Much slower (additional API calls per BioProject)

### Alternative 2: Accept Current Limitations
- Document which fields are/aren't populated
- Focus on title and description for analysis
- Use BioProject as index, not complete metadata repository
- **Trade-off:** Simplest, but limited metadata enrichment

### Alternative 3: Query BioProject in Different Way
- Use `efetch` instead of `esummary` for detailed XML
- Parse Project Organism data from XML attributes
- **Trade-off:** More complex parsing, may only add limited data

## Recommendation

**Status Quo:** Continue with current approach
- ✅ Clean implementation
- ✅ No performance impact
- ✅ Honest about data availability
- ✅ CSV accurately reflects what NCBI provides

The empty fields in CSV are **not a bug** - they reflect **NCBI API limitations**, not code issues. The fact that Organism_Name, Strain, etc. are often empty in the BioProject response itself suggests these fields weren't populated at submission time by many researchers.

---
**Conclusion:** The metadata fields that appear empty are NOT AVAILABLE from NCBI BioProject API. This is expected behavior.
