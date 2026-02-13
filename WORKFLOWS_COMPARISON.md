# Pyner Workflows Comparison

Quick guide to choose between the two main scripts.

---

## 📊 Quick Comparison

| Feature | `test_data_analyzer.sh` | `test_fetcher_integrator.sh` |
|---------|------------------------|------------------------------|
| **Generate Boolean Query** | ✅ Yes | ✅ Yes |
| **Search PubMed** | ✅ Yes | ✅ Yes |
| **Search BioProjects** | ✅ Yes | ✅ Yes |
| **AI Analysis** | ✅ **Yes** | ❌ No |
| **Relevance Scoring** | ✅ **Yes** | ❌ No |
| **Extract Organisms** | ✅ **Yes** | ❌ No |
| **Extract Tissues** | ✅ **Yes** | ❌ No |
| **Extract Conditions** | ✅ **Yes** | ❌ No |
| **Requires Ollama** | ✅ Yes (optional) | ❌ No |
| **Speed** | 🐢 Slower (AI analysis) | ⚡ Fast |
| **Final Output** | Classified table CSV | Raw CSV + JSON |

---

## 🎯 When to Use Each

### Use `test_data_analyzer.sh` when:

✅ You want **automatic analysis** of papers  
✅ You need **relevance scoring** for filtering  
✅ You want **structured metadata** extracted (organisms, tissues, conditions)  
✅ You're doing **systematic review** or **meta-analysis**  
✅ You have **Ollama running** (or can start it)  

**Example use case:** 
> "I need to find all papers about tomato drought stress in roots using RNA-Seq, automatically classify them by relevance, and extract all mentioned organisms and conditions for further analysis."

### Use `test_fetcher_integrator.sh` when:

✅ You just want **raw search results** quickly  
✅ You'll do **manual review** of papers  
✅ You **don't have Ollama** installed  
✅ You want to **explore** what's available first  
✅ You need **maximum speed**  

**Example use case:**
> "I want to quickly see what papers are out there about my topic, download the raw list, and review them manually later."

---

## 🔄 Complete Workflow

### Option 1: Full Automated Analysis (Recommended)

```bash
bash test_data_analyzer.sh
```

**Steps:**
1. Enter natural language query
2. Select database (PubMed/BioProject)
3. Generate boolean query with AI
4. Fetch papers from NCBI
5. **Analyze with AI (Ollama)**
6. Get classified table

**Output:**
```
Data_Analyzer/output/classified_papers_20260213_114404.csv
  ↓
PMID, Title, Relevance_Score, Is_Relevant, Organisms, Tissues, 
Conditions, Strategies, Year, Journal, DOI, Abstract_Preview
```

### Option 2: Quick Search Only

```bash
bash test_fetcher_integrator.sh
```

**Steps:**
1. Enter natural language query
2. Select database (PubMed/BioProject)
3. Generate boolean query with AI
4. Fetch papers from NCBI

**Output:**
```
pubmed_results_20260213_113559.csv
  ↓
PMID, Title, Authors, Year, Journal, Publication_Type, DOI, PMCID, URL
```

---

## 📈 Performance

### `test_data_analyzer.sh`

- **Query generation:** ~2-5 seconds
- **Fetching (PubMed):** ~1-10 seconds (depends on number of papers)
- **AI Analysis:** ~5-10 seconds per paper
- **Total for 10 papers:** ~2-3 minutes
- **Total for 100 papers:** ~10-15 minutes

### `test_fetcher_integrator.sh`

- **Query generation:** ~2-5 seconds
- **Fetching (PubMed):** ~1-10 seconds  
- **Total:** ~5-15 seconds

**Speed difference:** Analysis adds ~5-10 seconds per paper

---

## 🔧 Requirements Comparison

### `test_data_analyzer.sh`

**Required:**
- Python 3
- BioPython
- Ollama server running
- Qwen 2.5:14b model pulled

**Setup:**
```bash
# Install dependencies
pip install biopython

# Install and start Ollama
curl -fsSL https://ollama.com/install.sh | sh
ollama pull qwen2.5:14b
ollama serve  # Keep running in background
```

### `test_fetcher_integrator.sh`

**Required:**
- Python 3
- BioPython

**Setup:**
```bash
# Install dependencies
pip install biopython
```

---

## 💡 Pro Tips

### Hybrid Approach

You can run both sequentially:

```bash
# 1. Quick exploration (fast)
bash test_fetcher_integrator.sh
# Review: "Are these papers relevant? Do I need analysis?"

# 2. If yes, run full analysis on saved results
cd Data_Analyzer
python3 paper_analyzer.py ../pubmed_results_<timestamp>.json
```

### Skip Analysis in test_data_analyzer.sh

If Ollama is not available when running `test_data_analyzer.sh`, it will:
- ✅ Complete the search normally
- ⚠️ Skip the analysis step
- ℹ️ Show warning message
- ✅ Save raw results (same as test_fetcher_integrator.sh)

So you can use `test_data_analyzer.sh` even without Ollama, and it will work like `test_fetcher_integrator.sh`.

---

## 🎓 Use Case Examples

### Research Use Case: Systematic Review

**Goal:** Find all papers about drought stress in tomato roots, classify by relevance, extract all studied conditions.

**Script:** `test_data_analyzer.sh` ✅

**Why:** Automatic classification saves hours of manual review. Structured extraction makes it easy to group papers by specific conditions.

### Quick Exploration

**Goal:** "What papers exist about my topic? Let me browse them."

**Script:** `test_fetcher_integrator.sh` ✅

**Why:** Fast results, no need for Ollama, raw data ready for manual review.

### Meta-Analysis

**Goal:** Collect papers and extract all organisms, tissues, and techniques used across studies.

**Script:** `test_data_analyzer.sh` ✅

**Why:** Structured extraction across all papers makes meta-analysis trivial. Export to R/Python for further stats.

### Citation Mining

**Goal:** Get DOIs and PMIDs for a reference manager.

**Script:** `test_fetcher_integrator.sh` ✅

**Why:** Raw metadata is all you need, faster.

---

## 📊 Output Format Comparison

### test_data_analyzer.sh Output

**File:** `Data_Analyzer/output/classified_papers_<timestamp>.csv`

```csv
PMID,Title,Relevance_Score,Is_Relevant,Organisms,Tissues,Conditions,Strategies,Year,Journal,DOI,Abstract_Preview
41068586,"Transcriptome...",10,Yes,"Solanum lycopersicum","root","drought stress ; water stress","RNA-Seq ; qRT-PCR",2025,"BMC plant biology",10.1186/...,Plants possess...
```

**Columns:** 12 (including analysis)

### test_fetcher_integrator.sh Output

**File:** `pubmed_results_<timestamp>.csv`

```csv
PMID,Title,Authors,Year,Journal,Publication_Type,DOI,PMCID,URL,Abstract
41068586,"Transcriptome...","Murtaza M; Aqeel M; Waheeb SA",2025,"BMC plant biology","Journal Article; Meta-Analysis",10.1186/...,PMC12513158,https://...,Plants possess...
```

**Columns:** 10 (raw data only)

---

## 🚀 Recommendation

**Start with:** `test_data_analyzer.sh`

It includes everything `test_fetcher_integrator.sh` does, plus AI analysis. If Ollama is not available, it gracefully falls back to raw results.

**Only use test_fetcher_integrator.sh if:**
- You explicitly want raw results only
- You're running on a machine where Ollama cannot be installed
- You need maximum speed and don't care about analysis

---

**Questions?** See [README.md](README.md) for installation and [Data_Analyzer/README.md](Data_Analyzer/README.md) for analysis details.
