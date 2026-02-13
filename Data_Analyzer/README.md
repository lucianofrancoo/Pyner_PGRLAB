# Data Analyzer

**Automated paper classification and analysis using Ollama LLM**

Analyzes papers from Fetcher_NCBI output (PubMed results) and extracts structured information using AI-powered classification.

---

## 🎯 Purpose

After fetching papers with `Fetcher_NCBI`, this module:
1. **Evaluates relevance** of each paper to your query (0-10 scale)
2. **Extracts structured metadata**:
   - Organisms mentioned (e.g., "Solanum lycopersicum", "Arabidopsis thaliana")
   - Tissues/organs (e.g., "root", "leaf", "shoot")
   - Experimental conditions (e.g., "drought stress", "heat shock")
   - Experimental strategies (e.g., "RNA-Seq", "qRT-PCR")
3. **Generates classified table** in CSV format

---

## 📦 Requirements

### 1. Ollama + Qwen Model

Install Ollama and pull the model:

```bash
# Install Ollama (if not already installed)
curl -fsSL https://ollama.com/install.sh | sh

# Pull Qwen model
ollama pull qwen2.5:14b

# Start Ollama server
ollama serve
```

### 2. BioPython (for PMC full text fetching)

```bash
pip3 install biopython
```

**Note:** BioPython is required for fetching full text from PubMed Central when PMCID is available.

### 3. Python Dependencies

No additional dependencies required beyond standard library (requests is used, should be available).

---

## 🔬 PMC Full Text Enhancement

**NEW FEATURE:** When a paper has a PMCID (PubMed Central ID), the analyzer automatically fetches the **full text** from PubMed Central instead of using only the abstract.

### Why Full Text?

- **Abstracts** typically contain ~200-500 characters of summary information
- **Full text** (Methods + Results sections) contains ~3,000-6,000 characters with detailed experimental information
- **Better extraction** of organisms (strains, varieties), tissues, conditions, and techniques

### How it Works

1. Paper has PMCID → Fetch full text from PMC
2. Extract Methods + Results sections
3. Send to Ollama for classification (~15x more context than abstract)
4. Fallback to abstract if PMC fetch fails

### Example Improvement

**Paper: PMID 36904041**

| Source | Organisms Extracted | Tissues Extracted |
|--------|-------------------|------------------|
| Abstract only | Solanum lycopersicum | root ; shoot |
| PMC full text | Solanum lycopersicum ; Micro-Tom | root ; shoot ; leaf |

### Statistics

After analysis, you'll see:
```
PMC full text used: 3
Abstract only:       0
```

### Configuration

In `config.py`:
```python
USE_PMC_FULL_TEXT = True  # Enable/disable PMC fetching
MAX_FULL_TEXT_LENGTH = 15000  # Max characters to send to LLM
```

---

## 🚀 Quick Start

### Basic Usage

```bash
cd Data_Analyzer
python3 paper_analyzer.py <path_to_fetcher_json>
```

**Example:**
```bash
python3 paper_analyzer.py ../pubmed_results_20260213_113559.json
```

This will:
- Read the JSON output from Fetcher
- Analyze each paper with Ollama
- Generate `classified_papers_<timestamp>.csv` in `output/`

### Specify Output File

```bash
python3 paper_analyzer.py ../pubmed_results.json my_classified_papers.csv
```

### Test Ollama Connection

```bash
python3 paper_analyzer.py --test-connection
```

---

## 📊 Output Format

### CSV Columns

| Column | Description | Example |
|--------|-------------|---------|
| **PMID** | PubMed ID | 41068586 |
| **Title** | Paper title | "Transcriptome-based meta-analysis..." |
| **Relevance_Score** | Relevance to query (0-10) | 9 |
| **Is_Relevant** | Relevant? (Yes/No) | Yes |
| **Organisms** | Organisms mentioned | "Solanum lycopersicum ; tomato" |
| **Tissues** | Tissues/organs | "root ; leaf" |
| **Conditions** | Experimental conditions | "drought stress ; water deficit" |
| **Strategies** | Experimental techniques | "RNA-Seq ; qRT-PCR" |
| **Year** | Publication year | 2025 |
| **Journal** | Journal name | "BMC plant biology" |
| **DOI** | DOI identifier | 10.1186/s12870-025-07348-2 |
| **Abstract_Preview** | First 200 chars of abstract | "Plants possess various..." |

**Note:** Multiple values in a field are separated by ` ; ` (space-semicolon-space)

### Example Output

```csv
PMID,Title,Relevance_Score,Is_Relevant,Organisms,Tissues,Conditions,Strategies,Year,Journal,DOI
41068586,"Transcriptome-based meta-analysis...",9,Yes,"Solanum lycopersicum ; tomato","root ; leaf","drought stress ; water deficit","RNA-Seq ; meta-analysis",2025,"BMC plant biology",10.1186/s12870-025-07348-2
```

---

## ⚙️ Configuration

Edit `config.py` to customize:

```python
# Ollama settings
OLLAMA_BASE_URL = "http://localhost:11434"
OLLAMA_MODEL = "qwen2.5:14b"
OLLAMA_TIMEOUT = 120

# Relevance threshold (papers with score >= 5 marked as relevant)
RELEVANCE_THRESHOLD = 5

# Maximum abstract length for analysis (chars)
MAX_ABSTRACT_LENGTH = 3000

# Output formatting
MULTIVALUE_SEPARATOR = " ; "  # Separator for multiple values
```

---

## 🔍 How It Works

### Analysis Pipeline

```
Input JSON → For each paper:
                1. Extract title + abstract
                2. Send to Ollama with structured prompt
                3. Parse JSON response
                4. Validate and format data
                5. Add to results table
           → Save CSV
```

### Ollama Prompt Strategy

The system uses a carefully crafted prompt that:
- Provides the **user's original query** for relevance context
- Asks for **structured JSON output** with specific fields
- Uses **low temperature** (0.1) for consistent extraction
- Limits **response length** to avoid hallucination

### Relevance Scoring

- **0-2**: Completely irrelevant
- **3-4**: Tangentially related
- **5-6**: Somewhat relevant
- **7-8**: Relevant
- **9-10**: Highly relevant

Papers with score ≥ 5 are marked as "Yes" in `Is_Relevant` column.

---

## 📈 Performance

**Speed:** ~5-10 seconds per paper (with Qwen 14B model)

For 100 papers: ~10-15 minutes

**Accuracy:** Depends on:
- Abstract quality and completeness
- Specificity of terms used in paper
- Model used (14B is balance of speed/accuracy)

---

## 🔧 Troubleshooting

### "ERROR: Ollama is not available"

**Solution:**
```bash
# Check if Ollama is running
curl http://localhost:11434/api/tags

# If not, start it
ollama serve
```

### "Model not found"

**Solution:**
```bash
ollama pull qwen2.5:14b
```

### Slow analysis

**Options:**
1. Use a smaller model: `ollama pull qwen2.5:7b` (edit `config.py`)
2. Reduce `MAX_ABSTRACT_LENGTH` in `config.py`
3. Use GPU if available (Ollama will use it automatically)

### Empty/error results

Check the log file in `logs/analyzer_<timestamp>.log` for detailed error messages.

---

## 📁 Directory Structure

```
Data_Analyzer/
├── config.py              # Configuration
├── ollama_client.py      # Ollama API wrapper
├── paper_analyzer.py     # Main script
├── README.md             # This file
├── output/               # Generated CSV files
│   └── classified_papers_<timestamp>.csv
└── logs/                 # Analysis logs
    └── analyzer_<timestamp>.log
```

---

## 🔗 Integration with Full Workflow

**Complete Pipeline:**

```bash
# 1. Fetch papers
cd ../
bash test_fetcher_integrator.sh
# → Generates: pubmed_results_<timestamp>.json

# 2. Analyze papers
cd Data_Analyzer
python3 paper_analyzer.py ../pubmed_results_<timestamp>.json
# → Generates: output/classified_papers_<timestamp>.csv

# 3. Review results
cat output/classified_papers_<timestamp>.csv
```

---

## 💡 Tips

### Multiple Organisms/Tissues/Conditions

The analyzer can extract multiple values per category:

**Example:**
- **Organisms**: `"Arabidopsis thaliana ; Col-0 ecotype"`
- **Conditions**: `"drought stress ; 21 days ; PEG treatment"`

Simply split by ` ; ` to get individual values.

### Filtering Results

```bash
# Get only relevant papers (score >= 5)
awk -F',' 'NR==1 || $4=="Yes"' classified_papers.csv > relevant_only.csv

# Get papers about specific organism
grep "Arabidopsis" classified_papers.csv > arabidopsis_papers.csv

# Count papers by relevance
awk -F',' 'NR>1 {print $4}' classified_papers.csv | sort | uniq -c
```

---

## 🚀 Future Enhancements

- [ ] Parallel processing for faster analysis
- [ ] Support for BioProject results (not just PubMed)
- [ ] Custom extraction schemas
- [ ] Export to JSON/Excel formats
- [ ] Interactive refinement of classifications
- [ ] Integration with vector databases for similarity search

---

## 📝 License

Part of the Pyner project. See main repository for license information.

---

**Questions?** Check the main documentation in `docs/FETCHER_DOCUMENTATION.md`
