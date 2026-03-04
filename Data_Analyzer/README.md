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
ollama pull qwen3.5:9b

# Start Ollama server
ollama serve
```

### 2. BioPython (for PMC full text fetching)

```bash
pip3 install biopython
```

**Note:** BioPython is required for fetching full text from PubMed Central when PMCID is available.

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

---

## 📊 Output Format

### CSV Columns (Simplified view)

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

**Note:** The output includes **53 columns** of experimental metadata. See `EXPANDED_COLUMNS.md` for details.

---

## ⚙️ Configuration

Edit `config.py` to customize:

```python
# Ollama settings
OLLAMA_BASE_URL = "http://localhost:11434"
OLLAMA_MODEL = "qwen3.5:9b"
OLLAMA_TIMEOUT = 600
```

---

## 📈 Performance

**Speed:** ~5-15 seconds per paper (with Qwen 3.5 9B model)

**Accuracy:** Highly detailed extraction using the latest stable Qwen model.

---

## 🔧 Troubleshooting

### "Model not found"

**Solution:**
```bash
ollama pull qwen3.5:9b
```

---

## 📁 Directory Structure

```
Data_Analyzer/
├── config.py              # Configuration
├── ollama_client.py      # Ollama API wrapper
├── paper_analyzer.py     # Main script
├── README.md             # This file
├── output/               # Generated CSV files
└── logs/                 # Analysis logs
```

---

**Version:** 3.0.0-beta | **Date:** 2026-03-04
