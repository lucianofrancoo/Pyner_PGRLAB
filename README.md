# 🔬 PYNER — Intelligent Literature & Data Mining Pipeline (Beta)

An advanced system for searching and linking information from genomic projects (BioProject), experimental data (SRA), and scientific publications (PubMed/PMC) using natural language and AI-powered analysis.

---

## 🚀 Quick Start

### Main Execution
Simply run the central script:
```bash
bash pyner_miner.sh
```

**Integrated Workflow:**
1. **Mode:** Lite (fast) or Pro (deep analysis with 53 columns).
2. **Query:** Enter your question in natural language (e.g., "Tomato drought stress RNA-Seq").
3. **Database:** Choose between PubMed (literature) or BioProject (omics data).
4. **Boolean Generator:** AI automatically generates the NCBI boolean query.
5. **Mining & Analysis:** The system fetches data and, in Pro mode, classifies it using Ollama.
6. **Visualization:** Automatic generation of an interactive HTML paper network.

---

## 🏗️ Project Structure

The system is organized into 4 main modules:

- **`Query_generator`**: AI engine that translates natural language to NCBI boolean queries with synonym expansion.
- **`Fetcher_NCBI`**: Search orchestrator that links BioProjects, SRA, and PubMed.
- **`Data_Analyzer`**: Deep experimental metadata extractor (53 fields) using local LLMs (Ollama).
- **`Data_visualization`**: Interactive visual report generator (D3.js).

```
Pyner_PGRLAB/
├── pyner_miner.sh           # 🚀 Main entry point
├── Query_generator/         # 🤖 AI: Natural Language to Boolean
├── Fetcher_NCBI/            # 🔍 NCBI Metadata Fetching
├── Data_Analyzer/           # 📊 Deep Analysis & Classification
├── Data_visualization/      # 🔬 Paper Network Visualization (HTML)
├── output/                  # 📁 All results (TSV, JSON, HTML)
└── README.md                # 📖 This guide
```

---

## ✨ Key Features

### 📈 Interactive Visualization
At the end of each search, an `output/paper_network_TIMESTAMP.html` file is generated, allowing you to:
- View relationships between papers (by organism or journal).
- Filter by year, relevance, and journal in real-time.
- Toggle between **Light/Dark Mode**.
- **Drag and pin** nodes to organize your view (double-click to release).
- Auto-fit to screen to see all results.

### 📊 Pro Analysis (Ollama)
Automatic extraction of 53 metadata columns:
- **Biology:** Organisms, genotypes, tissues, cell types.
- **Experiment:** Conditions, molecules, time-course design, replication.
- **Results:** Relevance score (0-10) and detailed justification.

---

## 🔧 Requirements

- **Python 3.8+**
- **BioPython**: `pip install biopython`
- **Ollama** (for Pro/Query Gen modes): 
  1. Install Ollama: `curl -fsSL https://ollama.com/install.sh | sh`
  2. Pull a model (e.g., 4b, 9b, 2b, or 0.8b): `ollama pull qwen2.5:14b`
  3. Verify: `ollama list`


---
**Version:** 3.0.0-beta | **Date:** 2026-03-04 | **Status:** 🚀 Beta Ready
