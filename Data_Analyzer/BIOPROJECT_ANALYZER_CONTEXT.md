# PYNER MINER - BioProject Analyzer Context & Requirements

## 1. Project Context
**Project Name:** PYNER MINER
**Core Purpose:** An automated intelligence pipeline designed for biological and agricultural research. It transforms natural language queries into complex Boolean statements, queries NCBI databases, and uses Local Large Language Models (LLMs via Ollama) to extract highly structured metadata from scientific literature.

**Current Architecture:**
1. **Phase 1 (Query Expander):** Uses a "Zero-Bias" semantic approach. It combines a Canonical Biological Lexicon (MeSH + Taxonomy + Gene Dictionary alias shielding) and an Empirical Condition Lexicon (ECL) extracted from 23 million abstracts to build robust PubMed/NCBI Boolean queries.
2. **Phase 2 (Data Fetcher):** Executes the queries against either the `pubmed` or `bioproject` databases.
    *   *PubMed Mode:* Retrieves classic paper abstracts.
    *   *BioProject Mode:* Executes a "Cascade Search" (BioProject → SRA Experiments → BioSamples → linked PubMed IDs).
3. **Phase 3 (Data Analyzer):** Currently only works for `pubmed` mode. It feeds paper abstracts (and PMC full-texts if available) to a local LLM (`qwen2.5`) via `paper_analyzer.py` and `ollama_client.py`, returning a TSV with 53 structured metadata columns.

## 2. The Problem
Currently, when a user selects "BioProject" as the search target (`pyner_miner.sh -d bioproject`), Phase 2 successfully fetches the data and creates a JSON with a completely different structure than the PubMed JSON.
However, **Phase 3 (Data Analyzer) is hardcoded to expect PubMed structures** (it looks for titles and abstracts). When the BioProject JSON arrives, the analyzer fails or does nothing because it doesn't understand the schema.

## 3. The Objective: Creating the "BioProject Analyzer"
We want to create a **brand-new, independent analyzer** specifically for BioProjects to avoid polluting or breaking the highly functional `paper_analyzer.py` script.

**Input Data Structure:**
The input JSON from the BioProject Fetcher (`boolean_fetcher_integrated.py`) has this structure:
```json
{
  "metadata": {...},
  "results": [
    {
      "bioproject": "PRJNA12345",
      "title": "Transcriptome of Solanum lycopersicum under drought",
      "submission_date": "2023-05-12",
      "organism": "Solanum lycopersicum",
      "project_type": "Transcriptome or Gene expression",
      "description": "RNA-seq data from leaves exposed to 14 days of water deficit...",
      "sra_experiments_count": 12,
      "sra_hierarchy": {
        "SAMN12345": {
          "sample_id": "SAMN12345",
          "experiments": [
            {
              "experiment_id": "SRX123",
              "metadata": {
                "library_strategy": "RNA-Seq",
                "instrument": "Illumina HiSeq"
             },
              "runs": [{"accession": "SRR123", "spots": 25000000, "size": "2.5 GB"}]
            }
          ]
        }
      },
      "pmids": ["31234567"]
    }
  ]
}
```

## 4. Key Requirements for the New Agent
The AI agent tasked with building this script should follow these guidelines:

1.  **Architecture:** Create two new scripts: `bioproject_analyzer.py` (Orchestrator, similar to `paper_analyzer.py`) and a modified LLM client method (or separate class) to handle the new prompt.
2.  **LLM Task Focus:** The LLM's goal is no longer evaluating a "Paper." The LLM must evaluate the **Data Usability**.
    *   *Prompt logic:* "Based on the description of this BioProject and its SRA metadata (such as RNA-Seq, Illumina, etc.), how useful are these raw datasets for studying [USER QUERY]?"
3.  **Output Structure (New TSV schema):** We need to define new columns relevant to data-mining, such as:
    *   `BioProject_ID`, `Title`, `Organism`
    *   `Data_Type` (e.g., RNA-Seq, ChIP-Seq, WGS)
    *   `Total_Samples` (extracted from the SRA count)
    *   `Estimated_Volume` (Adding up the Gb from the runs, if possible)
    *   `Relevance_Score` (0-10, AI-generated)
    *   `Experimental_Design_Summary` (AI-generated based on description)
    *   `Has_Published_Paper` (Yes/No based on presence of PMIDs)
4.  **Integration:** Explain how to connect this new script inside the `pyner_miner.sh` bash wrapper, specifically within the `elif [ "$DATABASE" = "bioproject" ]; then` block.
