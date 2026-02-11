# Fetcher_NCBI

**Fetch biological sequencing data from NCBI SRA using boolean queries**

Fetcher_NCBI queries the NCBI Sequence Read Archive (SRA) database and retrieves metadata from sequencing studies. It integrates seamlessly with the Query Generator to transform natural language into NCBI queries and fetch relevant datasets.

---

## Features

✅ **Query NCBI SRA** with boolean search syntax  
✅ **Automatic deduplication** by BioProject ID  
✅ **Rate limiting** respects NCBI API guidelines (3 req/sec without key, 10 req/sec with key)  
✅ **Metadata extraction** from XML responses (organism, library strategy, tissue, etc.)  
✅ **JSON output** for easy integration with downstream analysis  
✅ **CLI and Python API** for flexible usage

---

## Quick Start

### 1. Configure Credentials

Edit [config.py](config.py) or set environment variables:

```bash
export NCBI_EMAIL="your.email@example.com"
export NCBI_API_KEY="your_api_key_here"  # Optional but recommended
```

Get your API key at: https://www.ncbi.nlm.nih.gov/account/settings/

### 2. Run a Search

```bash
# Direct query
python main.py --query "arabidopsis[Organism] AND drought"

# From file
python main.py --query-file query.txt

# Specify output and max results
python main.py -q "stress response" -m 500 -o results.json
```

### 3. Integration with Query Generator

Use the Query Generator to create optimized queries, then fetch results:

```bash
# Step 1: Generate query
cd ../Query_generator/phases/phase3
python generate_query.py "arabidopsis drought stress RNA-seq" > query.json

# Step 2: Fetch data
cd ../../../Fetcher_NCBI
python main.py --query-file ../Query_generator/phases/phase3/query.json
```

---

## Architecture

```
Fetcher_NCBI/
├── config.py              # Configuration and credentials
├── ncbi_fetcher.py        # Core fetching logic
├── main.py                # CLI interface
├── data/                  # Output data storage
│   ├── sra_results.json   # Default output file
│   └── bioproject_cache.json  # Deduplication cache
└── logs/                  # Execution logs
```

### Key Components

#### `config.py`
- Manages NCBI credentials (email, API key)
- Configurable rate limits and output paths
- Validation checks for missing credentials

#### `ncbi_fetcher.py`
- **`NCBIFetcher`**: Main class for querying NCBI
- **`BioProjectCache`**: Tracks processed BioProjects to avoid duplicates
- **`extract_metadata()`**: Parses XML responses and extracts structured data
- Uses BioPython's `Entrez` API for NCBI communication

#### `main.py`
- Command-line interface with argparse
- Supports direct queries or loading from files
- JSON and plain text query formats

---

## Usage Examples

### Basic Queries

```bash
# Simple organism search
python main.py -q "arabidopsis[Organism]"

# Organism + condition
python main.py -q "arabidopsis[Organism] AND drought"

# Multiple conditions with Boolean operators
python main.py -q "arabidopsis[Organism] AND (drought OR 'water stress')"

# Specific sequencing strategy
python main.py -q "arabidopsis[Organism] AND RNA-Seq[Strategy]"
```

### Advanced Options

```bash
# Limit results
python main.py -q "drought stress" --max-results 100

# Custom output location
python main.py -q "arabidopsis" -o /path/to/output.json

# Disable deduplication (fetch all results)
python main.py -q "stress" --no-deduplicate

# Clear cache and start fresh
python main.py -q "drought" --clear-cache
```

### Python API

```python
from ncbi_fetcher import NCBIFetcher

# Initialize fetcher
fetcher = NCBIFetcher(
    email="your.email@example.com",
    api_key="your_api_key"
)

# Fetch data
results = fetcher.fetch_all(
    query="arabidopsis[Organism] AND drought",
    max_results=1000,
    deduplicate=True
)

# Save results
fetcher.save_results("output.json")

# Access results
for study in results:
    print(f"BioProject: {study['bioproject']}")
    print(f"Organism: {study['organism']}")
    print(f"Strategy: {study['library_strategy']}")
```

---

## Output Format

Results are saved as JSON with the following structure:

```json
{
  "metadata": {
    "fetched_at": "2024-01-15T10:30:00",
    "total_results": 42,
    "statistics": {
      "total_found": 150,
      "processed": 42,
      "skipped_duplicate": 8,
      "skipped_error": 0,
      "unique_bioprojects": 42
    }
  },
  "results": [
    {
      "sra_id": "SRX123456",
      "bioproject": "PRJNA123456",
      "title": "Arabidopsis drought stress RNA-seq",
      "study_name": "GSE123456",
      "organism": "Arabidopsis thaliana",
      "library_strategy": "RNA-Seq",
      "library_source": "TRANSCRIPTOMIC",
      "library_selection": "cDNA",
      "library_layout": "PAIRED",
      "library_name": "leaf tissue",
      "biosample": "SAMN123456",
      "fetched_at": "2024-01-15T10:30:05"
    }
  ]
}
```

### Key Fields

- **`sra_id`**: NCBI SRA accession (e.g., SRX123456)
- **`bioproject`**: Unique study identifier (used for deduplication)
- **`organism`**: Scientific name (e.g., Arabidopsis thaliana)
- **`library_strategy`**: Sequencing type (RNA-Seq, WGS, ChIP-Seq, etc.)
- **`library_layout`**: PAIRED or SINGLE-end
- **`biosample`**: Sample identifier

---

## Deduplication

By default, Fetcher_NCBI tracks processed BioProjects to avoid fetching the same study multiple times.

### How It Works

1. Each BioProject ID is stored in `data/bioproject_cache.json`
2. When a duplicate is encountered, it's skipped
3. Cache persists across runs to avoid re-fetching old data

### Managing Cache

```bash
# Clear cache (reset deduplication)
python main.py -q "drought" --clear-cache

# Disable deduplication for one run
python main.py -q "drought" --no-deduplicate
```

---

## Rate Limiting

NCBI enforces strict rate limits:

| Configuration | Rate Limit | Requests/Second |
|--------------|------------|-----------------|
| **No API Key** | 3 requests/sec | ⚠️ Slow |
| **With API Key** | 10 requests/sec | ✅ Fast |

**Recommendation**: Always use an API key for faster fetching.

Get your key: https://www.ncbi.nlm.nih.gov/account/settings/

---

## Troubleshooting

### Configuration Warnings

If you see warnings about missing credentials:

```
⚠️  NCBI_EMAIL not configured
⚠️  NCBI_API_KEY not set
```

**Solution**: Edit [config.py](config.py) or set environment variables:

```bash
export NCBI_EMAIL="your.email@example.com"
export NCBI_API_KEY="your_api_key"
```

### No Results Found

If fetcher returns no results:

1. **Check query syntax**: Use NCBI's search builder to test: https://www.ncbi.nlm.nih.gov/sra/
2. **Try broader terms**: Remove specific filters
3. **Clear cache**: Old cache might skip all results
   ```bash
   python main.py -q "your query" --clear-cache
   ```

### Rate Limit Errors

If you hit rate limits:

1. **Add API key**: Increases limit from 3/sec to 10/sec
2. **Reduce max results**: Fetch fewer results per run
3. **Wait and retry**: NCBI blocks are temporary

---

## Integration Workflow

### Complete Pipeline: Query Generator → Fetcher_NCBI

```bash
# 1. Generate optimized query from natural language
cd Query_generator/phases/phase3
python api/main.py cli "arabidopsis drought stress rna-seq"

# Output: arabidopsis[Organism] AND ("drought"[All Fields] OR "water stress"[All Fields]) ...

# 2. Fetch data using generated query
cd ../../../Fetcher_NCBI
python main.py -q 'arabidopsis[Organism] AND ("drought"[All Fields] OR "water stress")' -o drought_data.json

# 3. Results ready for analysis
cat data/drought_data.json
```

---

## Development

### Dependencies

Install required packages:

```bash
pip install biopython
```

### Testing

```bash
# Test configuration
python config.py

# Test fetcher with small query
python ncbi_fetcher.py

# Test CLI
python main.py -q "test[Organism]" -m 5
```

### Logging

All operations are logged to `logs/fetcher_YYYYMMDD_HHMMSS.log`

Check logs for:
- API errors
- XML parsing issues
- Rate limiting warnings
- Duplicate skips

---

## Roadmap

- [ ] Support for additional NCBI databases (GEO, Nucleotide, Protein)
- [ ] Parallel fetching for faster retrieval
- [ ] CSV export in addition to JSON
- [ ] Filter results by date, publication status, etc.
- [ ] Integration with LLM analysis (port from Pyner_v0.2.py)

---

## Reference

Based on the original [Pyner_v0.2.py](../scripts_iniciales_beta/Pyner_v0.2.py) script.

**API Documentation**: https://www.ncbi.nlm.nih.gov/books/NBK25501/  
**SRA Search Guide**: https://www.ncbi.nlm.nih.gov/sra/docs/srasearch/
