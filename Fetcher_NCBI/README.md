# Fetcher_NCBI

**Fetch BioProjects from NCBI using boolean queries**

Fetcher_NCBI queries the NCBI BioProject database and retrieves project metadata. It integrates seamlessly with the Query Generator to transform natural language into NCBI queries and fetch relevant datasets.

---

## Features

✅ **Query NCBI BioProject** with boolean search syntax  
✅ **Automatic deduplication** by BioProject ID  
✅ **Rate limiting** respects NCBI API guidelines (3 req/sec without key, 10 req/sec with key)  
✅ **Metadata extraction** from BioProject summaries (title, organism, description)  
✅ **JSON/CSV output** for easy integration with downstream analysis  
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

# Specify output
python main.py -q "stress response" -o results.json
```

### 3. Integration with Query Generator

Use the Query Generator to create optimized queries, then fetch results:

```bash
# Step 1: Generate query
cd ../Query_generator/phases/phase3
python api/main.py --quick "arabidopsis drought stress RNA-seq"

# Step 2: Fetch data
cd ../../../Fetcher_NCBI
python main.py --query "PASTE_QUERY_HERE" -o results.json
```

---

## Architecture

```
Fetcher_NCBI/
├── config.py              # Configuration and credentials
├── ncbi_fetcher.py        # Core fetching logic
├── main.py                # CLI interface
├── data/                  # Output data storage
│   ├── bioproject_results.json   # Default output file
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
# Custom output location
python main.py -q "arabidopsis" -o /path/to/output.json

# Save as CSV table
python main.py -q "arabidopsis" -o /path/to/output.csv

# Disable deduplication (fetch all results)
python main.py -q "stress" --no-deduplicate

# Stop after collecting N unique BioProjects
python main.py -q "drought" --min-bioprojects 30

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
  deduplicate=True
)

# Save results
fetcher.save_results("output.json")

# Access results
for study in results:
    print(f"BioProject: {study['bioproject']}")
    print(f"Organism: {study['organism']}")
  print(f"Title: {study.get('title', '')}")
```

---

## Output Format

Results are saved as JSON or CSV. JSON structure:

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
      "sra_id": "12345678",
      "bioproject": "PRJNA123456",
      "title": "Arabidopsis drought stress RNA-seq",
      "description": "Study description...",
      "organism": "Arabidopsis thaliana",
      "fetched_at": "2024-01-15T10:30:05"
    }
  ]
}
```

### Key Fields

- **`sra_id`**: NCBI internal ID from BioProject search
- **`bioproject`**: Unique study identifier (used for deduplication)
- **`organism`**: Scientific name (e.g., Arabidopsis thaliana)
- **`title`**: Project title
- **`description`**: Project description

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

1. **Check query syntax**: Use NCBI's search builder to test: https://www.ncbi.nlm.nih.gov/bioproject/
2. **Try broader terms**: Remove specific filters
3. **Clear cache**: Old cache might skip all results
   ```bash
   python main.py -q "your query" --clear-cache
   ```

### Rate Limit Errors

If you hit rate limits:

1. **Add API key**: Increases limit from 3/sec to 10/sec
2. **Wait and retry**: NCBI blocks are temporary

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

# Or use the interactive integrator at repo root
cd ../..
bash test_fetcher_integrator.sh
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
python main.py -q "test[Organism]" -o data/test.csv
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
- [x] CSV export in addition to JSON
- [ ] Filter results by date, publication status, etc.
- [ ] Integration with LLM analysis (port from Pyner_v0.2.py)

---

## Reference

Based on the original [Pyner_v0.2.py](../scripts_iniciales_beta/Pyner_v0.2.py) script.

**API Documentation**: https://www.ncbi.nlm.nih.gov/books/NBK25501/  
**BioProject Search Guide**: https://www.ncbi.nlm.nih.gov/bioproject/
