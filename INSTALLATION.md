# Pyner Installation Guide

Quick guide to set up the environment for running Pyner workflows.

## Prerequisites

### 1. Python 3.8+

**Check if installed:**
```bash
python3 --version
```

**Install on Ubuntu/Debian:**
```bash
sudo apt update
sudo apt install python3 python3-pip
```

**Install on macOS:**
```bash
brew install python3
```

### 2. Python Dependencies

**Install BioPython (required for NCBI data fetching):**
```bash
pip install biopython
```

**Or install all dependencies at once:**
```bash
cd Pyner_PGRLAB
pip install -r Fetcher_NCBI/requirements.txt
```

### 3. Optional: FastAPI/Pydantic (only for API server mode)

The CLI workflow works without FastAPI. Only install if you plan to run the API server:

```bash
pip install fastapi uvicorn pydantic
```

## Quick Start

After installing dependencies, run the integrated workflow:

```bash
bash test_fetcher_integrator.sh
```

The script will automatically:
- ✅ Validate that Python 3 is available
- ✅ Check that BioPython is installed
- ✅ Verify all required files exist
- ❌ Show clear error messages if something is missing

## Configuration

### NCBI API Key (Recommended)

To avoid rate limits, configure your NCBI credentials:

1. Get an API key from: https://www.ncbi.nlm.nih.gov/account/settings/

2. Edit `Fetcher_NCBI/config.py`:
```python
NCBI_EMAIL = "your.email@example.com"
NCBI_API_KEY = "your_api_key_here"
```

Or set environment variables:
```bash
export NCBI_EMAIL="your.email@example.com"
export NCBI_API_KEY="your_api_key_here"
```

## Troubleshooting

### python3: command not found
**Solution:** Install Python 3 (see Prerequisites above)

### ModuleNotFoundError: No module named 'Bio'
**Solution:** Install BioPython:
```bash
pip install biopython
```

### ModuleNotFoundError: No module named 'fastapi'
**This is OK for CLI mode.** FastAPI is only needed for server mode.
The integrated workflow (`test_fetcher_integrator.sh`) works without it.

### Permission denied: test_fetcher_integrator.sh
**Solution:** Make the script executable:
```bash
chmod +x test_fetcher_integrator.sh
```

## What Gets Validated

When you run `test_fetcher_integrator.sh`, it checks:

1. ✅ **Python 3** - Must be available in PATH
2. ✅ **BioPython** - Required for NCBI API access
3. ✅ **Required files** - All Python scripts must exist
4. ❌ Shows clear error messages with installation instructions

## Verification

Test that everything works:

```bash
# Test BioPython installation
python3 -c "import Bio; print('BioPython OK')"

# Test the full workflow
bash test_fetcher_integrator.sh
```

## Complete Installation Script

```bash
# Clone repository
cd ~/
git clone https://github.com/lucianofrancoo/Pyner_PGRLAB.git
cd Pyner_PGRLAB

# Install dependencies
pip install biopython

# Configure NCBI (optional but recommended)
export NCBI_EMAIL="your.email@example.com"
export NCBI_API_KEY="your_api_key_here"

# Run workflow
bash test_fetcher_integrator.sh
```

## Minimal Dependencies

The absolute minimum to run the workflow:

```bash
pip install biopython
```

That's it! Everything else is optional.

## For Developers

If you want to run the API server or contribute code:

```bash
# Install all dev dependencies
pip install biopython fastapi uvicorn pydantic

# Start API server (optional)
cd Query_generator/phases/phase3
python3 api/main.py --server
```

---

**Questions?** Check the documentation in `docs/FETCHER_DOCUMENTATION.md`
