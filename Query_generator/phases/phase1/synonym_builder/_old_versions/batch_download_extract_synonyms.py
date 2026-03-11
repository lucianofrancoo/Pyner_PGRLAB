#!/usr/bin/env python3
"""
Batch download PubMed Baseline XMLs and extract MeSH synonyms.

For each XML file:
1. Download from NCBI FTP
2. Verify MD5
3. Extract MeSH terms with their synonyms (from Entry Terms in XML)
4. Merge synonyms (add new ones if term already exists)
5. Delete XML to save space
6. Continue with next file

Output: Dictionary of MeSH terms with their synonyms
"""

import os
import sys
import json
import gzip
import hashlib
import logging
import requests
import time
from pathlib import Path
from datetime import datetime
from lxml import etree
from collections import defaultdict

import config

# Setup logging
log_file = config.OUTPUT_DIR / f"batch_synonyms_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

# Configuration
NCBI_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"
DOWNLOAD_DIR = config.PHASE1_DIR / "pubmed_baseline"  # Use pubmed_baseline folder
DOWNLOAD_DIR.mkdir(exist_ok=True)

# Checkpoint and data files
CHECKPOINT_FILE = config.OUTPUT_DIR / "synonym_extraction_checkpoint.json"
SYNONYMS_OUTPUT = config.OUTPUT_DIR / "mesh_synonyms_from_baseline.json"


def load_checkpoint():
    """Load checkpoint to resume processing."""
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE, 'r') as f:
            return json.load(f)
    return {
        "last_processed": None,
        "processed_count": 0,
        "failed": [],
        "started_at": datetime.now().isoformat()
    }


def save_checkpoint(checkpoint):
    """Save checkpoint."""
    checkpoint["last_updated"] = datetime.now().isoformat()
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(checkpoint, f, indent=2)


def load_synonyms():
    """Load existing synonym dictionary."""
    if SYNONYMS_OUTPUT.exists():
        logger.info(f"Loading existing synonyms...")
        with open(SYNONYMS_OUTPUT, 'r') as f:
            data = json.load(f)
        logger.info(f"  Loaded {len(data):,} unique MeSH terms")
        return data
    return {}


def save_synonyms(synonyms_dict):
    """Save synonym dictionary."""
    logger.info(f"Saving synonyms...")
    with open(SYNONYMS_OUTPUT, 'w') as f:
        json.dump(synonyms_dict, f, indent=2, ensure_ascii=False)

    file_size_mb = SYNONYMS_OUTPUT.stat().st_size / (1024 * 1024)
    logger.info(f"  Saved {len(synonyms_dict):,} unique MeSH terms ({file_size_mb:.2f} MB)")


def calculate_md5(file_path):
    """Calculate MD5 checksum."""
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def download_file(url, dest_path, desc="file"):
    """Download file with progress."""
    try:
        logger.info(f"  Downloading {desc}...")
        response = requests.get(url, stream=True, timeout=300)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))
        downloaded = 0
        last_log = 0

        with open(dest_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)

                    if downloaded - last_log >= 10 * 1024 * 1024:
                        progress = (downloaded / total_size * 100) if total_size > 0 else 0
                        logger.info(f"    {downloaded/(1024*1024):.1f} MB / {total_size/(1024*1024):.1f} MB ({progress:.0f}%)")
                        last_log = downloaded

        logger.info(f"    ✓ Downloaded: {downloaded/(1024*1024):.1f} MB")
        return True
    except Exception as e:
        logger.error(f"    ✗ Download failed: {e}")
        return False


def extract_synonyms_from_xml(xml_path, synonyms_dict):
    """
    Extract MeSH terms and their synonyms from XML.

    Looks for Entry Terms which are the synonyms in PubMed XML.
    If term already exists, merges synonyms (adds new ones, no duplicates).
    """
    logger.info(f"  Extracting synonyms from {xml_path.name}...")

    stats = {
        "total_articles": 0,
        "articles_with_mesh": 0,
        "new_terms": 0,
        "updated_terms": 0,
        "new_synonyms_added": 0
    }

    try:
        context = etree.iterparse(str(xml_path), events=('end',), tag='PubmedArticle')

        for event, elem in context:
            stats["total_articles"] += 1

            # Check if has MeSH
            mesh_list = elem.find('.//MeshHeadingList')
            if mesh_list is not None:
                stats["articles_with_mesh"] += 1

                # Extract each MeSH heading
                for mesh_heading in mesh_list.findall('MeshHeading'):
                    descriptor = mesh_heading.find('DescriptorName')
                    if descriptor is None or descriptor.text is None:
                        continue

                    mesh_term = descriptor.text
                    key = mesh_term.lower()

                    # Extract synonyms from XML (Entry Terms are in the broader context)
                    # Note: Entry Terms are not in individual articles, they're in the MeSH ontology
                    # For now, we'll just collect unique terms and count occurrences
                    # Synonyms will come from analyzing text variations in Step 2B

                    # Create or update entry
                    if key not in synonyms_dict:
                        synonyms_dict[key] = {
                            "preferred_term": mesh_term,
                            "synonyms": [],  # Will be filled from abstract analysis later
                            "article_count": 0
                        }
                        stats["new_terms"] += 1
                    else:
                        stats["updated_terms"] += 1

                    # Increment count
                    synonyms_dict[key]["article_count"] += 1

            # Progress
            if stats["total_articles"] % 10000 == 0:
                logger.info(f"    Progress: {stats['total_articles']:,} articles, "
                          f"{len(synonyms_dict):,} unique terms")

            # Clean memory
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]

        logger.info(f"    ✓ Completed: {stats['total_articles']:,} articles")
        logger.info(f"      New terms: {stats['new_terms']:,}")
        logger.info(f"      Updated terms: {stats['updated_terms']:,}")
        logger.info(f"      Total unique terms: {len(synonyms_dict):,}")

        return stats

    except Exception as e:
        logger.error(f"    ✗ Error processing: {e}")
        return stats


def process_file(file_num, synonyms_dict):
    """Process a single baseline file."""
    filename = f"pubmed26n{file_num:04d}.xml.gz"

    logger.info(f"\n{'='*80}")
    logger.info(f"FILE {file_num}/1334: {filename}")
    logger.info(f"{'='*80}")

    # Paths
    gz_path = DOWNLOAD_DIR / filename
    md5_path = DOWNLOAD_DIR / f"{filename}.md5"
    xml_path = DOWNLOAD_DIR / filename.replace('.gz', '')

    success = False
    stats = None

    try:
        # 1. Download .gz
        if not download_file(NCBI_FTP_BASE + filename, gz_path, ".gz file"):
            return False, None

        # 2. Download MD5
        if not download_file(NCBI_FTP_BASE + f"{filename}.md5", md5_path, "MD5 file"):
            return False, None

        # 3. Verify MD5
        logger.info(f"  Verifying MD5...")
        with open(md5_path, 'r') as f:
            content = f.read().strip()
            if '=' in content:
                expected_md5 = content.split('=')[-1].strip()
            else:
                expected_md5 = content.split()[0]

        actual_md5 = calculate_md5(gz_path)

        if actual_md5 != expected_md5:
            logger.error(f"    ✗ MD5 mismatch! Expected: {expected_md5}, Got: {actual_md5}")
            return False, None

        logger.info(f"    ✓ MD5 verified")

        # 4. Decompress
        logger.info(f"  Decompressing...")
        with gzip.open(gz_path, 'rb') as f_in:
            with open(xml_path, 'wb') as f_out:
                f_out.write(f_in.read())
        logger.info(f"    ✓ Decompressed: {xml_path.stat().st_size/(1024*1024):.1f} MB")

        # 5. Extract synonyms
        stats = extract_synonyms_from_xml(xml_path, synonyms_dict)

        success = True
        logger.info(f"✓ Successfully processed {filename}")

    except Exception as e:
        logger.error(f"✗ Error processing {filename}: {e}")
        success = False

    finally:
        # Cleanup: Delete downloaded files
        logger.info(f"  Cleaning up...")
        for temp_file in [gz_path, md5_path, xml_path]:
            if temp_file.exists():
                temp_file.unlink()
        logger.info(f"    ✓ Files deleted from {DOWNLOAD_DIR.name}/")

    return success, stats


def main():
    """Main execution."""
    logger.info("=" * 80)
    logger.info("PUBMED BASELINE - MeSH SYNONYM EXTRACTION")
    logger.info("=" * 80)
    logger.info(f"Started: {datetime.now()}")
    logger.info(f"Log file: {log_file}")
    logger.info(f"Download to: {DOWNLOAD_DIR}")
    logger.info(f"Output: {SYNONYMS_OUTPUT}")
    logger.info(f"Checkpoint: {CHECKPOINT_FILE}")
    logger.info("")

    # Load checkpoint
    checkpoint = load_checkpoint()
    start_from = checkpoint["processed_count"] + 1

    # Load existing synonyms
    synonyms_dict = load_synonyms()

    logger.info(f"Resume from: File {start_from}/1334")
    logger.info(f"Already processed: {checkpoint['processed_count']} files")
    logger.info(f"Unique MeSH terms: {len(synonyms_dict):,}")
    logger.info("")

    start_time = time.time()

    # Process files 1 to 1334
    for file_num in range(start_from, 1335):
        success, stats = process_file(file_num, synonyms_dict)

        # Update checkpoint
        if success:
            checkpoint["last_processed"] = f"pubmed26n{file_num:04d}.xml.gz"
            checkpoint["processed_count"] = file_num
        else:
            checkpoint["failed"].append(f"pubmed26n{file_num:04d}.xml.gz")

        # Save checkpoint and synonyms after each file
        save_checkpoint(checkpoint)
        save_synonyms(synonyms_dict)

        # Calculate ETA
        elapsed = time.time() - start_time
        processed_this_session = file_num - start_from + 1
        avg_time = elapsed / processed_this_session if processed_this_session > 0 else 0
        remaining = 1334 - file_num
        eta_seconds = avg_time * remaining
        eta_hours = eta_seconds / 3600

        logger.info(f"\nProgress: {file_num}/1334 ({file_num/1334*100:.1f}%)")
        logger.info(f"Unique MeSH terms: {len(synonyms_dict):,}")
        logger.info(f"Avg time/file: {avg_time:.1f}s")
        logger.info(f"ETA: {eta_hours:.1f} hours")

        # Small delay between files
        if file_num < 1334:
            time.sleep(2)

    # Final summary
    total_time = time.time() - start_time
    logger.info("\n" + "=" * 80)
    logger.info("EXTRACTION COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Total time: {total_time/3600:.2f} hours")
    logger.info(f"Files processed: {checkpoint['processed_count']}")
    logger.info(f"Files failed: {len(checkpoint['failed'])}")
    logger.info(f"Unique MeSH terms: {len(synonyms_dict):,}")
    logger.info(f"Output file: {SYNONYMS_OUTPUT}")

    # Show top 10 most frequent terms
    logger.info("\nTop 10 most frequent MeSH terms:")
    sorted_terms = sorted(
        synonyms_dict.items(),
        key=lambda x: x[1]["article_count"],
        reverse=True
    )
    for i, (key, data) in enumerate(sorted_terms[:10], 1):
        logger.info(f"  {i:2d}. {data['preferred_term']:30s} - {data['article_count']:,} articles")

    logger.info("")
    logger.info("Note: Synonyms list is empty for now.")
    logger.info("Next step: Analyze abstracts to extract text variations (Step 2B)")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("\n\n⚠️  Interrupted by user")
        logger.info("✓ Progress saved in checkpoint")
        logger.info("✓ Run again to resume")
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
