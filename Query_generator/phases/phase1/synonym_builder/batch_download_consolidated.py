#!/usr/bin/env python3
"""
Batch download and process all 1334 PubMed Baseline XML files.
CONSOLIDATED VERSION: No duplicates, merges data as it processes.

For each MeSH term:
- Tracks how many articles mention it
- Extracts Title + Abstract + Keywords from each article
- Stores sample texts (limited to 100 per term to save memory)
- Automatically merges duplicates

Output format per term:
{
  "term_key": {
    "preferred_term": "Original MeSH Term",
    "mesh_synonyms": ["syn1", "syn2", ...],  # From Step 1
    "article_count": 12345,
    "sample_texts": [
      "Title text || Abstract text (500 chars) || keyword1 | keyword2 | keyword3",
      ...
    ]
  }
}

Updated: 2026-03-06 - Added Title and Keywords extraction
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
log_file = config.OUTPUT_DIR / f"batch_consolidated_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
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
TEMP_DIR = config.BASE_DIR / "temp_downloads"
TEMP_DIR.mkdir(exist_ok=True)

# Checkpoint and data files
CHECKPOINT_FILE = config.OUTPUT_DIR / "consolidated_checkpoint.json"
CONSOLIDATED_DATA_FILE = config.OUTPUT_DIR / "step2_mesh_consolidated.json"


def load_checkpoint():
    """Load checkpoint to resume processing."""
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE, 'r') as f:
            return json.load(f)
    return {
        "last_processed": None,
        "processed_count": 0,
        "failed": [],
        "total_articles": 0,
        "started_at": datetime.now().isoformat()
    }


def save_checkpoint(checkpoint):
    """Save checkpoint."""
    checkpoint["last_updated"] = datetime.now().isoformat()
    with open(CHECKPOINT_FILE, 'w') as f:
        json.dump(checkpoint, f, indent=2)


def load_mesh_synonyms():
    """Load MeSH synonyms from Step 1."""
    mesh_file = config.OUTPUT_DIR / "step1_mesh_synonyms.json"
    if mesh_file.exists():
        logger.info(f"Loading MeSH synonyms from Step 1...")
        with open(mesh_file, 'r') as f:
            mesh_data = json.load(f)
        logger.info(f"  Loaded {len(mesh_data):,} MeSH terms with synonyms")
        return mesh_data
    return {}


def load_consolidated_data():
    """Load existing consolidated data."""
    if CONSOLIDATED_DATA_FILE.exists():
        logger.info(f"Loading existing consolidated data...")
        with open(CONSOLIDATED_DATA_FILE, 'r') as f:
            data = json.load(f)
        logger.info(f"  Loaded {len(data):,} unique MeSH terms")
        return data
    return {}


def save_consolidated_data(data):
    """Save consolidated data."""
    logger.info(f"Saving consolidated data...")
    with open(CONSOLIDATED_DATA_FILE, 'w') as f:
        json.dump(data, f, indent=2)

    file_size_mb = CONSOLIDATED_DATA_FILE.stat().st_size / (1024 * 1024)
    logger.info(f"  Saved {len(data):,} unique MeSH terms ({file_size_mb:.2f} MB)")


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


def process_xml_and_merge(xml_path, consolidated_data, mesh_synonyms):
    """
    Process XML and merge into consolidated data structure.

    NO DUPLICATES: If MeSH term already exists, just increment count and add texts.
    Extracts: Title + Abstract + Keywords
    Includes MeSH synonyms from Step 1.
    """
    logger.info(f"  Processing {xml_path.name}...")

    stats = {
        "total_articles": 0,
        "with_mesh": 0,
        "with_text": 0,
        "with_title": 0,
        "with_abstract": 0,
        "with_keywords": 0,
        "mesh_terms_updated": 0,
        "new_mesh_terms": 0
    }

    try:
        context = etree.iterparse(str(xml_path), events=('end',), tag='PubmedArticle')

        for event, elem in context:
            stats["total_articles"] += 1

            # Extract MeSH terms
            mesh_terms = []
            for mesh in elem.xpath('.//MeshHeading/DescriptorName'):
                if mesh.text:
                    mesh_terms.append(mesh.text)

            # Extract title
            title_elem = elem.find('.//ArticleTitle')
            title = title_elem.text if title_elem is not None and title_elem.text else None

            # Extract abstract (handle multiple AbstractText elements)
            abstract_parts = []
            for abstract_elem in elem.xpath('.//Abstract/AbstractText'):
                if abstract_elem.text:
                    abstract_parts.append(abstract_elem.text)
            abstract = ' '.join(abstract_parts) if abstract_parts else None

            # Extract keywords
            keywords = []
            for kw in elem.xpath('.//KeywordList/Keyword'):
                if kw.text:
                    keywords.append(kw.text)

            # Statistics
            if mesh_terms:
                stats["with_mesh"] += 1
            if title:
                stats["with_title"] += 1
            if abstract and len(abstract) > 100:
                stats["with_abstract"] += 1
            if keywords:
                stats["with_keywords"] += 1

            # Build combined text (Title || Abstract || Keywords)
            text_parts = []
            if title:
                text_parts.append(title)
            if abstract and len(abstract) > 100:
                text_parts.append(abstract[:500])  # First 500 chars of abstract
            if keywords:
                text_parts.append(" | ".join(keywords))  # Keywords separated by |

            combined_text = " || ".join(text_parts)  # Separator between sections

            # Only process if we have MeSH terms and some text
            if mesh_terms and combined_text:
                stats["with_text"] += 1

                # For each MeSH term in this article
                for mesh_term in mesh_terms:
                    key = mesh_term.lower()  # Normalize for consistency

                    # Create entry if doesn't exist
                    if key not in consolidated_data:
                        # Get synonyms from Step 1 MeSH data
                        mesh_syns = []
                        if key in mesh_synonyms:
                            mesh_syns = mesh_synonyms[key].get("synonyms", [])

                        consolidated_data[key] = {
                            "preferred_term": mesh_term,  # Keep original case
                            "mesh_synonyms": mesh_syns,  # Add MeSH synonyms
                            "article_count": 0,
                            "sample_texts": []  # Changed from sample_abstracts
                        }
                        stats["new_mesh_terms"] += 1

                    # Increment article count
                    consolidated_data[key]["article_count"] += 1
                    stats["mesh_terms_updated"] += 1

                    # Add combined text to samples (limit to 100 per term)
                    if len(consolidated_data[key]["sample_texts"]) < 100:
                        consolidated_data[key]["sample_texts"].append(combined_text)

            # Progress
            if stats["total_articles"] % 10000 == 0:
                logger.info(f"    Progress: {stats['total_articles']:,} articles "
                          f"({stats['with_mesh']:,} with MeSH, "
                          f"{stats['with_text']:,} with text, "
                          f"{len(consolidated_data):,} unique terms so far)")

            # Clean memory
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]

        logger.info(f"    ✓ Completed: {stats['total_articles']:,} articles")
        logger.info(f"      With title: {stats['with_title']:,}")
        logger.info(f"      With abstract: {stats['with_abstract']:,}")
        logger.info(f"      With keywords: {stats['with_keywords']:,}")
        logger.info(f"      With text (any): {stats['with_text']:,}")
        logger.info(f"      New MeSH terms: {stats['new_mesh_terms']:,}")
        logger.info(f"      Updated terms: {stats['mesh_terms_updated']:,}")
        logger.info(f"      Total unique terms: {len(consolidated_data):,}")

        return stats

    except Exception as e:
        logger.error(f"    ✗ Error processing: {e}")
        return stats


def process_file(file_num, consolidated_data, mesh_synonyms):
    """Process a single baseline file and merge into consolidated data."""
    filename = f"pubmed26n{file_num:04d}.xml.gz"

    logger.info(f"\n{'='*80}")
    logger.info(f"FILE {file_num}/1334: {filename}")
    logger.info(f"{'='*80}")

    # Paths
    gz_path = TEMP_DIR / filename
    md5_path = TEMP_DIR / f"{filename}.md5"
    xml_path = TEMP_DIR / filename.replace('.gz', '')

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

        # 5. Process XML and merge into consolidated data
        stats = process_xml_and_merge(xml_path, consolidated_data, mesh_synonyms)

        success = True
        logger.info(f"✓ Successfully processed {filename}")

    except Exception as e:
        logger.error(f"✗ Error processing {filename}: {e}")
        success = False

    finally:
        # Cleanup temp files
        logger.info(f"  Cleaning up...")
        for temp_file in [gz_path, md5_path, xml_path]:
            if temp_file.exists():
                temp_file.unlink()
        logger.info(f"    ✓ Temp files deleted")

    return success, stats


def main():
    """Main execution."""
    logger.info("=" * 80)
    logger.info("PUBMED BASELINE BATCH DOWNLOAD & CONSOLIDATION")
    logger.info("=" * 80)
    logger.info(f"Started: {datetime.now()}")
    logger.info(f"Log file: {log_file}")
    logger.info(f"Output: {CONSOLIDATED_DATA_FILE}")
    logger.info(f"Checkpoint: {CHECKPOINT_FILE}")
    logger.info("")

    # Load checkpoint
    checkpoint = load_checkpoint()
    start_from = checkpoint["processed_count"] + 1

    # Load MeSH synonyms from Step 1
    mesh_synonyms = load_mesh_synonyms()

    # Load existing consolidated data
    consolidated_data = load_consolidated_data()

    logger.info(f"Resume from: File {start_from}/1334")
    logger.info(f"Already processed: {checkpoint['processed_count']} files")
    logger.info(f"Total articles so far: {checkpoint['total_articles']:,}")
    logger.info(f"Unique MeSH terms: {len(consolidated_data):,}")
    logger.info(f"MeSH synonyms available: {len(mesh_synonyms):,}")
    logger.info("")

    start_time = time.time()

    # Process files 1 to 1334
    for file_num in range(start_from, 1335):
        success, stats = process_file(file_num, consolidated_data, mesh_synonyms)

        # Update checkpoint
        if success:
            checkpoint["last_processed"] = f"pubmed26n{file_num:04d}.xml.gz"
            checkpoint["processed_count"] = file_num
            if stats:
                checkpoint["total_articles"] += stats["total_articles"]
        else:
            checkpoint["failed"].append(f"pubmed26n{file_num:04d}.xml.gz")

        # Save checkpoint and data after each file
        save_checkpoint(checkpoint)
        save_consolidated_data(consolidated_data)

        # Calculate ETA
        elapsed = time.time() - start_time
        processed_this_session = file_num - start_from + 1
        avg_time = elapsed / processed_this_session if processed_this_session > 0 else 0
        remaining = 1334 - file_num
        eta_seconds = avg_time * remaining
        eta_hours = eta_seconds / 3600

        logger.info(f"\nProgress: {file_num}/1334 ({file_num/1334*100:.1f}%)")
        logger.info(f"Total articles: {checkpoint['total_articles']:,}")
        logger.info(f"Unique MeSH terms: {len(consolidated_data):,}")
        logger.info(f"Avg time/file: {avg_time:.1f}s")
        logger.info(f"ETA: {eta_hours:.1f} hours")

        # Small delay between files
        if file_num < 1334:
            time.sleep(2)

    # Final summary
    total_time = time.time() - start_time
    logger.info("\n" + "=" * 80)
    logger.info("BATCH PROCESSING COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Total time: {total_time/3600:.2f} hours")
    logger.info(f"Files processed: {checkpoint['processed_count']}")
    logger.info(f"Files failed: {len(checkpoint['failed'])}")
    logger.info(f"Total articles: {checkpoint['total_articles']:,}")
    logger.info(f"Unique MeSH terms: {len(consolidated_data):,}")
    logger.info(f"Output file: {CONSOLIDATED_DATA_FILE}")

    # Show top 10 most frequent terms
    logger.info("\nTop 10 most frequent MeSH terms:")
    sorted_terms = sorted(
        consolidated_data.items(),
        key=lambda x: x[1]["article_count"],
        reverse=True
    )
    for i, (key, data) in enumerate(sorted_terms[:10], 1):
        logger.info(f"  {i:2d}. {data['preferred_term']:30s} - {data['article_count']:,} articles")

    logger.info("")
    logger.info("Next: Analyze abstracts to extract text variants")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("\n\n⚠️  Interrupted by user")
        logger.info("✓ Progress saved in checkpoint")
        logger.info("✓ Run again to resume from where you left off")
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
