#!/usr/bin/env python3
"""
Versión PARALELA de batch_download_consolidated.py
Usa ProcessPoolExecutor para descargar y procesar múltiples XMLs en paralelo.

Diferencias respecto al original:
  - MAX_WORKERS = 50 (configurable via --workers)
  - MAX_SAMPLES  = 1000 por término (configurable via --samples)
  - Cada worker tiene su propio temp dir para evitar colisiones
  - El merge al dict compartido lo hace el proceso principal (thread-safe)
  - Guarda checkpoint cada 50 archivos procesados (no en cada uno)
  - Retry automático con backoff exponencial si NCBI rate-limita

Uso:
  python3 batch_download_parallel.py               # 50 workers, 1000 samples
  python3 batch_download_parallel.py --workers 20  # 20 workers
  python3 batch_download_parallel.py --test 10     # solo los primeros 10 archivos

Output:
  output/step2_mesh_consolidated.json   (mismo formato que el original)
  output/consolidated_checkpoint.json   (checkpoint de progreso)
"""

import os
import sys
import json
import gzip
import hashlib
import logging
import argparse
import requests
import time
import traceback
import random
from pathlib import Path
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager, Value
from lxml import etree
import config

# ──────────────────────────────────────────────────────────────────────────────
# CONFIGURACIÓN
# ──────────────────────────────────────────────────────────────────────────────
NCBI_FTP_BASE      = "https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/"
TEMP_DIR           = config.BASE_DIR / "temp_downloads"

# Archivos de salida con sufijo distintivo para no sobreescribir el análisis
# original (batch_download_consolidated.py → 100 samples, secuencial).
# Este análisis: 50 workers en paralelo + 1000 samples por término.
CHECKPOINT_FILE    = config.OUTPUT_DIR / "consolidated_checkpoint_parallel_1000s.json"
CONSOLIDATED_FILE  = config.OUTPUT_DIR / "step2_mesh_consolidated_parallel_1000s.json"

DEFAULT_WORKERS    = 20   # OPTIMIZADO: Balance entre velocidad y estabilidad
DEFAULT_SAMPLES    = 1000
CHECKPOINT_EVERY   = 50   # Guardar cada N archivos completados
DOWNLOAD_TIMEOUT   = 600  # segundos
MAX_RETRIES        = 5
RETRY_BACKOFF      = [5, 15, 30, 60, 120]  # segundos entre reintentos

# Rate limiting: 1s delay POR WORKER
# Con 20 workers: cada uno descarga cada ~1-1.5s → ~20 requests/segundo sostenidas
# Esto es SEGURO para NCBI y mucho más rápido que 5 workers
RATE_LIMIT_DELAY = 1.0    # 1 segundo base por worker


# ──────────────────────────────────────────────────────────────────────────────
# LOGGING  (configurado en el proceso principal; workers usan print)
# ──────────────────────────────────────────────────────────────────────────────
def setup_logging():
    log_file = config.OUTPUT_DIR / f"batch_parallel_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger("parallel_downloader"), log_file


# ──────────────────────────────────────────────────────────────────────────────
# CHECKPOINT
# ──────────────────────────────────────────────────────────────────────────────
def load_checkpoint():
    if CHECKPOINT_FILE.exists():
        with open(CHECKPOINT_FILE) as f:
            return json.load(f)
    return {
        "last_processed": None,
        "processed_count": 0,
        "failed": [],
        "total_articles": 0,
        "started_at": datetime.now().isoformat(),
    }


def save_checkpoint(checkpoint):
    checkpoint["last_updated"] = datetime.now().isoformat()
    with open(CHECKPOINT_FILE, "w") as f:
        json.dump(checkpoint, f, indent=2)


# ──────────────────────────────────────────────────────────────────────────────
# MESH SYNONYMS
# ──────────────────────────────────────────────────────────────────────────────
def load_mesh_synonyms():
    mesh_file = config.OUTPUT_DIR / "step1_mesh_synonyms.json"
    if mesh_file.exists():
        with open(mesh_file) as f:
            return json.load(f)
    return {}


# ──────────────────────────────────────────────────────────────────────────────
# CONSOLIDATED DATA
# ──────────────────────────────────────────────────────────────────────────────
def load_consolidated_data():
    if CONSOLIDATED_FILE.exists():
        with open(CONSOLIDATED_FILE) as f:
            return json.load(f)
    return {}


def save_consolidated_data(data, logger=None):
    msg = f"Saving {len(data):,} MeSH terms to {CONSOLIDATED_FILE.name}..."
    if logger:
        logger.info(msg)
    else:
        print(msg)
    with open(CONSOLIDATED_FILE, "w") as f:
        json.dump(data, f, indent=2)
    size_mb = CONSOLIDATED_FILE.stat().st_size / (1024 * 1024)
    msg2 = f"  ✓ Saved ({size_mb:.1f} MB)"
    if logger:
        logger.info(msg2)
    else:
        print(msg2)


# ──────────────────────────────────────────────────────────────────────────────
# HELPERS (usados dentro de los workers — sin logger global)
# ──────────────────────────────────────────────────────────────────────────────
def _md5(file_path: Path) -> str:
    h = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path, retries: int = MAX_RETRIES) -> bool:
    """
    Descarga con reintentos, backoff exponencial y rate limiting por worker.

    Rate limiting strategy:
    - Cada worker espera RATE_LIMIT_DELAY + jitter antes de descargar
    - Con 20 workers y 1s delay: ~20 requests/segundo sostenidas
    - Jitter (0-0.5s) distribuye las peticiones naturalmente
    - NCBI puede manejar ~50-100 req/s → 20 req/s es muy seguro
    """
    for attempt in range(retries):
        try:
            # Rate limiting POR WORKER con jitter para distribución
            delay = RATE_LIMIT_DELAY + random.uniform(0, 0.5)
            time.sleep(delay)

            resp = requests.get(url, stream=True, timeout=DOWNLOAD_TIMEOUT)
            resp.raise_for_status()
            with open(dest, "wb") as f:
                for chunk in resp.iter_content(chunk_size=65536):
                    if chunk:
                        f.write(chunk)
            return True
        except Exception as e:
            wait = RETRY_BACKOFF[min(attempt, len(RETRY_BACKOFF) - 1)]
            print(f"  [Worker] ⚠ Download attempt {attempt+1}/{retries} failed: {e}. "
                  f"Retrying in {wait}s...")
            time.sleep(wait)
    return False


# ──────────────────────────────────────────────────────────────────────────────
# WORKER — corre en un proceso separado, NO comparte memoria con el principal
# ──────────────────────────────────────────────────────────────────────────────
def process_file_worker(file_num: int, mesh_synonyms: dict, max_samples: int) -> tuple:
    """
    Descarga, verifica, descomprime y parsea un archivo XML de PubMed Baseline.

    Returns:
        (file_num, local_data, stats, error_msg)
        local_data: dict {term_key: {"preferred_term": ..., "mesh_synonyms": [...],
                                      "article_count": N, "sample_texts": [...]}}
        stats: dict con métricas del archivo
        error_msg: None si todo OK, string si hubo error
    """
    pid = os.getpid()
    filename = f"pubmed26n{file_num:04d}.xml.gz"

    # Temp dir exclusivo para este worker (evita colisiones entre procesos)
    worker_tmp = TEMP_DIR / f"worker_{pid}_{file_num}"
    worker_tmp.mkdir(parents=True, exist_ok=True)

    gz_path  = worker_tmp / filename
    md5_path = worker_tmp / f"{filename}.md5"
    xml_path = worker_tmp / filename.replace(".gz", "")

    local_data = {}
    stats = {
        "file_num": file_num,
        "total_articles": 0,
        "with_mesh": 0,
        "with_text": 0,
        "new_mesh_terms": 0,
    }

    try:
        # 1. Download .gz
        if not _download(NCBI_FTP_BASE + filename, gz_path):
            return file_num, {}, stats, f"Failed to download {filename}"

        # 2. Download MD5
        if not _download(NCBI_FTP_BASE + f"{filename}.md5", md5_path):
            return file_num, {}, stats, f"Failed to download MD5 for {filename}"

        # 3. Verify MD5
        with open(md5_path) as f:
            content = f.read().strip()
        expected_md5 = content.split("=")[-1].strip() if "=" in content else content.split()[0]
        actual_md5   = _md5(gz_path)
        if actual_md5 != expected_md5:
            return file_num, {}, stats, f"MD5 mismatch for {filename}"

        # 4. Decompress
        with gzip.open(gz_path, "rb") as fin:
            with open(xml_path, "wb") as fout:
                fout.write(fin.read())

        # 5. Parse XML
        context = etree.iterparse(str(xml_path), events=("end",), tag="PubmedArticle")

        for event, elem in context:
            stats["total_articles"] += 1

            # MeSH terms
            mesh_terms = [
                m.text for m in elem.xpath(".//MeshHeading/DescriptorName")
                if m.text
            ]

            # Title
            title_elem = elem.find(".//ArticleTitle")
            title = title_elem.text if title_elem is not None and title_elem.text else None

            # Abstract
            abstract_parts = [
                a.text for a in elem.xpath(".//Abstract/AbstractText") if a.text
            ]
            abstract = " ".join(abstract_parts) if abstract_parts else None

            # Keywords
            keywords = [kw.text for kw in elem.xpath(".//KeywordList/Keyword") if kw.text]

            if mesh_terms:
                stats["with_mesh"] += 1

            # Build combined text
            text_parts = []
            if title:
                text_parts.append(title)
            if abstract and len(abstract) > 100:
                text_parts.append(abstract[:500])
            if keywords:
                text_parts.append(" | ".join(keywords))

            combined_text = " || ".join(text_parts)

            if mesh_terms and combined_text:
                stats["with_text"] += 1

                for mesh_term in mesh_terms:
                    key = mesh_term.lower()

                    if key not in local_data:
                        mesh_syns = mesh_synonyms.get(key, {}).get("synonyms", [])
                        local_data[key] = {
                            "preferred_term": mesh_term,
                            "mesh_synonyms": mesh_syns,
                            "article_count": 0,
                            "sample_texts": [],
                        }
                        stats["new_mesh_terms"] += 1

                    local_data[key]["article_count"] += 1

                    # Add sample text up to max_samples limit
                    if len(local_data[key]["sample_texts"]) < max_samples:
                        local_data[key]["sample_texts"].append(combined_text)

            # Free memory
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]

        return file_num, local_data, stats, None

    except Exception as e:
        err = f"{filename}: {e}\n{traceback.format_exc()}"
        return file_num, local_data, stats, err

    finally:
        # Clean up temp files
        for f in [gz_path, md5_path, xml_path]:
            try:
                if f.exists():
                    f.unlink()
            except Exception:
                pass
        try:
            worker_tmp.rmdir()
        except Exception:
            pass


# ──────────────────────────────────────────────────────────────────────────────
# MERGE — corre en el proceso principal
# ──────────────────────────────────────────────────────────────────────────────
def merge_partial_data(consolidated: dict, partial: dict, max_samples: int) -> int:
    """
    Fusiona el dict parcial devuelto por un worker en el dict consolidado.
    Retorna el número de términos actualizados.
    """
    updated = 0
    for key, pdata in partial.items():
        if key not in consolidated:
            consolidated[key] = {
                "preferred_term": pdata["preferred_term"],
                "mesh_synonyms":  pdata["mesh_synonyms"],
                "article_count":  0,
                "sample_texts":   [],
            }

        consolidated[key]["article_count"] += pdata["article_count"]
        updated += 1

        # Append samples up to the global max_samples limit
        existing = consolidated[key]["sample_texts"]
        remaining = max_samples - len(existing)
        if remaining > 0:
            existing.extend(pdata["sample_texts"][:remaining])

    return updated


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Parallel PubMed Baseline downloader and consolidator"
    )
    parser.add_argument("--workers", type=int,  default=DEFAULT_WORKERS,
                        help=f"Number of parallel workers (default: {DEFAULT_WORKERS})")
    parser.add_argument("--samples", type=int,  default=DEFAULT_SAMPLES,
                        help=f"Max sample texts per MeSH term (default: {DEFAULT_SAMPLES})")
    parser.add_argument("--test",    type=int,  default=None, metavar="N",
                        help="Only process first N files (for testing)")
    parser.add_argument("--retry-failed", action="store_true",
                        help="Retry only the files that failed in a previous run")
    args = parser.parse_args()

    MAX_WORKERS = args.workers
    MAX_SAMPLES = args.samples
    TEST_LIMIT  = args.test
    RETRY_FAILED = args.retry_failed

    logger, log_file = setup_logging()

    logger.info("=" * 80)
    logger.info("PUBMED BASELINE — PARALLEL BATCH DOWNLOAD & CONSOLIDATION")
    logger.info("=" * 80)
    logger.info(f"Workers:      {MAX_WORKERS}")
    logger.info(f"Max samples:  {MAX_SAMPLES} per MeSH term")
    if RETRY_FAILED:
        logger.info(f"Mode:         RETRY FAILED FILES")
    elif TEST_LIMIT:
        logger.info(f"Test mode:    YES — first {TEST_LIMIT} files")
    else:
        logger.info(f"Mode:         FULL RUN (all 1334 files)")
    logger.info(f"Log file:     {log_file}")
    logger.info(f"Output:       {CONSOLIDATED_FILE}")
    logger.info("")

    # Load state
    checkpoint        = load_checkpoint()
    mesh_synonyms     = load_mesh_synonyms()
    consolidated_data = load_consolidated_data()

    # Determine which files to process
    if RETRY_FAILED:
        failed_nums = []
        for fname in checkpoint.get("failed", []):
            try:
                # Extract number from filename like pubmed26n0240.xml.gz
                num = int(fname.replace("pubmed26n", "").replace(".xml.gz", ""))
                failed_nums.append(num)
            except ValueError:
                pass
        if not failed_nums:
            logger.info("✅ No failed files in checkpoint. Nothing to retry.")
            return
        files_todo = sorted(failed_nums)
        # Clear failed list so retried files start fresh
        checkpoint["failed"] = [f for f in checkpoint["failed"]
                                 if f not in [f"pubmed26n{n:04d}.xml.gz" for n in failed_nums]]
        logger.info(f"Retrying {len(files_todo)} previously failed files: {files_todo}")
    else:
        start_from = checkpoint["processed_count"] + 1
        end_at     = (min(start_from + TEST_LIMIT - 1, 1334)
                      if TEST_LIMIT else 1334)
        files_todo = list(range(start_from, end_at + 1))

    logger.info(f"Already processed: {checkpoint['processed_count']} files")
    logger.info(f"Unique MeSH terms so far: {len(consolidated_data):,}")
    logger.info(f"MeSH synonyms loaded: {len(mesh_synonyms):,}")
    logger.info("")

    if not files_todo:
        logger.info("✅ All files already processed. Nothing to do.")
        return

    # ── Parallel execution ────────────────────────────────────────────────────
    TEMP_DIR.mkdir(exist_ok=True)
    start_time      = time.time()
    files_todo      = list(range(start_from, end_at + 1))
    total_todo      = len(files_todo)
    completed_count = 0
    total_articles  = checkpoint["total_articles"]

    logger.info(f"Submitting {total_todo} files to {MAX_WORKERS} workers...")

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Submit all at once — executor manages the queue
        futures = {
            executor.submit(process_file_worker, n, mesh_synonyms, MAX_SAMPLES): n
            for n in files_todo
        }

        for future in as_completed(futures):
            file_num = futures[future]
            try:
                file_num, partial_data, stats, error = future.result()
            except Exception as exc:
                logger.error(f"Worker for file {file_num} raised exception: {exc}")
                checkpoint["failed"].append(f"pubmed26n{file_num:04d}.xml.gz")
                continue

            if error:
                logger.error(f"✗ File {file_num:04d}: {error}")
                checkpoint["failed"].append(f"pubmed26n{file_num:04d}.xml.gz")
            else:
                # Merge into consolidated dict (single-threaded, safe)
                merge_partial_data(consolidated_data, partial_data, MAX_SAMPLES)
                total_articles            += stats["total_articles"]
                checkpoint["last_processed"] = f"pubmed26n{file_num:04d}.xml.gz"
                checkpoint["processed_count"] = max(
                    checkpoint["processed_count"], file_num
                )
                checkpoint["total_articles"] = total_articles

            completed_count += 1

            # Progress report
            elapsed    = time.time() - start_time
            rate       = completed_count / elapsed if elapsed > 0 else 0
            remaining  = total_todo - completed_count
            eta_secs   = remaining / rate if rate > 0 else 0

            logger.info(
                f"[{completed_count}/{total_todo}] "
                f"File {file_num:04d} ✓ | "
                f"MeSH terms: {len(consolidated_data):,} | "
                f"Articles: {total_articles:,} | "
                f"Rate: {rate:.1f} files/s | "
                f"ETA: {eta_secs/60:.1f} min"
            )

            # Save checkpoint + data periodically
            if completed_count % CHECKPOINT_EVERY == 0 or completed_count == total_todo:
                save_checkpoint(checkpoint)
                save_consolidated_data(consolidated_data, logger)

    # ── Final save ────────────────────────────────────────────────────────────
    save_checkpoint(checkpoint)
    save_consolidated_data(consolidated_data, logger)

    total_time = time.time() - start_time
    logger.info("")
    logger.info("=" * 80)
    logger.info("BATCH PROCESSING COMPLETE")
    logger.info("=" * 80)
    logger.info(f"Total time:       {total_time/3600:.2f} hours ({total_time/60:.1f} min)")
    logger.info(f"Files processed:  {checkpoint['processed_count']}")
    logger.info(f"Files failed:     {len(checkpoint['failed'])}")
    logger.info(f"Total articles:   {total_articles:,}")
    logger.info(f"Unique MeSH:      {len(consolidated_data):,}")
    logger.info(f"Output:           {CONSOLIDATED_FILE}")
    logger.info(f"Workers used:     {MAX_WORKERS}")
    logger.info(f"Samples/term:     {MAX_SAMPLES}")

    if checkpoint["failed"]:
        logger.warning(f"\nFailed files ({len(checkpoint['failed'])}):")
        for f in checkpoint["failed"]:
            logger.warning(f"  - {f}")

    # Top-10 MeSH terms
    logger.info("\nTop 10 most frequent MeSH terms:")
    sorted_terms = sorted(
        consolidated_data.items(),
        key=lambda x: x[1]["article_count"],
        reverse=True,
    )
    for i, (key, data) in enumerate(sorted_terms[:10], 1):
        logger.info(f"  {i:2d}. {data['preferred_term']:35s} {data['article_count']:,} articles")

    logger.info("\nNext step: python3 step2b_find_variants.py")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠  Interrupted by user")
        print("✓ Progress was saved in checkpoint (every 50 files)")
        print("✓ Run again to resume from where you left off")
    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)
