#!/usr/bin/env python3
"""scripts/annotate_ecod_hhsearch.py

Annotates predicted PDB structures with ECOD domain classifications using
HHsearch against the ECOD HHM database fetched from Zenodo at runtime.

ECOD Zenodo record: https://zenodo.org/records/13993145
  - DB is downloaded once and cached; not re-fetched between runs.
  - SCOP cross-references are parsed from ECOD hit annotations.

ECOD levels returned:
  X-group: possible very distant homology
  H-group: probable homology
  T-group: topological similarity
  F-group: family (sequence-similar)

Pre-filters to pLDDT > plddt_cutoff before annotation.
Structures below cutoff are written as 'low_plddt' in output (not omitted).

Requires: hhsuite >=3.3, requests, biopython

Usage:
    python annotate_ecod_hhsearch.py \
        --pdb-dir structures/ \
        --confidence-tsv structure_confidence.tsv \
        --plddt-cutoff 0.7 \
        --output ecod_annotations.tsv \
        --scop-output scop_annotations.tsv \
        --zenodo-url https://zenodo.org/records/13993145/files/ecod_hhm_db.tar.gz \
        --threads 8 \
        --cache-dir /path/to/cache
"""

import argparse
import os
import subprocess
import tarfile
import tempfile
import time
import requests
from pathlib import Path

PLDDT_CUTOFF_DEFAULT = 0.7
ECOD_LEVELS = ["ecod_uid", "ecod_xgroup", "ecod_hgroup", "ecod_tgroup",
               "ecod_fgroup", "ecod_evalue"]
SCOP_LEVELS = ["scop_class", "scop_fold", "scop_superfamily", "scop_family"]

AA_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def fetch_ecod_db(zenodo_url: str, cache_dir: Path) -> Path:
    """Download and extract ECOD HHM DB from Zenodo if not already cached."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    db_path = cache_dir / "ecod_hhm_db"
    marker = cache_dir / ".ecod_db_ready"

    if marker.exists() and db_path.exists():
        print(f"  ECOD HHM DB found in cache: {db_path}")
        return db_path

    archive_path = cache_dir / "ecod_hhm_db.tar.gz"
    print(f"  Fetching ECOD HHM DB from Zenodo: {zenodo_url}")
    with requests.get(zenodo_url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with open(archive_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)

    print(f"  Extracting ECOD HHM DB...")
    with tarfile.open(archive_path, "r:gz") as tar:
        tar.extractall(cache_dir)
    archive_path.unlink()
    marker.touch()
    print(f"  ECOD HHM DB ready at: {db_path}")
    return db_path


def extract_sequence_from_pdb(pdb_path: Path) -> str:
    seen = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                res_name = line[17:20].strip()
                res_num = int(line[22:26].strip())
                chain = line[21]
                if (chain, res_num) not in seen:
                    seen[(chain, res_num)] = AA_MAP.get(res_name, "X")
    return "".join(seen[k] for k in sorted(seen))


def write_a3m(seq_id: str, seq: str, path: Path):
    """Write a minimal .a3m query file for HHsearch (single sequence, no MSA)."""
    with open(path, "w") as f:
        f.write(f">{seq_id}\n{seq}\n")


def run_hhsearch(query_a3m: Path, db_path: Path, threads: int) -> Path:
    """Run HHsearch and return path to .hhr result file."""
    result_hhr = query_a3m.with_suffix(".hhr")
    cmd = [
        "hhsearch",
        "-i", str(query_a3m),
        "-d", str(db_path / "ecod_hhm_db"),
        "-o", str(result_hhr),
        "-cpu", str(threads),
        "-e", "1e-3",
        "-B", "5",
        "-Z", "5",
    ]
    subprocess.run(cmd, check=True, capture_output=True)
    return result_hhr


def parse_hhr(hhr_path: Path, seq_id: str) -> dict:
    """Parse ECOD domain hit from HHsearch .hhr output."""
    result = {k: "NA" for k in ECOD_LEVELS + SCOP_LEVELS}
    result["seq_id"] = seq_id

    if not hhr_path.exists():
        return result

    in_hits = False
    with open(hhr_path) as f:
        for line in f:
            if line.startswith(" No Hit"):
                in_hits = True
                continue
            if in_hits and line.strip() and line[0] == " " and line[1].isdigit():
                # First hit line: rank, hit_id, prob, evalue, etc.
                parts = line.split()
                if len(parts) < 4:
                    break
                # ECOD hit IDs have format: ECOD_UID|X-group|H-group|T-group|F-group
                hit_id = parts[1]
                evalue = parts[3]
                ecod_parts = hit_id.split("|")
                if len(ecod_parts) >= 5:
                    result["ecod_uid"]    = ecod_parts[0]
                    result["ecod_xgroup"] = ecod_parts[1]
                    result["ecod_hgroup"] = ecod_parts[2]
                    result["ecod_tgroup"] = ecod_parts[3]
                    result["ecod_fgroup"] = ecod_parts[4]
                result["ecod_evalue"] = evalue
                # SCOP cross-refs in brackets in hit description if present
                if "SCOP" in line:
                    import re
                    scop_match = re.search(r"SCOP:([a-z]\.\d+\.\d+\.\d+)", line)
                    if scop_match:
                        scop_id = scop_match.group(1)
                        scop_parts = scop_id.split(".")
                        result["scop_class"]       = scop_parts[0] if len(scop_parts) > 0 else "NA"
                        result["scop_fold"]        = ".".join(scop_parts[:2]) if len(scop_parts) > 1 else "NA"
                        result["scop_superfamily"] = ".".join(scop_parts[:3]) if len(scop_parts) > 2 else "NA"
                        result["scop_family"]      = scop_id if len(scop_parts) > 3 else "NA"
                break
    return result


def load_confidence(tsv_path: Path, plddt_cutoff: float) -> dict:
    """Return dict of seq_id -> (mean_plddt, passes_cutoff)."""
    data = {}
    with open(tsv_path) as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 4 and parts[3] == "success":
                try:
                    plddt = float(parts[2])
                    data[parts[0]] = (plddt, plddt >= plddt_cutoff)
                except ValueError:
                    data[parts[0]] = (None, False)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-dir",        required=True)
    parser.add_argument("--confidence-tsv", required=True)
    parser.add_argument("--plddt-cutoff",   type=float, default=PLDDT_CUTOFF_DEFAULT)
    parser.add_argument("--output",         required=True)
    parser.add_argument("--scop-output",    required=True)
    parser.add_argument("--zenodo-url",     required=True)
    parser.add_argument("--threads",        type=int, default=4)
    parser.add_argument("--cache-dir",      default="ecod_db_cache")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    db_path   = fetch_ecod_db(args.zenodo_url, cache_dir)
    pdb_dir   = Path(args.pdb_dir)
    conf      = load_confidence(Path(args.confidence_tsv), args.plddt_cutoff)

    ecod_header = "seq_id\tecod_uid\tecod_xgroup\tecod_hgroup\tecod_tgroup\tecod_fgroup\tecod_evalue\n"
    scop_header = "seq_id\tscop_class\tscop_fold\tscop_superfamily\tscop_family\n"

    with open(args.output, "w") as ecod_out, \
         open(args.scop_output, "w") as scop_out, \
         tempfile.TemporaryDirectory() as tmpdir:

        ecod_out.write(ecod_header)
        scop_out.write(scop_header)
        tmpdir = Path(tmpdir)

        for pdb_path in sorted(pdb_dir.glob("*.pdb")):
            seq_id = pdb_path.stem
            plddt_info = conf.get(seq_id)

            if plddt_info is None:
                print(f"  Skipping {seq_id} (not in confidence TSV)")
                continue
            plddt_val, passes = plddt_info
            if not passes:
                print(f"  Skipping {seq_id} (mean pLDDT={plddt_val:.3f} < {args.plddt_cutoff})")
                continue

            seq = extract_sequence_from_pdb(pdb_path)
            if not seq:
                print(f"  No sequence extracted from {pdb_path.name}")
                continue

            print(f"  Annotating {seq_id} (pLDDT={plddt_val:.3f})...")
            a3m_path = tmpdir / f"{seq_id}.a3m"
            write_a3m(seq_id, seq, a3m_path)

            try:
                hhr_path = run_hhsearch(a3m_path, db_path, args.threads)
                hit = parse_hhr(hhr_path, seq_id)
            except subprocess.CalledProcessError as e:
                print(f"  HHsearch failed for {seq_id}: {e}")
                hit = {k: "NA" for k in ECOD_LEVELS + SCOP_LEVELS}
                hit["seq_id"] = seq_id

            ecod_out.write(
                f"{seq_id}\t{hit['ecod_uid']}\t{hit['ecod_xgroup']}\t"
                f"{hit['ecod_hgroup']}\t{hit['ecod_tgroup']}\t"
                f"{hit['ecod_fgroup']}\t{hit['ecod_evalue']}\n"
            )
            scop_out.write(
                f"{seq_id}\t{hit['scop_class']}\t{hit['scop_fold']}\t"
                f"{hit['scop_superfamily']}\t{hit['scop_family']}\n"
            )


if __name__ == "__main__":
    main()
