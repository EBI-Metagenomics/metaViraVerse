#!/usr/bin/env python3
"""scripts/annotate_ecod_hhsearch.py

Annotates predicted PDB structures with ECOD domain classifications using
HHsearch against the ECOD HHM database fetched from Zenodo at runtime.

ECOD Zenodo record: https://zenodo.org/records/13993145
  - DB downloaded once and cached; not re-fetched between runs.
  - SCOP cross-references parsed from HHsearch hit descriptions.

ECOD levels:
  X-group: possible very distant homology
  H-group: probable homology
  T-group: topological similarity
  F-group: family (sequence-similar)

Structures below plddt_cutoff are written with status='low_plddt' (not skipped)
so the merge step retains full sequence coverage.

Requires: hhsuite >= 3.3, requests

Usage:
    python annotate_ecod_hhsearch.py \
        --pdb-dir        structures/ \
        --confidence-tsv structure_confidence.tsv \
        --plddt-cutoff   0.7 \
        --output         ecod_annotations.tsv \
        --scop-output    scop_annotations.tsv \
        --zenodo-url     https://zenodo.org/records/13993145/files/ecod_hhm_db.tar.gz \
        --threads        8 \
        --cache-dir      /path/to/cache
"""

import argparse
import os
import re
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

import requests

PLDDT_CUTOFF_DEFAULT = 0.7

AA_MAP = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

ECOD_COLS  = ["ecod_uid", "ecod_xgroup", "ecod_hgroup",
              "ecod_tgroup", "ecod_fgroup", "ecod_evalue"]
SCOP_COLS  = ["scop_class", "scop_fold", "scop_superfamily", "scop_family"]
NA_ECOD    = {k: "NA" for k in ECOD_COLS}
NA_SCOP    = {k: "NA" for k in SCOP_COLS}


# ── Database fetch ────────────────────────────────────────────

def fetch_ecod_db(zenodo_url: str, cache_dir: Path) -> Path:
    """Download and extract ECOD HHM DB from Zenodo if not cached."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    marker  = cache_dir / ".ecod_db_ready"
    db_path = cache_dir / "ecod_hhm_db"

    if marker.exists() and db_path.exists():
        print(f"  ECOD HHM DB found in cache: {db_path}", file=sys.stderr)
        return db_path

    archive = cache_dir / "ecod_hhm_db.tar.gz"
    print(f"  Fetching ECOD HHM DB from Zenodo ({zenodo_url})...", file=sys.stderr)
    with requests.get(zenodo_url, stream=True, timeout=600) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with open(archive, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                fh.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    print(f"  Downloading... {pct:.1f}%", end="\r", file=sys.stderr)
    print("", file=sys.stderr)

    print("  Extracting ECOD HHM DB...", file=sys.stderr)
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(cache_dir)
    archive.unlink()
    marker.touch()
    print(f"  ECOD HHM DB ready: {db_path}", file=sys.stderr)
    return db_path


# ── PDB parsing ───────────────────────────────────────────────

def extract_sequence_from_pdb(pdb_path: Path) -> str:
    """Extract one-letter amino acid sequence from ATOM records."""
    seen = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            res_name = line[17:20].strip()
            chain    = line[21]
            try:
                res_num = int(line[22:26].strip())
            except ValueError:
                continue
            key = (chain, res_num)
            if key not in seen:
                seen[key] = AA_MAP.get(res_name, "X")
    return "".join(seen[k] for k in sorted(seen))


# ── HHsearch ─────────────────────────────────────────────────

def write_a3m(seq_id: str, seq: str, path: Path) -> None:
    """Minimal single-sequence A3M query (no MSA needed for domain scan)."""
    path.write_text(f">{seq_id}\n{seq}\n")


def run_hhsearch(query_a3m: Path, db_path: Path, threads: int) -> Path:
    hhr = query_a3m.with_suffix(".hhr")
    cmd = [
        "hhsearch",
        "-i",   str(query_a3m),
        "-d",   str(db_path / "ecod_hhm_db"),
        "-o",   str(hhr),
        "-cpu", str(threads),
        "-e",   "1e-3",
        "-B",   "5",
        "-Z",   "5",
        "-v",   "0",          # suppress verbose output
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"HHsearch failed: {result.stderr[:300]}")
    return hhr


# ── HHR parsing ───────────────────────────────────────────────

def parse_hhr(hhr_path: Path, seq_id: str) -> tuple[dict, dict]:
    """
    Parse top ECOD hit from HHsearch .hhr output.
    ECOD hit IDs have format: ECOD_UID|X-group|H-group|T-group|F-group
    Returns (ecod_dict, scop_dict).
    """
    ecod = {"seq_id": seq_id, **NA_ECOD}
    scop = {"seq_id": seq_id, **NA_SCOP}

    if not hhr_path.exists():
        return ecod, scop

    in_hits = False
    with open(hhr_path) as f:
        for line in f:
            if line.startswith(" No Hit"):
                in_hits = True
                continue
            if not in_hits:
                continue
            # Hit table lines start with rank digit after a space
            if line and line[0] == " " and len(line) > 4 and line[1].isdigit():
                parts = line.split()
                if len(parts) < 4:
                    break
                hit_id = parts[1]
                evalue = parts[3]
                # Parse ECOD pipe-delimited hit ID
                ep = hit_id.split("|")
                if len(ep) >= 5:
                    ecod.update({
                        "ecod_uid":    ep[0],
                        "ecod_xgroup": ep[1],
                        "ecod_hgroup": ep[2],
                        "ecod_tgroup": ep[3],
                        "ecod_fgroup": ep[4],
                        "ecod_evalue": evalue,
                    })
                # Parse optional SCOP cross-ref from description
                scop_m = re.search(r"SCOP:([a-z]\.\d+\.\d+\.\d+)", line)
                if scop_m:
                    sid = scop_m.group(1)
                    sp  = sid.split(".")
                    scop.update({
                        "scop_class":       sp[0]               if len(sp) > 0 else "NA",
                        "scop_fold":        ".".join(sp[:2])    if len(sp) > 1 else "NA",
                        "scop_superfamily": ".".join(sp[:3])    if len(sp) > 2 else "NA",
                        "scop_family":      sid                 if len(sp) > 3 else "NA",
                    })
                break   # only top hit needed

    return ecod, scop


# ── Confidence loading ────────────────────────────────────────

def load_confidence(tsv_path: Path, cutoff: float) -> dict:
    """Return {seq_id: (plddt_float, passes_bool)} for all 'success' rows."""
    data = {}
    with open(tsv_path) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            seq_id, _, plddt_str, status = parts[0], parts[1], parts[2], parts[3]
            if status != "success":
                continue
            try:
                plddt = float(plddt_str)
                data[seq_id] = (plddt, plddt >= cutoff)
            except ValueError:
                data[seq_id] = (None, False)
    return data


# ── Main ─────────────────────────────────────────────────────

def main() -> None:
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

    ecod_header = ("seq_id\tecod_uid\tecod_xgroup\tecod_hgroup\t"
                   "ecod_tgroup\tecod_fgroup\tecod_evalue\tstatus\n")
    scop_header = "seq_id\tscop_class\tscop_fold\tscop_superfamily\tscop_family\n"

    with open(args.output, "w")      as ecod_out, \
         open(args.scop_output, "w") as scop_out, \
         tempfile.TemporaryDirectory() as tmpdir:

        ecod_out.write(ecod_header)
        scop_out.write(scop_header)
        tmp = Path(tmpdir)

        for pdb_path in sorted(pdb_dir.glob("*.pdb")):
            seq_id     = pdb_path.stem
            plddt_info = conf.get(seq_id)

            # Sequences not in confidence TSV (prediction failed entirely)
            if plddt_info is None:
                ecod_out.write(f"{seq_id}\t" + "\t".join(["NA"] * 6) + "\tnot_predicted\n")
                scop_out.write(f"{seq_id}\t" + "\t".join(["NA"] * 4) + "\n")
                continue

            plddt_val, passes = plddt_info

            if not passes:
                print(f"  {seq_id}: pLDDT={plddt_val:.3f} < {args.plddt_cutoff} — skipping HHsearch",
                      file=sys.stderr)
                ecod_out.write(f"{seq_id}\t" + "\t".join(["NA"] * 6) + "\tlow_plddt\n")
                scop_out.write(f"{seq_id}\t" + "\t".join(["NA"] * 4) + "\n")
                continue

            seq = extract_sequence_from_pdb(pdb_path)
            if not seq:
                print(f"  {seq_id}: no ATOM records found — skipping", file=sys.stderr)
                ecod_out.write(f"{seq_id}\t" + "\t".join(["NA"] * 6) + "\tno_atoms\n")
                scop_out.write(f"{seq_id}\t" + "\t".join(["NA"] * 4) + "\n")
                continue

            print(f"  Annotating {seq_id} (len={len(seq)}, pLDDT={plddt_val:.3f})...",
                  file=sys.stderr)
            a3m = tmp / f"{seq_id}.a3m"
            write_a3m(seq_id, seq, a3m)

            try:
                hhr  = run_hhsearch(a3m, db_path, args.threads)
                ecod, scop = parse_hhr(hhr, seq_id)
                status = "annotated" if ecod["ecod_uid"] != "NA" else "no_hit"
            except RuntimeError as e:
                print(f"  HHsearch error for {seq_id}: {e}", file=sys.stderr)
                ecod  = {"seq_id": seq_id, **NA_ECOD}
                scop  = {"seq_id": seq_id, **NA_SCOP}
                status = "hhsearch_failed"

            ecod_out.write(
                f"{seq_id}\t{ecod['ecod_uid']}\t{ecod['ecod_xgroup']}\t"
                f"{ecod['ecod_hgroup']}\t{ecod['ecod_tgroup']}\t"
                f"{ecod['ecod_fgroup']}\t{ecod['ecod_evalue']}\t{status}\n"
            )
            scop_out.write(
                f"{seq_id}\t{scop['scop_class']}\t{scop['scop_fold']}\t"
                f"{scop['scop_superfamily']}\t{scop['scop_family']}\n"
            )

    print("ECOD annotation complete.", file=sys.stderr)


if __name__ == "__main__":
    main()
