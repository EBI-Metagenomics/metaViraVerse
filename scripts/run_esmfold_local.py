#!/usr/bin/env python3
"""scripts/run_esmfold_local.py

Local ESMFold inference via HuggingFace transformers.
Requires: torch>=2.0, transformers>=4.36, accelerate>=0.24

GPU strongly recommended (~16 GB VRAM; A100/V100).
  GPU:  ~1–4 s per protein
  CPU:  ~2–5 min per protein (only for very small test sets)

First run downloads ~2.5 GB model weights to HF_HOME.
Set HF_HOME to scratch space on HPC to avoid home quota issues:
  export HF_HOME=/scratch/$USER/hf_cache

pLDDT from outputs.plddt is on 0–100 scale; normalised to 0–1 in output TSV.

Usage:
    python run_esmfold_local.py \
        --input  filtered_reps.faa \
        --outdir structures/ \
        --summary structure_confidence.tsv
"""

import argparse
import re
import sys
from pathlib import Path


def parse_fasta(path: Path) -> dict:
    seqs, current = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def load_model():
    import torch
    from transformers import EsmForProteinFolding, AutoTokenizer
    print("Loading ESMFold model (first run downloads ~2.5 GB weights)...",
          file=sys.stderr)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model     = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", low_cpu_mem_usage=True
    )
    device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cpu":
        print("WARNING: No GPU detected — inference will be very slow.",
              file=sys.stderr)
    model = model.to(device)
    model.eval()
    print(f"  Model ready on {device}", file=sys.stderr)
    return tokenizer, model, device


def predict(seq_id: str, seq: str, tokenizer, model, device) -> tuple[str, float]:
    """Returns (pdb_str, mean_plddt_0_to_1)."""
    import torch
    # Replace non-standard amino acids with X (HF tokeniser requirement)
    seq_clean = re.sub(r"[UZOB]", "X", seq)
    inputs = tokenizer([seq_clean], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        outputs = model(**inputs)
    pdb_str    = model.output_to_pdb(outputs)[0]
    mean_plddt = outputs.plddt.mean().item() / 100.0  # normalise 0–100 → 0–1
    return pdb_str, mean_plddt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",   required=True)
    parser.add_argument("--outdir",  required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(Path(args.input))

    if not seqs:
        print("No sequences in input — writing empty summary.", file=sys.stderr)
        Path(args.summary).write_text("seq_id\tlength\tmean_plddt\tstatus\n")
        return

    tokenizer, model, device = load_model()

    with open(args.summary, "w") as fh:
        fh.write("seq_id\tlength\tmean_plddt\tstatus\n")
        for seq_id, seq in seqs.items():
            pdb_path = outdir / f"{seq_id}.pdb"
            if pdb_path.exists():
                print(f"  Skipping {seq_id} (already predicted)")
                continue
            try:
                pdb_str, plddt = predict(seq_id, seq, tokenizer, model, device)
                pdb_path.write_text(pdb_str)
                fh.write(f"{seq_id}\t{len(seq)}\t{plddt:.3f}\tsuccess\n")
                print(f"  {seq_id}: pLDDT={plddt:.3f}")
            except Exception as e:
                print(f"  FAILED {seq_id}: {e}", file=sys.stderr)
                fh.write(f"{seq_id}\t{len(seq)}\tNA\tfailed\n")


if __name__ == "__main__":
    main()
