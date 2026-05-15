#!/usr/bin/env python3
"""scripts/run_esmfold_local.py

Local ESMFold inference using HuggingFace transformers.
Requires: pip install transformers torch accelerate

Recommended for sequences >400aa or batches >100 proteins.
GPU required for practical use (~16 GB VRAM; A100/V100).
  ~1-4s per protein on GPU vs minutes on CPU.

Usage:
    python run_esmfold_local.py \
        --input filtered_reps.faa \
        --outdir structures/ \
        --summary structure_confidence.tsv
"""

import argparse
import re
from pathlib import Path
import torch


def load_model():
    from transformers import EsmForProteinFolding, AutoTokenizer
    print("Loading ESMFold model (first run downloads ~2.5 GB weights)...")
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1",
        low_cpu_mem_usage=True,
    )
    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = model.to(device)
    model.eval()
    print(f"  Model loaded on {device}")
    return tokenizer, model, device


def predict_structure(seq_id, seq, tokenizer, model, device):
    """Returns PDB string and mean pLDDT (0-1 scale)."""
    seq_clean = re.sub(r"[UZOB]", "X", seq)
    inputs = tokenizer([seq_clean], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    with torch.no_grad():
        outputs = model(**inputs)
    pdb_str = model.output_to_pdb(outputs)[0]
    # outputs.plddt is on 0-100 scale; normalise to 0-1
    mean_plddt = outputs.plddt.mean().item() / 100.0
    return pdb_str, mean_plddt


def parse_fasta(path):
    seqs = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            elif current:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    seqs = parse_fasta(Path(args.input))
    tokenizer, model, device = load_model()

    with open(args.summary, "w") as summary:
        summary.write("seq_id\tlength\tmean_plddt\tstatus\n")
        for seq_id, seq in seqs.items():
            pdb_path = outdir / f"{seq_id}.pdb"
            if pdb_path.exists():
                continue
            try:
                pdb_str, mean_plddt = predict_structure(
                    seq_id, seq, tokenizer, model, device
                )
                pdb_path.write_text(pdb_str)
                summary.write(f"{seq_id}\t{len(seq)}\t{mean_plddt:.3f}\tsuccess\n")
                print(f"  {seq_id}: pLDDT={mean_plddt:.3f}")
            except Exception as e:
                print(f"  FAILED {seq_id}: {e}")
                summary.write(f"{seq_id}\t{len(seq)}\tNA\tfailed\n")


if __name__ == "__main__":
    main()
