#!/usr/bin/env python3
"""Taxonomic assignment of viral contigs from their ViPhOG annotations.

Input: the annotation table of viral_contigs_annotation.py (one row per protein: Contig,
CDS_ID, ..., Best_hit ('No hit' without a ViPhOG hit), ..., Label = taxon of the hit).

For each contig, the labels of its proteins are resolved to lineages at the ViPhOG ranks
(genus, subfamily, family, order) in the NCBI taxonomy (ete3 database). From genus upwards,
the first rank where a taxon reaches the threshold (--taxthres x its factor) and the contig
is not larger than the taxon's mean + 2 SD CDS count gives the contig's lineage. Below the
threshold, the taxon closest to it is kept for that rank (a tie leaves the rank empty).

Each label is looked up in the NCBI database once, and the table is read in one pass, so
the run time grows linearly with the number of proteins.
"""

import argparse
import csv
import gzip
import logging
import math
import os
import re
import sys
from collections import Counter
from functools import lru_cache

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


# Some taxa are discontinued, we should exclude them (https://github.com/EBI-Metagenomics/emg-viral-pipeline/issues/113):
EXCLUDE_TAXA = {
    "Allolevivirus",
    "Autographivirinae",
    "Buttersvirus",
    "Caudovirales",
    "Chungbukvirus",
    "Incheonvirus",
    "Leviviridae",
    "Levivirus",
    "Mandarivirus",
    "Pbi1virus",
    "Phicbkvirus",
    "Radnorvirus",
    "Sitaravirus",
    "Vidavervirus",
    "Myoviridae",
    "Siphoviridae",
    "Podoviridae",
    "Viunavirus",
    "Orthohepevirus",
    "Klosneuvirus",
    "Hendrixvirus",
    "Rubulavirus",
    "Avulavirus",
    "Catovirus",
    "Nucleorhabdovirus",
    "Gammalipothrixvirus",
    "Peduovirinae",
    "Sedoreovirinae",
}

# superkingdom;kingdom;phylum;subphylum;class;order;suborder;family;subfamily;genus
TAXA_ORDER = [
    "superkingdom",
    "kingdom",
    "phylum",
    "subphylum",
    "class",
    "order",
    "suborder",
    "family",
    "subfamily",
    "genus",
]
VIPHOG_RANKS = {"genus", "subfamily", "family", "order"}

# values pandas.read_csv read as missing (the earlier version of this script used it)
MISSING_VALUES = {"", "NA", "N/A", "n/a", "NaN", "nan", "-nan", "-NaN", "null", "NULL", "None", "<NA>", "#N/A", "#NA"}


def open_text(path):
    with open(path, "rb") as fh:
        is_gzip = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", newline="") if is_gzip else open(path, newline="")


def read_factors(path):
    """taxon -> {avg_cds, std_cds, factor} (CSV with a header; name first, the last four
    columns avg_cds, std_cds, -, factor)."""
    with open_text(path) as infile:
        reader = csv.reader(infile)
        next(reader)
        return {
            name: {"avg_cds": float(avg_cds), "std_cds": float(std_cds), "factor": float(mult_factor)}
            for name, *_, avg_cds, std_cds, _, mult_factor in reader
        }


def read_annotations(path):
    """contig -> (number of proteins, number of proteins with a hit, labels of the hits), in file order."""
    contigs = {}
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = {"Contig", "Best_hit", "Label"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"Error: {path} has no column(s) {', '.join(sorted(missing))}")
        for row in reader:
            entry = contigs.setdefault(row["Contig"], [0, 0, []])
            entry[0] += 1
            if row["Best_hit"] != "No hit":
                entry[1] += 1
            label = (row["Label"] or "").strip()
            if label not in MISSING_VALUES:
                entry[2].append(label)
    return contigs


class Lineages:
    """NCBI lineages of taxon names, each name looked up once."""

    def __init__(self, ncbi, exclude):
        self.ncbi = ncbi
        self.exclude = exclude
        self.viphog_lineage = lru_cache(maxsize=None)(self._viphog_lineage)
        self.output_lineage = lru_cache(maxsize=None)(self._output_lineage)

    def _ranks_and_names(self, taxid):
        lineage = self.ncbi.get_lineage(taxid)
        return self.ncbi.get_rank(lineage), self.ncbi.get_taxid_translator(lineage)

    def _viphog_lineage(self, label):
        """{rank: name} at the ViPhOG ranks, without excluded taxa; None for an unknown label."""
        translated = self.ncbi.get_name_translator([label])
        if not translated:
            logging.warning(f"Taxon label '{label}' not found in NCBI db")
            return None
        taxid = translated[label][0]
        try:
            ranks, names = self._ranks_and_names(taxid)
        except ValueError:
            logging.warning(f"Cannot retrieve lineage for taxid {taxid} ({label})")
            return {}
        return {
            rank: names[t]
            for t, rank in ranks.items()
            if rank in VIPHOG_RANKS and names[t] not in self.exclude
        }

    def _output_lineage(self, taxon):
        """The names at the output ranks (TAXA_ORDER) of a taxon's lineage."""
        taxid = self.ncbi.get_name_translator([taxon])[taxon][0]
        ranks, names = self._ranks_and_names(taxid)
        by_rank = {rank: names[t] for t, rank in ranks.items()}
        return tuple(by_rank.get(rank, "") for rank in TAXA_ORDER)


def assign_rank(taxa, total_prot, tax_thres, factors):
    """Decision at one rank: ('lineage', taxon) to assign the contig the taxon's lineage,
    ('rank', taxon) to keep the taxon (or '') for this rank only, or None when every
    taxon over the threshold fails the CDS-size filter."""
    total_hits = len(taxa)
    under_thres = []
    over_thres = []
    for taxon, count in Counter(taxa).items():
        prop_hits = count / total_hits
        effective_thres = tax_thres * (factors[taxon]["factor"] if taxon in factors else 1)
        if prop_hits < effective_thres:
            under_thres.append((taxon, math.ceil(effective_thres * total_hits) - count))
        else:
            over_thres.append((taxon, prop_hits))

    if not over_thres:
        # All candidates fall below the threshold: keep the one closest to it (smallest
        # number of missing hits). A tie is ambiguous: leave the rank unassigned rather
        # than pick a winner by file order.
        under_thres.sort(key=lambda x: x[1])
        best_diff = under_thres[0][1]
        tied = [t for t in under_thres if t[1] == best_diff]
        return ("rank", "" if len(tied) > 1 else under_thres[0][0])

    for taxon, _ in sorted(over_thres, key=lambda x: x[1], reverse=True):
        f = factors.get(taxon)
        if f and total_prot <= f["avg_cds"] + 2 * f["std_cds"]:
            return ("lineage", taxon)
    return None


def contig_tax(contigs, lineages, tax_thres, factors):
    """Yield [contig, lineage] (';'-separated TAXA_ORDER names, '' where unassigned)."""
    logging.info(f"Assigning taxonomy for {len(contigs)} contigs")
    assigned = unassigned_no_hits = unassigned_below_thres = 0

    for contig, (total_prot, annot_prot, labels) in contigs.items():
        if annot_prot == 0:
            unassigned_no_hits += 1
            yield [contig, ";".join([""] * len(TAXA_ORDER))]
            continue

        hit_lineages = [lin for lin in map(lineages.viphog_lineage, labels) if lin]
        contig_lineage = []
        contig_assigned = False
        for rank in TAXA_ORDER[::-1]:
            taxa = [lin[rank] for lin in hit_lineages if rank in lin]
            if not taxa:
                contig_lineage.append("")
                continue
            decision = assign_rank(taxa, total_prot, tax_thres, factors)
            if decision is None:
                contig_lineage.append("")
            elif decision[0] == "rank":
                contig_lineage.append(decision[1])
            else:
                contig_lineage = list(lineages.output_lineage(decision[1]))
                contig_assigned = True
                break

        if contig_assigned:
            assigned += 1
        else:
            unassigned_below_thres += 1
        yield [contig, ";".join(contig_lineage)]

    logging.info(
        f"Summary: {assigned}/{len(contigs)} contigs assigned | {unassigned_no_hits} no ViPhOG hits | "
        f"{unassigned_below_thres} hits but unassigned"
    )


def main(argv=None, ncbi=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-d", "--db", dest="ncbi_db", help="path to the ete3 processed NCBI taxonomy db.")
    parser.add_argument("-i", "--input", dest="input_file", required=True,
                        help="Annotation table generated with viral_contigs_annotation.py (plain or gzip)")
    parser.add_argument("-f", "--factor", dest="factor", required=True,
                        help="File with the ViPhOG-specific multiplication factors and CDS counts per taxon "
                             "(viphogs_cds_per_taxon_cummulative.csv)")
    parser.add_argument("--taxthres", dest="tax_thres", type=float, default=0.6,
                        help="Minimum proportion of annotated genes required for taxonomic assignment (default: 0.6)")
    parser.add_argument("-o", "--outdir", dest="outdir", default=".",
                        help="Directory for the output file (default: cwd)")
    parser.add_argument("--version4", dest="version4", action="store_true",
                        help="Exclude the taxa that are deprecated in ViPhOGs version 4.")
    args = parser.parse_args(argv)

    logging.info(f"Processing input file: {args.input_file}")
    logging.info(f"Settings: tax_thres={args.tax_thres:.2f}, version4={args.version4}")

    if ncbi is None:
        from ete3 import NCBITaxa
        ncbi = NCBITaxa(dbfile=args.ncbi_db)
    lineages = Lineages(ncbi, EXCLUDE_TAXA if args.version4 else set())
    factors = read_factors(args.factor)
    contigs = read_annotations(args.input_file)

    out_file = re.split(r"\.[a-z]+$", os.path.basename(args.input_file))[0]
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, out_file + "_taxonomy.tsv"), "w", newline="") as output_file:
        writer = csv.writer(output_file, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["contig_ID", "taxonomy"])
        writer.writerows(contig_tax(contigs, lineages, args.tax_thres, factors))


if __name__ == "__main__":
    main()
