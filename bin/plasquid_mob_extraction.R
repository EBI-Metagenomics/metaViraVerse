#!/usr/bin/env Rscript
#
# Usage: plasquid_mob_extraction.R <proteins.faa> <mob_table.tsv>
#
# Extracts the protein sequences accepted by plasquid_filter_mob.R and
# renames each to "<orf_id>_<MOB_family>" (e.g. "ctg1_5_MOBF").
#
# Output: mob_seqs.faa

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: plasquid_mob_extraction.R <proteins.faa> <mob_table.tsv>")
}
proteins_file  <- args[1]
mob_table_file <- args[2]

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
.this_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
})
source(file.path(.this_dir, "plasquid_utils.R"))

proteins <- readAAStringSet(proteins_file)
names(proteins) <- fasta_id(names(proteins))

mob_hits <- read_tsv(mob_table_file, show_col_types = FALSE)

# Exact-match lookup (not grep()/regex): ORF ids like "ctg1_5" are a
# substring of "ctg1_50", "ctg1_51", etc., so a substring match would
# rename every one of them identically instead of just the intended hit.
mob_seqs <- proteins[names(proteins) %in% mob_hits$Mob_det]
mob_hits <- mob_hits[match(names(mob_seqs), mob_hits$Mob_det), ]
names(mob_seqs) <- paste0(mob_hits$Mob_det, "_", mob_hits$query_name)

writeXStringSet(mob_seqs, "mob_seqs.faa")
