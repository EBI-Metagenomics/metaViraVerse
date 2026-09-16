#!/usr/bin/env Rscript
#
# Usage: plasquid_rip_extraction.R <proteins.faa> <filtered_classif.tsv> <rep_domains.tsv>
#
# Pulls out the protein sequences accepted as replication initiator proteins
# (RIPs) by either line of evidence -- REPSEARCH's domain-based calls
# (rep_domains.tsv) or INCSEARCH's Inc-group classification of the same ORF
# (filtered_classif.tsv, restricted to protein-based rows: RNA/whole-contig
# rows never carry an ORF id, so their query_name has no "_" in it) -- and
# renames each "<Rep_ORF>#<Rep_type>" (merging repeat calls for the same ORF
# with "#" too, e.g. two REPSEARCH domain calls for the same protein).
#
# Output: rip_seqs.faa

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plasquid_rip_extraction.R <proteins.faa> <filtered_classif.tsv> <rep_domains.tsv>")
}
proteins_file      <- args[1]
inc_classif_file   <- args[2]
rep_domains_file   <- args[3]

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

inc_classif <- read_tsv(inc_classif_file, show_col_types = FALSE) %>%
  filter(grepl("_", query_name)) %>% # protein-based rows only (RNA/contig rows have no ORF id)
  transmute(contig, Rep_ORF = query_name, Rep_type = Inc_det)

rep_domains <- read_tsv(rep_domains_file, show_col_types = FALSE) %>%
  transmute(contig, Rep_ORF, Rep_type)

rip_calls <- bind_rows(rep_domains, inc_classif) %>%
  filter(!is.na(contig)) %>%
  arrange(contig, Rep_ORF) %>%
  group_by(Rep_ORF) %>%
  summarise(across(everything(), ~ paste(.x, collapse = "#")), .groups = "drop")

rip_seqs <- proteins[names(proteins) %in% rip_calls$Rep_ORF]
rip_calls <- rip_calls[match(names(rip_seqs), rip_calls$Rep_ORF), ]
names(rip_seqs) <- paste(rip_calls$Rep_ORF, rip_calls$Rep_type, sep = "#")

writeXStringSet(rip_seqs, "rip_seqs.faa")
