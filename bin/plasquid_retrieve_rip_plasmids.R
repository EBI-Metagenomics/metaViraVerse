#!/usr/bin/env Rscript
#
# Usage: plasquid_retrieve_rip_plasmids.R <filtered_classif.tsv> <rep_domains.tsv> \
#          <mob_table.tsv> <assembly.fna>
#
# Combines the three independent plasmid-evidence tables -- INCSEARCH's
# Inc-group classification, REPSEARCH's replicon (RIP) domain calls, and
# MOBSEARCH's mobilisation-gene calls -- into one per-contig report, and
# extracts the nucleotide sequence of every contig with at least one hit
# from any of the three.
#
# Output: plasmids_contigs.fasta, plasmid_report.tsv (Contig, RIP_domain,
# MOB_group, Inc_group, contig_length)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plasquid_retrieve_rip_plasmids.R <filtered_classif.tsv> <rep_domains.tsv> <mob_table.tsv> <assembly.fna>")
}
inc_classif_file <- args[1]
rep_domains_file <- args[2]
mob_table_file   <- args[3]
assembly_file    <- args[4]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(Biostrings))
.this_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
})
source(file.path(.this_dir, "plasquid_utils.R"))

inc_classif <- read_tsv(inc_classif_file, show_col_types = FALSE)
rep_domains <- na.omit(read_tsv(rep_domains_file, show_col_types = FALSE))
mob_table   <- read_tsv(mob_table_file, show_col_types = FALSE)
assembly    <- readDNAStringSet(assembly_file)

# ------------------------------------------------------------
# One row per contig per evidence source, its hits' family/domain names
# comma-joined (a contig can carry more than one Rep/Inc/MOB call).
# ------------------------------------------------------------

per_contig_label <- function(tab, label_col, out_name) {
  if (nrow(tab) == 0) return(tibble(contig = character(0), "{out_name}" := character(0)))
  tab %>%
    group_by(contig) %>%
    summarise("{out_name}" := paste(.data[[label_col]], collapse = ","), .groups = "drop")
}

rip_domain <- per_contig_label(rep_domains, "Rep_type",   "RIP_domain")
mob_group  <- per_contig_label(mob_table,   "query_name", "MOB_group")
inc_group  <- per_contig_label(inc_classif, "Inc_det",    "Inc_group")

report <- purrr::reduce(list(rip_domain, mob_group, inc_group), full_join, by = "contig") %>%
  dplyr::rename(Contig = contig)

# ------------------------------------------------------------
# Pull out each reported contig's nucleotide sequence and length
# ------------------------------------------------------------

seq_id     <- fasta_id(names(assembly))
seq_contig <- orf_to_contig(seq_id)
names(assembly) <- seq_id

plasmid_seqs <- assembly[seq_contig %in% report$Contig]
plasmid_contigs <- orf_to_contig(names(plasmid_seqs))

report <- report %>%
  mutate(contig_length = purrr::map_chr(Contig, function(ct) {
    lens <- width(plasmid_seqs)[plasmid_contigs == ct]
    if (length(lens) == 0) NA_character_ else paste(lens, collapse = ",")
  }))

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

writeXStringSet(plasmid_seqs, "plasmids_contigs.fasta")
write_delim(report, "plasmid_report.tsv", delim = "\t")
