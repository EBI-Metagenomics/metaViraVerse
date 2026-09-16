#!/usr/bin/env Rscript
#
# Usage: plasquid_filter_rip.R <single_dom_rip.tsv> <domain_architecture.RDS> \
#          <repfilter_db.RDS> <protein_to_contig.tsv>
#
# Second stage of REPSEARCH: takes plasquid_dom_arch.R's per-RIP domain
# calls and applies curated, domain-specific bit-score (and, for some
# domains, length) cutoffs to decide which candidates are real replication
# initiator proteins (RIPs). Three independent lines of evidence are
# combined: single-domain RIPs passing a per-Pfam-domain score cutoff,
# single-domain RIPs passing a stricter score+length rule for domains known
# to need one, and multi-domain RIPs whose resolved domain architecture
# matches a curated reference architecture list (repfilter_db).
#
# Output: rep_domains.tsv (Rep_type, contig, Rep_ORF -- one row per accepted RIP)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plasquid_filter_rip.R <single_dom_rip.tsv> <domain_architecture.RDS> <repfilter_db.RDS> <protein_to_contig.tsv>")
}
single_dom_file      <- args[1]
architecture_file    <- args[2]
repfilter_db_file    <- args[3]
protein_contig_file  <- args[4]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
.this_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
})
source(file.path(.this_dir, "plasquid_utils.R"))

single_dom    <- read_tsv(single_dom_file, show_col_types = FALSE)
architectures <- readRDS(architecture_file)     # named list: multi-domain RIP id -> resolved architecture
reference_architectures <- readRDS(repfilter_db_file) # curated list of architectures accepted as real RIPs
protein_contig_map <- read_protein_contig_map(protein_contig_file)

empty_hits <- function() tibble(Rep_type = character(0), contig = character(0), Rep_ORF = character(0))

# ------------------------------------------------------------
# Single-domain RIPs, filtered by a fixed bit-score cutoff per Pfam domain
# ------------------------------------------------------------

single_domain_cutoffs <- tibble::tribble(
  ~queryname,       ~min_score,
  "IncFII_repA",     10.0,
  "RepA_C",           45,
  "RepA_N",           38,
  "RepC",             76,
  "Replicase",        76,
  "Rop",              77.1,
  "RPA",              67.8,
  "RP-C",             45,
  "TrfA",             87,
  "Bac_RepA_C",       30,
  "RepB-RCR_reg",     24,
  "RP-C_C",           34
)

fixed_cutoff_hits <- pmap_dfr(single_domain_cutoffs, function(queryname, min_score) {
  single_dom %>%
    filter(queryname == .env$queryname, score > .env$min_score) %>%
    transmute(Rep_type = .env$queryname, contig = protein_to_contig(RIP, protein_contig_map), Rep_ORF = RIP)
})

# ------------------------------------------------------------
# Single-domain RIPs needing a score+length rule instead of score alone
# ------------------------------------------------------------

length_gated_cutoffs <- tibble::tribble(
  ~queryname,   ~min_score, ~min_tlen, ~max_tlen,
  "PriCT_1",    49,         420,       500,
  "Rep_1",      37,         130,       Inf,
  "Rep_1",      27,         -Inf,      130,
  "Rep_3",      45,         -Inf,      Inf,
  "RepL",       85,         90,        Inf,
  "Rep_trans",  27,         -Inf,      130
)

length_gated_hits <- pmap_dfr(length_gated_cutoffs, function(queryname, min_score, min_tlen, max_tlen) {
  single_dom %>%
    filter(queryname == .env$queryname, score > .env$min_score, tlen > .env$min_tlen, tlen < .env$max_tlen) %>%
    transmute(Rep_type = .env$queryname, contig = protein_to_contig(RIP, protein_contig_map), Rep_ORF = RIP)
})

# ------------------------------------------------------------
# Multi-domain RIPs whose resolved architecture matches the reference list
# ------------------------------------------------------------

matched_architectures <- architectures[architectures %in% reference_architectures]
multi_domain_hits <- if (length(matched_architectures) > 0) {
  rip_ids <- names(matched_architectures)
  tibble(
    Rep_type = "Conserved Domain Arch",
    contig   = protein_to_contig(rip_ids, protein_contig_map),
    Rep_ORF  = rip_ids
  )
} else {
  empty_hits()
}

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

rep_domains <- bind_rows(fixed_cutoff_hits, length_gated_hits, multi_domain_hits)
write_delim(rep_domains, "rep_domains.tsv", delim = "\t")
