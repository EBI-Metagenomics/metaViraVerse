#!/usr/bin/env Rscript
#
# Usage: plasquid_dom_arch.R <domtblout_file>
#
# First stage of REPSEARCH's replicon-domain-architecture analysis: reads the
# hmmsearch --domtblout hits of predicted proteins against the RIP
# (replication initiator protein) profile database, and for each candidate
# RIP works out whether it hit a single domain or several. For multi-domain
# hits, overlapping domains are resolved by keeping the higher-scoring one at
# each position, producing an ordered domain "architecture" per RIP that
# `plasquid_filter_rip.R` later matches against a reference architecture list.
#
# Output: multi_dom_rip.tsv, single_dom_rip.tsv (feeds plasquid_filter_rip.R),
# domain_architecture.RDS (a named list of RIP -> architecture, also consumed
# by plasquid_filter_rip.R).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: plasquid_dom_arch.R <domtblout_file>")
}
domtblout_file <- args[1]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
.this_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
})
source(file.path(.this_dir, "plasquid_utils.R"))

# ------------------------------------------------------------
# Read hmmsearch --domtblout (empty input -> zero-row tibble, not an error:
# finding no RIP domains at all in a batch is a valid outcome, not a failure)
# ------------------------------------------------------------

tib <- read_hmmer_domtblout(domtblout_file) %>%
  rename(
    RIP        = target_name,
    queryname  = query_name,
    qaccession = query_accession,
    num        = dom_num,
    of         = dom_of,
    Evalue     = evalue,
    hmmfrom    = hmm_from,
    hmmto      = hmm_to,
    alifrom    = ali_from,
    alito      = ali_to
  )

# ------------------------------------------------------------
# Single-domain output: RIPs with exactly one domain hit
# ------------------------------------------------------------

single_dom <- tib %>%
  filter(num == 1, of == 1) %>%
  transmute(RIP, tlen, queryname, qaccession, qlen, Evalue, score, hmmfrom, hmmto, alifrom, alito)

# ------------------------------------------------------------
# Resolve a multi-domain RIP's architecture: order domains by alignment
# start, and where two domains overlap keep only the higher-scoring one.
# ------------------------------------------------------------

resolve_architecture <- function(df) {
  df <- df %>% arrange(alifrom, alito)
  n <- nrow(df)

  if (n == 1) {
    return(c(as.character(df$queryname[1]), "single_domain"))
  }

  architecture <- character(0)
  i <- 1
  while (i <= n) {
    current <- df[i, ]

    if (i == n) {
      architecture <- c(architecture, as.character(current$queryname))
      break
    }

    next_domain <- df[i + 1, ]
    overlap <- current$alito >= next_domain$alifrom

    if (!overlap) {
      architecture <- c(architecture, as.character(current$queryname), "not_over")
      i <- i + 1
    } else {
      winner <- if (current$score >= next_domain$score) current else next_domain
      architecture <- c(architecture, as.character(winner$queryname), "overlapped")
      i <- i + 2 # skip both overlapping domains
    }
  }

  architecture
}

# ------------------------------------------------------------
# Process each RIP with more than one domain hit
# ------------------------------------------------------------

rip_ids <- unique(tib$RIP)
architectures <- vector("list", length(rip_ids))
names(architectures) <- rip_ids
multi_dom_rip <- character(0)

for (ri in rip_ids) {
  tab1 <- tib %>% filter(RIP == ri) %>% arrange(alifrom, alito)
  if (nrow(tab1) <= 1) next

  multi_dom_rip <- c(multi_dom_rip, as.character(ri))
  architectures[[ri]] <- resolve_architecture(tab1)
}
architectures <- architectures[names(architectures) %in% multi_dom_rip]

# ------------------------------------------------------------
# Output (always written, even when no RIP hits were found at all)
# ------------------------------------------------------------

saveRDS(architectures, "domain_architecture.RDS")
write_delim(tibble(RIP = multi_dom_rip), "multi_dom_rip.tsv", delim = "\t")
write_delim(single_dom, "single_dom_rip.tsv", delim = "\t")
