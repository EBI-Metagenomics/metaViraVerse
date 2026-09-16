#!/usr/bin/env Rscript
#
# Usage: plasquid_filter_mob.R <mob_candidates.tsv>
#
# Filters MOBSEARCH's hmmsearch --domtblout hits (predicted proteins against
# the MOB relaxase/mobilisation-protein profile database) down to the
# proteins that pass a curated, MOB-family-specific bit-score cutoff.
#
# Output: mob_table.tsv (Mob_det, tlen, query_name, score, alifrom, alito,
# contig -- one row per accepted protein hit)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: plasquid_filter_mob.R <mob_candidates.tsv>")
}
domtblout_file <- args[1]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
.this_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
})
source(file.path(.this_dir, "plasquid_utils.R"))

hits <- read_hmmer_domtblout(domtblout_file) %>%
  rename(Mob_det = target_name)

mob_family_cutoffs <- tibble::tribble(
  ~query_name, ~min_score,
  "MOBB",       92.7,
  "MOBC",       96.6,
  "MOBF",      332.2,
  "MOBH",       81.1,
  "MOBP1",      74.5,
  "MOBP2",     480,
  "MOBP3",     105,
  "MOBQ",       61,
  "MOBT",      115.9,
  "MOBV",       68,
  "MOBM",      584
)

mob_table <- pmap_dfr(mob_family_cutoffs, function(query_name, min_score) {
  hits %>%
    filter(.data$query_name == .env$query_name, score >= .env$min_score) %>%
    distinct(Mob_det, query_name, .keep_all = TRUE) %>%
    select(Mob_det, tlen, query_name, score, ali_from, ali_to)
}) %>%
  rename(alifrom = ali_from, alito = ali_to)

if (nrow(mob_table) > 0) {
  mob_table <- mob_table %>% mutate(contig = orf_to_contig(Mob_det))
}

write_delim(mob_table, "mob_table.tsv", delim = "\t")
