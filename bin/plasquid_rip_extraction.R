#!/usr/bin/env Rscript
#
# Usage: plasquid_rip_extraction.R <proteins.faa> <filtered_classif.tsv> \
#          <rep_domains.tsv> <mob_table.tsv> [mpf_contigs.tsv]
#
# Three things, all keyed off the same per-protein evidence tables:
#
#  1. Pulls out the protein sequences accepted as replication initiator
#     proteins (RIPs) by either line of evidence -- REPSEARCH's domain-based
#     calls (rep_domains.tsv) or INCSEARCH's Inc-group classification of the
#     same ORF (filtered_classif.tsv, restricted to protein-based rows:
#     RNA/whole-contig rows never carry an ORF id, so their query_name has no
#     "_" in it) -- and renames each "<Rep_ORF>#<Rep_type>" (merging repeat
#     calls for the same ORF with "#" too, e.g. two REPSEARCH domain calls
#     for the same protein).
#
#  2. Builds a per-protein evidence table (Contig, Protein, RIP_domain,
#     MOB_group, Inc_group) across all three plasquid searches -- the
#     protein-level counterpart to EXTRACT_PLASMIDS_DATA's per-contig
#     plasmid_report.tsv, useful when a contig's RIP/MOB/Inc calls sit on
#     different proteins and that distinction matters. One row per protein
#     with at least one hit from any of the three searches; a protein
#     missing a given kind of evidence gets NA in that column, and repeat
#     hits of the same kind on one protein are comma-joined.
#
#  3. Classifies every contig with plasquid evidence into the standard
#     three-tier plasmid mobility scheme and writes summary counts to
#     mobility_stats.json:
#       - conjugative:     has a relaxase (MOB) AND mating-pair-formation
#                           (MPF) / T4SS conjugation-machinery evidence
#       - mobilizable:     has a relaxase (MOB) but no MPF evidence
#       - non_mobilizable: no relaxase (MOB) detected
#     CAVEAT: plaSquid's MOBSEARCH only detects the relaxase (MOB) gene --
#     this pipeline has no MPF/T4SS-detection step, so no contig can
#     currently be placed in "conjugative" from plaSquid's own evidence
#     alone; a genuinely conjugative plasmid is reported as "mobilizable"
#     instead. The optional [mpf_contigs.tsv] argument (a one-column
#     `contig` TSV of MPF-positive contigs) lets a future caller supply that
#     missing evidence -- no current caller of this script does. MOB-suite's
#     `mob_typer` already reports this same three-tier classification
#     directly on the same representative contigs and is not MPF-blind; see
#     its output for a complete call today.
#
# Output: rip_seqs.faa, protein_report.tsv, mobility_stats.json

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: plasquid_rip_extraction.R <proteins.faa> <filtered_classif.tsv> <rep_domains.tsv> <mob_table.tsv>")
}
proteins_file    <- args[1]
inc_classif_file <- args[2]
rep_domains_file <- args[3]
mob_table_file   <- args[4]
mpf_contigs_file <- if (length(args) >= 5) args[5] else NA_character_

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
.this_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
})
source(file.path(.this_dir, "plasquid_utils.R"))

proteins <- readAAStringSet(proteins_file)
names(proteins) <- fasta_id(names(proteins))

inc_classif <- read_tsv(inc_classif_file, show_col_types = FALSE)
rep_domains <- read_tsv(rep_domains_file, show_col_types = FALSE)
mob_table   <- read_tsv(mob_table_file, show_col_types = FALSE)

# Protein-based Inc rows only (RNA/whole-contig rows have no ORF id, so their
# query_name has no "_" in it and can't be attributed to a single protein).
inc_classif_proteins <- inc_classif %>% filter(grepl("_", query_name))

# ------------------------------------------------------------
# 1. RIP protein sequences (unchanged behaviour: REPSEARCH domain calls +
#    INCSEARCH's protein-based Inc calls, both treated as RIP evidence)
# ------------------------------------------------------------

rip_calls <- bind_rows(
  rep_domains %>% transmute(contig, Rep_ORF, Rep_type),
  inc_classif_proteins %>% transmute(contig, Rep_ORF = query_name, Rep_type = Inc_det)
) %>%
  filter(!is.na(contig)) %>%
  arrange(contig, Rep_ORF) %>%
  group_by(Rep_ORF) %>%
  summarise(across(everything(), ~ paste(.x, collapse = "#")), .groups = "drop")

rip_seqs <- proteins[names(proteins) %in% rip_calls$Rep_ORF]
rip_calls_matched <- rip_calls[match(names(rip_seqs), rip_calls$Rep_ORF), ]
names(rip_seqs) <- paste(rip_calls_matched$Rep_ORF, rip_calls_matched$Rep_type, sep = "#")

writeXStringSet(rip_seqs, "rip_seqs.faa")

# ------------------------------------------------------------
# 2. Per-protein evidence table across all three searches
# ------------------------------------------------------------

per_protein_evidence <- function(tab, protein_col, value_col, out_name) {
  if (nrow(tab) == 0) {
    return(tibble(contig = character(0), Protein = character(0), "{out_name}" := character(0)))
  }
  tab %>%
    filter(!is.na(contig), !is.na(.data[[protein_col]])) %>%
    group_by(contig, Protein = .data[[protein_col]]) %>%
    summarise("{out_name}" := paste(.data[[value_col]], collapse = ","), .groups = "drop")
}

rip_evidence <- per_protein_evidence(rep_domains, "Rep_ORF", "Rep_type", "RIP_domain")
mob_evidence <- per_protein_evidence(mob_table, "Mob_det", "query_name", "MOB_group")
inc_evidence <- per_protein_evidence(inc_classif_proteins, "query_name", "Inc_det", "Inc_group")

protein_report <- purrr::reduce(
  list(rip_evidence, mob_evidence, inc_evidence),
  full_join,
  by = c("contig", "Protein")
) %>%
  dplyr::rename(Contig = contig) %>%
  arrange(Contig, Protein)

write_delim(protein_report, "protein_report.tsv", delim = "\t")

# ------------------------------------------------------------
# 3. Mobility classification + summary stats (JSON)
# ------------------------------------------------------------

#' Classify every contig with plasquid evidence into conjugative /
#' mobilizable / non_mobilizable, per the CAVEAT in the header comment.
#'
#' The contig population matches EXTRACT_PLASMIDS_DATA's plasmid_report.tsv:
#' every contig appearing in any of the three evidence tables (not just
#' MOB-positive ones), so non-mobilizable plasmids are counted too.
classify_mobility <- function(rep_domains, inc_classif, mob_table, mpf_contigs = character(0)) {
  contigs     <- unique(na.omit(c(rep_domains$contig, inc_classif$contig, mob_table$contig)))
  mob_contigs <- unique(na.omit(mob_table$contig))

  tibble(contig = contigs) %>%
    mutate(
      has_mob  = contig %in% mob_contigs,
      has_mpf  = contig %in% mpf_contigs,
      category = case_when(
        has_mob & has_mpf  ~ "conjugative",
        has_mob & !has_mpf ~ "mobilizable",
        TRUE               ~ "non_mobilizable"
      )
    )
}

#' Write classify_mobility()'s per-contig table as summary counts/percentages
#' to a JSON file (hand-rolled, not jsonlite: the shape is small and fixed,
#' and this avoids depending on a package the plasquid container may not
#' carry -- no other R script in this pipeline uses jsonlite).
write_mobility_stats_json <- function(classification, mpf_data_available, path) {
  total  <- nrow(classification)
  counts <- table(factor(classification$category, levels = c("conjugative", "mobilizable", "non_mobilizable")))
  pct    <- function(n) if (total > 0) round(100 * n / total, 2) else 0

  category_block <- function(name, note = NULL) {
    n <- unname(counts[[name]])
    fields <- sprintf('"count": %d, "percent": %s', n, pct(n))
    if (!is.null(note)) fields <- paste0(fields, sprintf(', "note": "%s"', note))
    paste0("{", fields, "}")
  }

  conjugative_note <- if (mpf_data_available) {
    NULL
  } else {
    paste(
      "plaSquid detects the MOB relaxase gene only, not MPF/T4SS conjugation",
      "machinery; no contig can be placed here without independently supplied",
      "MPF evidence. A genuinely conjugative plasmid is reported as",
      "'mobilizable' instead. See MOB-suite's mob_typer output for a complete",
      "conjugative/mobilizable/non_mobilizable call on the same representative",
      "contigs."
    )
  }

  json <- paste0(
    "{\n",
    '  "total_contigs": ', total, ",\n",
    '  "mpf_data_available": ', tolower(as.character(mpf_data_available)), ",\n",
    '  "categories": {\n',
    '    "conjugative": ',     category_block("conjugative", conjugative_note), ",\n",
    '    "mobilizable": ',     category_block("mobilizable"), ",\n",
    '    "non_mobilizable": ', category_block("non_mobilizable"), "\n",
    "  }\n",
    "}\n"
  )
  writeLines(json, path)
}

mpf_contigs <- if (!is.na(mpf_contigs_file) && file.exists(mpf_contigs_file)) {
  read_tsv(mpf_contigs_file, show_col_types = FALSE)$contig
} else {
  character(0) # no MPF-detection step exists in this pipeline yet -- see CAVEAT above
}

mobility <- classify_mobility(rep_domains, inc_classif, mob_table, mpf_contigs)
write_mobility_stats_json(mobility, mpf_data_available = length(mpf_contigs) > 0, "mobility_stats.json")
