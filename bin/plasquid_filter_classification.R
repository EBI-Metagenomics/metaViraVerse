#!/usr/bin/env Rscript
#
# Usage: plasquid_filter_classification.R <classification_table.tsv>
#
# Resolves known cross-reactivity between Inc-group profiles: IncFII-family
# and IncZ profiles cross-react, as do Col and Inc13 profiles. When a contig
# has both members of one of these pairs, the secondary (IncZ / Inc13) call
# is considered a false positive of the primary one and dropped.
#
# NOTE: the original version of this script tested for the family-name
# substring (e.g. "IncFII") against the `query_name` column, which holds the
# hit's ORF/contig id (e.g. "ctg1_5"), not its Inc-family classification --
# so the substring never matched and this filter never actually removed
# anything. It now tests the intended column, `Inc_det`.
#
# Output: filtered_classif.tsv

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: plasquid_filter_classification.R <classification_table.tsv>")
}
classification_file <- args[1]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))

classification <- read_tsv(classification_file, show_col_types = FALSE)

#' Drop rows whose Inc_det value was seen co-occurring with `primary` on any
#' contig that also has a `secondary`-matching call (cross-reactivity: keep
#' the primary call, drop the specific Inc_det values it was seen confused
#' with, anywhere in the table -- matching the original per-value, batch-wide
#' removal rather than scoping the drop back down to just that one contig).
drop_cross_reactive <- function(tab, primary, secondary) {
  cross_reactive_values <- tab %>%
    group_by(contig) %>%
    filter(any(grepl(primary, Inc_det, fixed = TRUE))) %>%
    filter(grepl(secondary, Inc_det, fixed = TRUE)) %>%
    ungroup() %>%
    pull(Inc_det) %>%
    unique()

  tab %>% filter(!(Inc_det %in% cross_reactive_values))
}

filtered <- classification %>%
  drop_cross_reactive("IncFII", "IncZ") %>%
  drop_cross_reactive("Col", "Inc13")

write_delim(filtered, "filtered_classif.tsv", delim = "\t")
