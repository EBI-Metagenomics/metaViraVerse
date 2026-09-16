#!/usr/bin/env Rscript

# Shared helpers for the plasquid_* scripts (modules/plasquid/*.nf).
#
# Meant to be `source()`d, not executed directly. Each plasquid_*.R script
# locates this file next to itself (it lives alongside them on $PATH, as
# staged by Nextflow's bin/ handling) and sources it -- see
# `source_plasquid_utils()` usage at the top of every plasquid_*.R script.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))

#' Locate and source a sibling script next to the currently-running Rscript.
#'
#' Works regardless of how the script was invoked (bare name resolved via
#' PATH, relative path, or absolute path) because it reads the running
#' script's own resolved path from `--file=` rather than relying on the
#' working directory or a PATH lookup of the sibling file itself.
source_sibling <- function(filename) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  this_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1])))
  } else {
    "." # interactive / Rscript -e session, e.g. when unit-testing helpers directly
  }
  source(file.path(this_dir, filename))
}

#' Column names for a HMMER3 `--domtblout` file, in on-disk field order.
#' See http://eddylab.org/software/hmmer/Userguide.pdf (domtblout format).
HMMER_DOMTBLOUT_COLS <- c(
  "target_name", "target_accession", "tlen",
  "query_name", "query_accession", "qlen",
  "evalue", "score", "bias",
  "dom_num", "dom_of",
  "dom_c_evalue", "dom_i_evalue", "dom_score", "dom_bias",
  "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
  "acc"
)
HMMER_DOMTBLOUT_NUMERIC_COLS <- c(
  "tlen", "qlen", "evalue", "score", "bias",
  "dom_num", "dom_of", "dom_c_evalue", "dom_i_evalue", "dom_score", "dom_bias",
  "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"
)

#' Read a HMMER3 `--domtblout` file into a correctly-typed tibble.
#'
#' Handles the two footguns that the original per-module parsers each hit:
#'  - scores/e-values/etc. are parsed as text and never coerced to numeric,
#'    so `>=` cutoff comparisons silently become lexicographic string
#'    comparisons (e.g. "9.1" >= "61" is TRUE) instead of numeric ones;
#'  - a domtblout with zero passing hits (all-comment, i.e. no data lines --
#'    a normal outcome, not an error condition) used to crash `colnames<-`
#'    with a column-count mismatch instead of yielding an empty table.
#'
#' @param path Path to the domtblout file.
#' @param keep_description If TRUE, also return the free-text "description
#'   of target" field (columns 23+, whitespace-collapsed into one string).
#' @return Tibble with `HMMER_DOMTBLOUT_COLS` (plus `description` if
#'   requested), one row per domain hit. Zero rows (not an error) when the
#'   file has no hits.
read_hmmer_domtblout <- function(path, keep_description = FALSE) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]

  out_cols <- HMMER_DOMTBLOUT_COLS
  if (keep_description) out_cols <- c(out_cols, "description")

  if (length(lines) == 0) {
    empty <- as_tibble(setNames(
      rep(list(character(0)), length(out_cols)), out_cols
    ))
    return(coerce_domtblout_types(empty))
  }

  fields <- strsplit(trimws(lines), "\\s+")
  rows <- lapply(fields, function(x) {
    if (length(x) < 22) {
      stop("Malformed domtblout line in ", path, ": expected at least 22 fields, got ", length(x))
    }
    fixed <- x[1:22]
    if (keep_description) {
      description <- if (length(x) >= 23) paste(x[23:length(x)], collapse = " ") else ""
      c(fixed, description)
    } else {
      fixed
    }
  })

  tab <- as_tibble(
    as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE),
    .name_repair = "minimal"
  )
  colnames(tab) <- out_cols
  coerce_domtblout_types(tab)
}

coerce_domtblout_types <- function(tab) {
  tab %>% mutate(across(any_of(HMMER_DOMTBLOUT_NUMERIC_COLS), as.numeric))
}

#' Read a headerless whitespace-delimited table (e.g. Infernal --tblout),
#' given explicit column names.
#'
#' Unlike `read_table()`/`read.table()` alone, this stays correctly shaped
#' (all requested columns, zero rows) when the file has no data rows -- with
#' nothing to sample, `read_table()` would otherwise guess just one column,
#' and assigning the real (longer) set of column names to it errors out.
#'
#' @param path Path to the file.
#' @param col_names Column names, in on-disk field order.
#' @param comment Comment-line prefix to skip (default "#", as used by both
#'   HMMER and Infernal output).
read_whitespace_table <- function(path, col_names, comment = "#") {
  lines <- readLines(path, warn = FALSE)
  has_data <- any(!grepl(paste0("^\\s*", comment), lines) & nzchar(trimws(lines)))

  if (!has_data) {
    return(as_tibble(setNames(rep(list(character(0)), length(col_names)), col_names)))
  }

  tab <- read_table(path, comment = comment, col_names = FALSE, show_col_types = FALSE)
  colnames(tab) <- col_names
  tab
}

#' Contig/sequence id from a *nucleotide/RNA* hit id, e.g. a cmsearch target
#' name that already *is* the (possibly renamed) contig id. NOT safe for
#' protein/ORF ids -- see `protein_to_contig()` below.
orf_to_contig <- function(x) sub("_.*", "", x)

#' Read the (protein_id, contig) crosswalk produced by the
#' MAP_PROTEIN_TO_CONTIG module (derived from the representative GFF).
read_protein_contig_map <- function(path) {
  read_tsv(path, show_col_types = FALSE)
}

#' Resolve protein/ORF hit ids to their parent contig id via an explicit
#' crosswalk (see `read_protein_contig_map()`), instead of string-stripping
#' the id itself.
#'
#' This pipeline's protein FASTA/domtblout ids keep their original
#' (pre-rename) form, e.g. "MGYG000517684_26|plasmid-1:6000_1", while the
#' matching nucleotide contig has since been renamed to a short accession,
#' e.g. "seq15" -- the two no longer share a prefix, so
#' `sub("_.*", "", protein_id)` (i.e. `orf_to_contig()`) cannot recover it;
#' it silently returns a wrong/unmatched id, which downstream turns into
#' dropped rows or a hard crash indexing the real assembly by that id.
#'
#' Ids with no entry in the map resolve to NA rather than a guess.
protein_to_contig <- function(protein_ids, map) {
  map$contig[match(protein_ids, map$protein_id)]
}

#' First whitespace-delimited token of a FASTA header/description, i.e. the
#' sequence id Biostrings would use. Vectorised equivalent of the
#' `strsplit(x, " ")[[1]][1]` loops previously repeated across several
#' plasquid_*.R scripts.
fasta_id <- function(x) sub(" .*", "", x)
