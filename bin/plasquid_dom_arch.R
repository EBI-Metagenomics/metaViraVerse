#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript script.R <domtblout_file>")
}

teb <- args[1]

library(tidyverse)


# ------------------------------------------------------------
# Read HMMER --domtblout
# ------------------------------------------------------------

# HMMER domtblout has 22 fixed fields followed by an optional
# free-text description. We only need the first 22 fields.
#
# Reading line-by-line avoids problems caused by spaces in the
# description field.

lines <- readLines(teb, warn = FALSE)

# Remove comments and empty lines
lines <- lines[
  !grepl("^\\s*#", lines) &
  nzchar(trimws(lines))
]

if (length(lines) == 0) {
  stop("No data found in: ", teb)
}

# Split on whitespace
fields <- strsplit(trimws(lines), "\\s+")

# Keep only the first 22 fields
fields <- lapply(fields, function(x) {
  if (length(x) < 22) {
    stop(
      "Malformed domtblout line: expected at least 22 fields, got ",
      length(x)
    )
  }
  x[1:22]
})

tab <- as.data.frame(
  do.call(rbind, fields),
  stringsAsFactors = FALSE
)

colnames(tab) <- c(
  "RIP",
  "taccession",
  "tlen",
  "queryname",
  "qaccession",
  "qlen",
  "Evalue",
  "score",
  "bias",
  "num",
  "of",
  "cEvalue",
  "iEvalue",
  "domscore",
  "dombias",
  "hmmfrom",
  "hmmto",
  "alifrom",
  "alito",
  "from_env",
  "to_env",
  "acc"
)


# ------------------------------------------------------------
# Convert relevant columns to appropriate types
# ------------------------------------------------------------

tib <- as_tibble(tab) %>%
  mutate(
    tlen     = as.numeric(tlen),
    qlen     = as.numeric(qlen),
    Evalue   = as.numeric(Evalue),
    score    = as.numeric(score),
    bias     = as.numeric(bias),
    num      = as.integer(num),
    of       = as.integer(of),
    cEvalue  = as.numeric(cEvalue),
    iEvalue  = as.numeric(iEvalue),
    domscore = as.numeric(domscore),
    dombias  = as.numeric(dombias),
    hmmfrom  = as.integer(hmmfrom),
    hmmto    = as.integer(hmmto),
    alifrom  = as.integer(alifrom),
    alito    = as.integer(alito),
    from_env = as.integer(from_env),
    to_env   = as.integer(to_env),
    acc      = as.numeric(acc)
  )


# ------------------------------------------------------------
# Single-domain output
# ------------------------------------------------------------

single_dom <- tib %>%
  filter(num == 1, of == 1) %>%
  transmute(
    RIP,
    tlen,
    queryname,
    qaccession,
    qlen,
    Evalue,
    score,
    hmmfrom,
    hmmto,
    alifrom,
    alito
  )


# ------------------------------------------------------------
# Function to resolve overlapping domains
# ------------------------------------------------------------

resolve_architecture <- function(df) {

  # Sort by alignment start
  df <- df %>%
    arrange(alifrom, alito)

  n <- nrow(df)

  if (n == 1) {
    return(
      c(
        as.character(df$queryname[1]),
        "single_domain"
      )
    )
  }

  architecture <- character(0)

  i <- 1

  while (i <= n) {

    # Current domain
    current <- df[i, ]

    # If this is the final domain
    if (i == n) {

      architecture <- c(
        architecture,
        as.character(current$queryname)
      )

      break
    }

    # Next domain
    next_domain <- df[i + 1, ]

    # Determine whether domains overlap
    #
    # Current domain ends after next domain starts
    overlap <- current$alito >= next_domain$alifrom

    if (!overlap) {

      # No overlap
      architecture <- c(
        architecture,
        as.character(current$queryname),
        "not_over"
      )

      i <- i + 1

    } else {

      # Overlap: retain domain with higher score
      if (current$score >= next_domain$score) {

        architecture <- c(
          architecture,
          as.character(current$queryname),
          "overlapped"
        )

      } else {

        architecture <- c(
          architecture,
          as.character(next_domain$queryname),
          "overlapped"
        )
      }

      # Skip both overlapping domains
      i <- i + 2
    }
  }

  architecture
}


# ------------------------------------------------------------
# Process each RIP
# ------------------------------------------------------------

RIP <- unique(tib$RIP)

l1 <- vector("list", length(RIP))
names(l1) <- RIP

multi_dom_rip <- character(0)

for (i in seq_along(RIP)) {

  ri <- RIP[i]

  tab1 <- tib %>%
    filter(RIP == ri) %>%
    arrange(alifrom, alito)

  n_domains <- nrow(tab1)

  # ----------------------------------------------------------
  # Single domain
  # ----------------------------------------------------------

  if (n_domains == 1) {

    next
  }

  # ----------------------------------------------------------
  # Multiple domains
  # ----------------------------------------------------------

  multi_dom_rip <- c(
    multi_dom_rip,
    as.character(ri)
  )

  l1[[ri]] <- resolve_architecture(tab1)
}


# ------------------------------------------------------------
# Remove empty entries
# ------------------------------------------------------------

l1 <- l1[names(l1) %in% multi_dom_rip]


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

saveRDS(
  l1,
  "domain_architecture.RDS"
)

write_delim(
  tibble(RIP = multi_dom_rip),
  "multi_dom_rip.tsv",
  delim = "\t"
)

write_delim(
  single_dom,
  "single_dom_rip.tsv",
  delim = "\t"
)