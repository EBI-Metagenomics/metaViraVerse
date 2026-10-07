#!/usr/bin/env Rscript
#
# Usage: plasquid_inc_classification.R <inc_candidates.tsv> <rna_candidates.tsv> \
#          <protein_to_contig.tsv>
#
# Classifies plasmid incompatibility (Inc) groups from two independent
# searches: hmmsearch hits of predicted proteins against Inc-associated
# protein profiles (inc_candidates.tsv, a HMMER --domtblout), and cmsearch
# hits of whole contigs against Inc-associated RNA (Rep/Col) covariance
# models (rna_candidates.tsv, an Infernal --tblout). Each Inc family has its
# own curated bit-score cutoff.
#
# Output: classification_table.tsv (Inc_det, query_name, score, tlen, contig)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: plasquid_inc_classification.R <inc_candidates.tsv> <rna_candidates.tsv> <protein_to_contig.tsv>")
}
domtblout_file      <- args[1]
cmsearch_file       <- args[2]
protein_contig_file <- args[3]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
.this_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
})
source(file.path(.this_dir, "plasquid_utils.R"))

# ------------------------------------------------------------
# Protein-based Inc determinants (hmmsearch --domtblout)
# ------------------------------------------------------------

protein_hits <- read_hmmer_domtblout(domtblout_file) %>%
  rename(Inc_det = query_name)
protein_contig_map <- read_protein_contig_map(protein_contig_file)

protein_cutoffs <- tibble::tribble(
  ~Inc_det,             ~min_score,
  "Inc11_B",             350.2,
  "Inc11",                193.7,
  "Inc13_A",              491.0,
  "Inc13_B",              448.4,
  "Inc13_C",              530.2,
  "Inc18",                450.6,
  "Inc1",                 532.6,
  "Inc1_B",               712.9,
  "Inc4_A",               568.4,
  "Inc4B-9-10-14",        515.6,
  "Inc7_A",               489.0,
  "Inc7_B",               551.3,
  "Inc8",                 657.4,
  "IncAC",                655.7,
  "IncB-O-K-I",           573.9,
  "IncFI_RepB",           438,
  "IncFI_RepE",           517.8,
  "IncGU",                860.4,
  "IncHIA",               451.0,
  "IncHIB",               567,
  "IncHI2",               852.1,
  "IncLM",                466.2,
  "IncN",                 257.0,
  "IncP2",                392.1,
  "IncP7",                674.8,
  "IncP9",                236.4,
  "IncP",                 593.5,
  "IncQ",                 469.7,
  "IncR",                 658.9,
  "IncT",                 675.5,
  "IncW",                 449.4,
  "IncX",                 443.1,
  "IncZ",                 515.0
)

protein_based <- pmap_dfr(protein_cutoffs, function(Inc_det, min_score) {
  protein_hits %>%
    filter(.data$Inc_det == .env$Inc_det, score >= .env$min_score) %>%
    transmute(Inc_det, tlen, query_name = target_name, score, contig = protein_to_contig(target_name, protein_contig_map))
})

# ------------------------------------------------------------
# RNA-based Inc determinants (cmsearch --tblout)
# ------------------------------------------------------------

cmsearch_hits <- read_whitespace_table(cmsearch_file, c(
  "query_name", "accession", "Inc_det", "accession2", "mdl",
  "mdl_from", "mdl_to", "seq_from", "seq_to", "strand", "trunc",
  "pass", "gc", "bias", "score", "E-value", "inc", "DoT"
)) %>%
  mutate(across(c(mdl_from, mdl_to, seq_from, seq_to, pass, gc, bias, score, `E-value`), as.numeric))

rna_cutoffs <- tibble::tribble(
  ~Inc_det,                              ~min_score, ~mdl_from_must_be_1,
  "Col156",                               138.2,      FALSE,
  "Col3M",                                106.9,      FALSE,
  "Col440I",                              115.4,      TRUE,
  "Col440II",                             248.8,      FALSE,
  "Col8282",                              193.6,      FALSE,
  "Col(BS512)",                           245.0,      FALSE,
  "ColE10",                               179.6,      FALSE,
  "Col(IMGS31)",                          197.2,      FALSE,
  "Col(IRGK)",                            190.3,      FALSE,
  "ColKP3",                               298.4,      FALSE,
  "Col(MP18)",                            182.9,      FALSE,
  "ColpVC",                               178.0,      FALSE,
  "ColRNAI",                              113.2,      FALSE,
  "Col(SD853)",                           173,        FALSE,
  "Col(Ye4449)",                          188.8,      FALSE,
  "Col(MG828)",                           182,        FALSE,
  "IncFII_1_pKP91",                       210.5,      FALSE,
  "IncFII(29)_1_pUTI89",                  266.6,      FALSE,
  "IncFII(p96A)_1_p96A",                  587.6,      FALSE,
  "IncFII(pCoo)_1_pCoo",                  264,        FALSE,
  "IncFII(pCRY)_1_pCRY",                  608.5,      FALSE,
  "IncFII(pCTU2)_1_pCTU2",                636.1,      FALSE,
  "IncFII(pECLA)_1_pECLA",                786.8,      FALSE,
  "IncFII(pHN7A8)_1_pHN7A8",              262.5,      FALSE,
  "IncFII(pKPX1)",                        585.0,      FALSE,
  "IncFII(pMET)_1_pMET1",                 629.7,      FALSE,
  "IncFII(pRSB107)_1_pRSB107",            261.8,      FALSE,
  "IncFII(pSE11)_1_pSE11",                274.2,      FALSE,
  "IncFII(pYVa12790)_1_pYVa12790",        697.6,      FALSE,
  "IncFII(S)_1",                          238.4,      FALSE,
  "IncFII(Y)_1_ps",                       217.4,      FALSE,
  "IncFII_p14_Yersenia",                  205.9,      FALSE
)

rna_based <- pmap_dfr(rna_cutoffs, function(Inc_det, min_score, mdl_from_must_be_1) {
  hits <- cmsearch_hits %>% filter(.data$Inc_det == .env$Inc_det, score >= .env$min_score)
  if (mdl_from_must_be_1) hits <- hits %>% filter(mdl_from == 1)
  hits %>% transmute(Inc_det, query_name, score)
}) %>%
  mutate(tlen = NA_real_, contig = query_name)

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

classification_table <- bind_rows(rna_based, protein_based)
write_delim(classification_table, "classification_table.tsv", delim = "\t")
