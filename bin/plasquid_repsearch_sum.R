#!/usr/bin/env Rscript
#
# Usage: plasquid_repsearch_sum.R <plasmid_report.tsv> <assembly.fna>
#
# Final REPSEARCH-chain step: pulls the accepted plasmid contigs' own
# nucleotide sequences out of the assembly, in plasmid_report.tsv's order,
# and renames the report's columns to their published names.
#
# Output: result.tsv, result.fasta

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: plasquid_repsearch_sum.R <plasmid_report.tsv> <assembly.fna>")
}
plasmid_report_file <- args[1]
assembly_file        <- args[2]

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(readr))
.this_dir <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else "."
})
source(file.path(.this_dir, "plasquid_utils.R"))

assembly <- readDNAStringSet(assembly_file)
names(assembly) <- fasta_id(names(assembly))
contigs <- orf_to_contig(names(assembly)) # sequence id -> its parent contig id

report <- read_delim(plasmid_report_file, delim = "\t", col_types = cols(contig_length = "c"))

# Guard against a Contig with no matching sequence in the assembly (should not
# happen in normal operation, since both come from the same representative set,
# but a plain match()-based index would otherwise hard-crash the whole run on an
# NA subscript instead of just dropping the offending row).
found <- report$Contig %in% contigs
if (any(!found)) {
  warning(
    "Dropping ", sum(!found), " row(s) with no matching sequence in ", assembly_file, ": ",
    paste(unique(report$Contig[!found]), collapse = ", ")
  )
  report <- report[found, ]
}

plasmid_seqs <- assembly[match(report$Contig, contigs)]
names(plasmid_seqs) <- report$Contig
writeXStringSet(plasmid_seqs, "result.fasta")

colnames(report) <- c("Contig", "RIP_domain", "MOB_group", "Rep_type", "contig_length")
write_delim(report, "result.tsv", delim = "\t")
