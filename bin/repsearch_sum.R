#!/usr/bin/env Rscript

library(Biostrings)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

tbb = args[1]
dnn = args[2]

dna <- readDNAStringSet(dnn)

nms <- sub(" .*", "", names(dna))
ctg <- sub("_.*", "", nms)
names(dna) <- nms

tab <- read_delim(tbb, delim = "\t", col_types = cols(contig_length = "c"))

cnt <- tab$Contig
idx <- match(cnt, ctg)

dnr <- dna[idx]
names(dnr)<-tab$Contig
writeXStringSet(dnr, "Result.fasta")

colnames(tab) <- c("Contig", "RIP_domain", "MOB_group", "Rep_type", "contig_length")
write_delim(tab, "Result.tsv", delim = "\t")
