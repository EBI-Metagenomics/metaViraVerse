#!/usr/bin/env Rscript
install.packages(c("data.table", "stringr", "networkD3", "htmlwidgets"))

# Load required packages
suppressMessages({
  library(data.table)
  library(stringr)
  library(networkD3)
})

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)

input_file <- NULL
for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--input") {
    input_file <- args[i + 1]
  }
}

if (is.null(input_file)) {
  stop("❌ Error: Missing required argument --input\nUsage: Rscript taxonomy_sankey.R --input <gff_file>")
}

if (!file.exists(input_file)) {
  stop(paste("❌ File not found:", input_file))
}

# ---- Read GFF ----
cat("📖 Reading GFF file...")
gff <- fread(input_file, sep = "\t", header = FALSE, comment.char = "#", quote = "")
if (ncol(gff) < 9) stop("❌ GFF file does not have 9 columns.")

# ---- Extract taxonomy strings ----
cat("🔍 Extracting taxonomy info...")
taxonomy_raw <- str_extract(gff$V9, "taxonomy=[^;]+")
taxonomy_clean <- str_replace(taxonomy_raw, "taxonomy=", "")
taxonomy_clean <- taxonomy_clean[!is.na(taxonomy_clean)]

if (length(taxonomy_clean) == 0) stop("❌ No taxonomy= entries found in the GFF file.")

# ---- Split taxonomy into levels ----
taxonomy_split <- str_split(taxonomy_clean, ";\\s*")

# Build edge list for Sankey
edges <- data.table()
for (tax_path in taxonomy_split) {
  tax_path <- str_trim(tax_path)
  if (length(tax_path) > 1) {
    for (i in seq_len(length(tax_path) - 1)) {
      edges <- rbind(edges, data.table(source = tax_path[i], target = tax_path[i + 1]))
    }
  }
}

# ---- Aggregate and prepare Sankey data ----
edges_count <- edges[, .N, by = .(source, target)]
nodes <- data.table(name = unique(c(edges_count$source, edges_count$target)))

edges_count[, source_id := match(source, nodes$name) - 1]
edges_count[, target_id := match(target, nodes$name) - 1]

# ---- Build Sankey ----
cat("🎨 Generating Sankey plot...")
sankey <- sankeyNetwork(
  Links = edges_count,
  Nodes = nodes,
  Source = "source_id",
  Target = "target_id",
  Value = "N",
  NodeID = "name",
  fontSize = 12,
  nodeWidth = 25
)

htmlwidgets::saveWidget(sankey, "taxonomy_sankey.html", selfcontained = TRUE)
cat("✅ Sankey plot saved as taxonomy_sankey.html\n")
