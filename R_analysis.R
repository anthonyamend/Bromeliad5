## ---- INSTALL & LOAD REQUIRED PACKAGES ----

# Install Bioconductor manager (for Bioconductor packages)
install.packages("BiocManager", repos = "https://cloud.r-project.org")
require(BiocManager)

# Install specific Bioconductor packages if not already installed
BiocManager::install(c("microbiome", "ComplexHeatmap"), update = FALSE)

# Install microViz package from a personal R-universe repository
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)

# Load all required libraries for analysis
# These libraries support sequence analysis, statistics, plotting, and microbiome workflows
library(BiocManager)
library(Biostrings)
library(car)
library(doBy)
library(dplyr)
library(foreach)
library(ggplot2)
library(ggpubr)
library(here)
library(iNEXT)
library(microseq)
library(multcompView)
library(pheatmap)
library(phyloseq)
library(plotrix)
library(readxl)
library(reshape2)
library(scales)
library(smartsnp)
library(SQMtools)
library(tidyverse)
library(trekcolors)
library(vegan)

# Optional: unload a package (example shown commented out)
# detach(package:fossil)

# Confirm working directory using 'here' package
here()

# Optional: list available color palettes in trekcolors
# lcars_colors()

# Use `grateful` (if installed) to generate citation-friendly list of packages
scan_packages(pkgs = "All")
cite_packages(output = "paragraph", out.dir = ".")

## ---- LOAD SQUEEZEMETA DATA (commented) ----

# Example: Load and subset SqueezeMeta output by domain "Bacteria"
# These lines are commented out but illustrate the process for loading and subsetting:
# setwd("~/Dropbox/.../Squeezemeta_outputs")
# LYBR1D = loadSQM("LYBR1D.zip")
# LYBR1D = subsetTax(LYBR1D, "superkingdom", "Bacteria")
# ...
# BacOnly = combineSQMlite(LYBR1D, LYBR1Detrit, ..., LYBR5W1)
# rm(LYBR1D, LYBR1Detrit, ...)  # Clean up workspace
# saveRDS(BacOnly, "BacOnly")

## ---- LOAD METADATA ----

# Read in the sample metadata CSV (e.g. with trophic level or treatment info)
BRO5_metadata <- read.csv("~/Dropbox/.../BRO5_metadata.csv", row.names = 1)

## ---- LOAD & PROCESS FUNCTIONAL DATA ----

# Load saved RDS object containing combined bacterial SqueezeMeta output
BacOnly <- readRDS("~/Dropbox/.../BacOnly")

# Convert functional KEGG abundances to a data frame (transpose from genes × samples to samples × genes)
KEGG.bac <- as.data.frame(t(BacOnly$functions$KEGG$abund))

# Fix sample name inconsistency (sample 8)
rn <- BacOnly$misc$samples
rn[[8]] <- "LYBR2W2"
rownames(KEGG.bac) <- rn

# Read KO name annotations and convert to tax_table format
FUNKEGGNAMES <- as.matrix(read.csv("keggnames.csv", row.names = 1))
Funtax <- tax_table(FUNKEGGNAMES)

# Merge into a phyloseq object: OTU table (KEGG), sample metadata, and KO tax table
FunPS <- merge_phyloseq(
  otu_table(KEGG.bac, taxa_are_rows = FALSE),
  sample_data(BRO5_metadata),
  Funtax
)

# Save phyloseq object and KEGG name file
write.csv(BacOnly$misc$KEGG_names, "KEGGnames.CSV")
saveRDS(FunPS, "FunPS.rds")

## ---- ALTERNATE VERSION WITH EXPLICIT KO NAMES ----

# Create a version with KO as taxa names (useful for abundance visualization)
KEGG.bac.abund <- as.data.frame(t(BacOnly$functions$KEGG$abund))
rownames(KEGG.bac.abund) <- rn

# Add KO IDs as an explicit column in taxonomy
FUNKEGGNAMES <- as.matrix(read.csv("keggnames.csv", row.names = 1))
FUNKEGGNAMES <- cbind(KO = rownames(FUNKEGGNAMES), FUNKEGGNAMES)
colnames(FUNKEGGNAMES) <- c("KO", "KOnames")
Funtax <- tax_table(FUNKEGGNAMES)

# Merge into phyloseq object
FunPSabund <- merge_phyloseq(
  otu_table(KEGG.bac.abund, taxa_are_rows = FALSE),
  sample_data(BRO5_metadata),
  Funtax
)

saveRDS(FunPSabund, "FunPSabund.rds")

## ---- TPM-NORMALIZED VERSION ----

# Same structure but with TPM-normalized values
KEGG.bac.tpm <- as.data.frame(t(BacOnly$functions$KEGG$tpm))
rownames(KEGG.bac.tpm) <- rn

# Use previously defined KO names and create phyloseq object
FunPStps <- merge_phyloseq(
  otu_table(KEGG.bac.tpm, taxa_are_rows = FALSE),
  sample_data(BRO5_metadata),
  Funtax
)

saveRDS(FunPStps, "FunPSabund.rds")
write.csv(BacOnly$misc$KEGG_names, "KEGGnames.CSV")

## ---- READ OTU PHYLOSEQ OBJECT ----

# Load phyloseq object representing OTU table (from DADA2 or other pipeline)
ASV.tax <- readRDS("bro5.otus.rds")

# Optionally re-save to a known location
# saveRDS(ASV.tax, "bro5.otus.rds")
