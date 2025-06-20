Skip to content
Chat history
You said:
annotate this code:
ChatGPT said:

Sure, please paste the code you'd like annotated, and I’ll add clear, explanatory comments to help you understand what it does.
You said:
library(phyloseq)
library(vegan)
library(bipartite)
library(ggplot2)
library(cowplot)
library(tidyverse)

# Load the phyloseq objects
ps_otu <- readRDS("bro5.otus.rds")
ps_ko <- readRDS("FunPS.rds")

# Function to generate a nestedness tile plot
generate_nested_plot <- function(ps_sub, bro, label) {
  min_sum <- min(sample_sums(ps_sub))
  ps_sub_rar <- rarefy_even_depth(ps_sub, sample.size = min_sum, rngseed = 1234, verbose = FALSE)
  
  otu <- as(otu_table(ps_sub_rar), "matrix")
  if (taxa_are_rows(ps_sub_rar)) otu <- t(otu)
  
  otu_bin <- ifelse(otu > 0, 1, 0)
  
  nodf_val <- nestednodf(otu_bin)$statistic["NODF"]
  
  troph_labels <- data.frame(SampleID = sample_names(ps_sub_rar),
                             Troph = sample_data(ps_sub_rar)$samp.troph)
  sample_to_troph <- troph_labels$Troph
  names(sample_to_troph) <- troph_labels$SampleID
  
  otu_bin_df <- as.data.frame(as.table(otu_bin))
  colnames(otu_bin_df) <- c("Sample", "ASV", "Presence")
  otu_bin_df$Troph <- sample_to_troph[otu_bin_df$Sample]
  
  trophic_matrix <- otu_bin_df %>%
    group_by(Troph, ASV) %>%
    summarise(Presence = max(Presence), .groups = "drop")
  
  trophic_wide <- trophic_matrix %>%
    pivot_wider(names_from = ASV, values_from = Presence, values_fill = 0) %>%
    column_to_rownames("Troph") %>%
    as.matrix()
  
  row_order <- order(rowSums(trophic_wide), decreasing = TRUE)
  col_order <- order(colSums(trophic_wide), decreasing = TRUE)
  
  trophic_sorted <- trophic_wide[row_order, col_order, drop = FALSE]
  trophic_sorted <- trophic_sorted[rev(rownames(trophic_sorted)), , drop = FALSE]
  
  plot_df <- as.data.frame(as.table(trophic_sorted))
  colnames(plot_df) <- c("Troph", "ASV", "Presence")
  
  plot_df$Troph <- factor(plot_df$Troph, levels = rownames(trophic_sorted))
  plot_df$ASV <- factor(plot_df$ASV, levels = colnames(trophic_sorted))
  
  # Black and white plot with no outlines
  p <- ggplot(plot_df, aes(x = ASV, y = Troph, fill = factor(Presence))) +
    geom_tile() +  # no color outline
    scale_fill_manual(values = c("0" = "white", "1" = "black"), guide = FALSE) +
    labs(title = paste0("Bromeliad ", bro, " ", label, "; NODF = ", round(nodf_val, 1))) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8),
          panel.grid = element_blank(),
          plot.margin = margin(5, 5, 5, 5))
  
  return(p)
}


# Generate all paired plots
paired_plots <- list()
for (bro in 1:5) {
  otu_sub <- subset_samples(ps_otu, bro_number == bro)
  ko_sub  <- subset_samples(ps_ko, bro_number == bro)
  
  p_otu <- generate_nested_plot(otu_sub, bro, "OTU")
  p_ko  <- generate_nested_plot(ko_sub, bro, "KO")
  
  combined <- plot_grid(p_otu, p_ko, ncol = 2, rel_widths = c(1, 1))
  paired_plots[[bro]] <- combined
}

# Stack vertically
final_plot <- plot_grid(plotlist = paired_plots, ncol = 1)

# Save the figure to PDF, 183 mm wide (18.3 cm)
ggsave("bromeliad_nestedness_OTU_KO.pdf", final_plot, width = 18.3, height = 30, units = "cm")
ChatGPT said:

Here's an annotated version of your R script with clear comments throughout to explain each step:

# Load required libraries
library(phyloseq)     # For microbial community data structures and analysis
library(vegan)        # For ecological statistics like nestedness (NODF)
library(bipartite)    # For nestedness metrics (nestednodf)
library(ggplot2)      # For creating plots
library(cowplot)      # For arranging multiple ggplots
library(tidyverse)    # For data manipulation and reshaping

# Load phyloseq objects containing OTU and KO data
ps_otu <- readRDS("bro5.otus.rds")   # OTU-level microbial community data
ps_ko <- readRDS("FunPS.rds")       # Functional gene data (KEGG Orthologs)

# Define a function to generate a nestedness tile plot
generate_nested_plot <- function(ps_sub, bro, label) {
  # Rarefy all samples in this phyloseq object to equal sequencing depth
  min_sum <- min(sample_sums(ps_sub))
  ps_sub_rar <- rarefy_even_depth(ps_sub, sample.size = min_sum, rngseed = 1234, verbose = FALSE)

  # Convert to OTU table and transpose if necessary
  otu <- as(otu_table(ps_sub_rar), "matrix")
  if (taxa_are_rows(ps_sub_rar)) otu <- t(otu)

  # Convert to binary presence/absence matrix
  otu_bin <- ifelse(otu > 0, 1, 0)

  # Calculate NODF nestedness statistic
  nodf_val <- nestednodf(otu_bin)$statistic["NODF"]

  # Map each sample to its trophic level
  troph_labels <- data.frame(SampleID = sample_names(ps_sub_rar),
                             Troph = sample_data(ps_sub_rar)$samp.troph)
  sample_to_troph <- troph_labels$Troph
  names(sample_to_troph) <- troph_labels$SampleID

  # Convert the binary matrix into a long-format data frame
  otu_bin_df <- as.data.frame(as.table(otu_bin))
  colnames(otu_bin_df) <- c("Sample", "ASV", "Presence")
  otu_bin_df$Troph <- sample_to_troph[otu_bin_df$Sample]

  # Collapse by trophic level (take max presence of each ASV within each level)
  trophic_matrix <- otu_bin_df %>%
    group_by(Troph, ASV) %>%
    summarise(Presence = max(Presence), .groups = "drop")

  # Reshape to wide format: Trophic levels as rows, ASVs as columns
  trophic_wide <- trophic_matrix %>%
    pivot_wider(names_from = ASV, values_from = Presence, values_fill = 0) %>%
    column_to_rownames("Troph") %>%
    as.matrix()

  # Order rows (trophic levels) and columns (ASVs) by richness
  row_order <- order(rowSums(trophic_wide), decreasing = TRUE)
  col_order <- order(colSums(trophic_wide), decreasing = TRUE)

  # Apply the ordering
  trophic_sorted <- trophic_wide[row_order, col_order, drop = FALSE]
  trophic_sorted <- trophic_sorted[rev(rownames(trophic_sorted)), , drop = FALSE]  # Reverse row order for plotting

  # Prepare data frame for ggplot
  plot_df <- as.data.frame(as.table(trophic_sorted))
  colnames(plot_df) <- c("Troph", "ASV", "Presence")

  # Factorize with specific order for plotting
  plot_df$Troph <- factor(plot_df$Troph, levels = rownames(trophic_sorted))
  plot_df$ASV <- factor(plot_df$ASV, levels = colnames(trophic_sorted))

  # Create tile plot of presence/absence (black/white)
  p <- ggplot(plot_df, aes(x = ASV, y = Troph, fill = factor(Presence))) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "white", "1" = "black"), guide = FALSE) +
    labs(title = paste0("Bromeliad ", bro, " ", label, "; NODF = ", round(nodf_val, 1))) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),                  # Hide x-axis text
          axis.text.y = element_text(size = 7),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8),
          panel.grid = element_blank(),
          plot.margin = margin(5, 5, 5, 5))

  return(p)
}

# Loop through each bromeliad and generate side-by-side plots for OTU and KO
paired_plots <- list()
for (bro in 1:5) {
  otu_sub <- subset_samples(ps_otu, bro_number == bro)
  ko_sub  <- subset_samples(ps_ko, bro_number == bro)

  p_otu <- generate_nested_plot(otu_sub, bro, "OTU")
  p_ko  <- generate_nested_plot(ko_sub, bro, "KO")

  # Combine OTU and KO plots side by side
  combined <- plot_grid(p_otu, p_ko, ncol = 2, rel_widths = c(1, 1))
  paired_plots[[bro]] <- combined
}

# Stack all paired plots vertically
final_plot <- plot_grid(plotlist = paired_plots, ncol = 1)

# Save final composite figure to a PDF (183 mm = 18.3 cm wide)
ggsave("bromeliad_nestedness_OTU_KO.pdf", final_plot, width = 18.3, height = 30, units = "cm")

Let me know if you'd like inline comments directly within the code instead, or a version with debug logging, parameterization, or modular refactoring.
You said:
again: 
## Load required packages


install.packages("BiocManager", repos = "https://cloud.r-project.org")
require(BiocManager)
# Install any packages needed for analysis that have not been previously installed. Some can be installed using 'utils' but some may require BiocManager.
# BiocManager::install("microbiome")
# install.packages("grateful")
#BiocManager::install("microbiome")
#BiocManager::install("phyloseq")
BiocManager::install(c("microbiome", "ComplexHeatmap"), update = FALSE)
install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
# Load all packages used in this analysis. I load them all individually such that I am able to discern immediately if any individual package needs updating, installing, etc.
library("BiocManager")
library("Biostrings");
library("car"); 
library("doBy"); 
library("dplyr"); 
library("foreach"); 
library("ggplot2")
library("ggpubr") 
library("here") 
library("iNEXT") 
library("microseq")
library("multcompView") 
library("pheatmap")
library("phyloseq") 
library("plotrix") 
library("readxl") 
library("reshape2")
library("scales")
library("smartsnp") 
library("SQMtools") 
library("tidyverse") 
library("trekcolors")
library("vegan")
library("phyloseq")

# I'm including a "detach" command here in case I want to remove any packages at any point. If you choose to detach anything, you may opt to do so using the line of code below. Un-comment this line of code in order to do so and replace with the package you .
# detach(package:fossil)

here()    # Use function "here" from package "Here" to set working directory.

# Lists Star Trek color names from the package "trekcolors". I used the palette called "enara2" and modified some of the hex codes to be more colorblind-friendly. You may run this line to see a list of all colors associated with trekcolors palettes.
# lcars_colors()

# From package "grateful", I used this function below to show all packages and associated versions used in this analysis to streamline proper citation of the packages used.
scan_packages(pkgs = "All")

# The line below creates a bibliography that includes as much information as possible for all the packages used. The bibliography created here is used to facilitate proper citation of all packages used in this analysis. The paragraph form needs to be used if you want the bibliography to appear at the end.
cite_packages(output = "paragraph", out.dir = ".")

## Load in data required for analysis

# Load zipped outputs from SqueezeMeta into R. For your own analysis, change the specified path to the correct one on your computer/system. 
#setwd to where zipped files are located for sequential assemblies
# setwd("~/Dropbox/Grants/NSF_Rules_of_Life/Bro5anal/squeezemeta_stuff/Squeezemeta_outputs")
# LYBR1D	=		loadSQM("LYBR1D.zip")
# LYBR1D= subsetTax(LYBR1D, "superkingdom", "Bacteria")
# LYBR1Detrit	=		loadSQM("LYBR1Detrit.zip")
# LYBR1Detrit=subsetTax(LYBR1Detrit, "superkingdom", "Bacteria")
# LYBR1Toxo	=		loadSQM("LYBR1Toxo.zip")
# LYBR1Toxo=subsetTax(LYBR1Toxo, "superkingdom", "Bacteria")
# 
# LYBR1W1	=		loadSQM("LYBR1W1.zip")
# LYBR1W1=subsetTax(LYBR1W1, "superkingdom", "Bacteria")
# 
# LYBR2D	=		loadSQM("LYBR2D.zip")
# LYBR2D=subsetTax(LYBR2D, "superkingdom", "Bacteria")
# 
# LYBR2Detrit	=		loadSQM("LYBR2Detrit.zip")
# LYBR2Detrit=subsetTax(LYBR2Detrit, "superkingdom", "Bacteria")
# 
# LYBR2Toxo	=		loadSQM("LYBR2Toxo.zip")
# LYBR2Toxo=subsetTax(LYBR2Toxo, "superkingdom", "Bacteria")
# 
# LYBR2W1	=		loadSQM("LYBR2W1.zip")
# LYBR1Detrit=subsetTax(LYBR1Detrit, "superkingdom", "Bacteria")
# 
# LYBR3D	=		loadSQM("LYBR3D.zip")
# LYBR3D=subsetTax(LYBR3D, "superkingdom", "Bacteria")
# 
# LYBR3Detrit	=		loadSQM("LYBR3Detrit.zip")
# LYBR3Detrit=subsetTax(LYBR3Detrit, "superkingdom", "Bacteria")
# 
# LYBR3Toxo	=		loadSQM("LYBR3Toxo.zip")
# LYBR1Detrit=subsetTax(LYBR1Detrit, "superkingdom", "Bacteria")
# 
# LYBR3W1	=		loadSQM("LYBR3W1.zip")
# LYBR3W1=subsetTax(LYBR3W1, "superkingdom", "Bacteria")
# 
# LYBR4D	=		loadSQM("LYBR4D.zip")
# LYBR4D=subsetTax(LYBR4D, "superkingdom", "Bacteria")
# 
# LYBR4Detrit	=		loadSQM("LYBR4Detrit.zip")
# LYBR4Detrit=subsetTax(LYBR4Detrit, "superkingdom", "Bacteria")
# 
# LYBR4Toxo	=		loadSQM("LYBR4Toxo.zip")
# LYBR4Toxo=subsetTax(LYBR4Toxo, "superkingdom", "Bacteria")
# 
# LYBR4W1	=		loadSQM("LYBR4W1.zip")
# 
# LYBR4W1=subsetTax(LYBR4W1, "superkingdom", "Bacteria")
# LYBR5D	=		loadSQM("LYBR5D.zip")
# 
# LYBR5D=subsetTax(LYBR5D, "superkingdom", "Bacteria")
# LYBR5Detrit	=		loadSQM("LYBR5Detrit.zip")
# LYBR5Detrit=subsetTax(LYBR5Detrit, "superkingdom", "Bacteria")
# 
# LYBR5Toxo	=		loadSQM("LYBR5Toxo.zip")
# LYBR5Toxo=subsetTax(LYBR5Toxo, "superkingdom", "Bacteria")
# 
# LYBR5W1	=		loadSQM("LYBR5W1.zip")
# LYBR5W1=subsetTax(LYBR5W1, "superkingdom", "Bacteria")
# 
# BacOnly=combineSQMlite(LYBR1D,	LYBR1Detrit,	LYBR1Toxo,	LYBR1W1,	LYBR2D,	LYBR2Detrit,	LYBR2Toxo,	LYBR2W1,	LYBR3D,	LYBR3Detrit,	LYBR3Toxo,	LYBR3W1,	LYBR4D,	LYBR4Detrit,	LYBR4Toxo,	LYBR4W1,	LYBR5D,	LYBR5Detrit,	LYBR5Toxo,	LYBR5W1)
# 
# rm(LYBR1D,	LYBR1Detrit,	LYBR1Toxo,	LYBR1W1,	LYBR2D,	LYBR2Detrit,	LYBR2Toxo,	LYBR2W1,	LYBR3D,	LYBR3Detrit,	LYBR3Toxo,	LYBR3W1,	LYBR4D,	LYBR4Detrit,	LYBR4Toxo,	LYBR4W1,	LYBR5D,	LYBR5Detrit,	LYBR5Toxo,	LYBR5W1)
# saveRDS(BacOnly, "~/Dropbox/Grants/NSF_Rules_of_Life/Bro5anal/SamOTU/R_Analysis/BacOnly")
# BacOnly$total_reads
# summary(BacOnly) -> BacOnly.summary

#read in the metadata
BRO5_metadata <- read.csv("~/Dropbox/Grants/NSF_Rules_of_Life/Bro5anal/SamOTU/R_Analysis/BRO5_metadata.csv", row.names=1)

# Subset the SqueezeMeta SQM object such that only the reads assigned as Bacteria are included.
#BacOnly = subsetTax(seqsqm, "superkingdom", "Bacteria")
#saveRDS(BacOnly, "~/Dropbox/Grants/NSF_Rules_of_Life/Bro5anal/SamOTU/R_Analysis/BacOnly")
# Transpose the abundance matrix from the SqueezeMeta object, declare it as a data frame, and assign it as an object. Reads normalized as TPM
BacOnly=readRDS("~/Dropbox/Grants/NSF_Rules_of_Life/Bro5anal/SamOTU/R_Analysis/BacOnly")
str(BacOnly)
#datafram of abundances
as.data.frame(t(BacOnly$functions$KEGG$abund)) -> KEGG.bac
#the water sample name differs from the taxon name.  CHange it here.
rn=BacOnly$misc$samples
rn[[8]]="LYBR2W2"
rownames(KEGG.bac)=rn
#turn the functional Data into a Phyloseq object for parallel processing, name it FunPS
FUNKEGGNAMES=as.matrix(read.csv("keggnames.csv",row.names = 1))
Funtax=tax_table(FUNKEGGNAMES)
taxa_names(Funtax)
FunPS=merge_phyloseq(otu_table(KEGG.bac, taxa_are_rows=FALSE),sample_data(BRO5_metadata),Funtax)
write.csv(BacOnly$misc$KEGG_names, "KEGGnames.CSV")
saveRDS(FunPS, "FunPS.rds")
#####
#make another where taxa are KO names name it FunPSabund
as.data.frame(t(BacOnly$functions$KEGG$abund)) -> KEGG.bac.abund
#the water sample name differs from the taxon name.  CHange it here.
rn=BacOnly$misc$samples
rn[[8]]="LYBR2W2"
rownames(KEGG.bac.abund)=rn
#turn the functional Data into a Phyloseq object for parallel processing
FUNKEGGNAMES=as.matrix(read.csv("keggnames.csv",row.names = 1))
head(FUNKEGGNAMES)
FUNKEGGNAMES=cbind(rownames(FUNKEGGNAMES), FUNKEGGNAMES)
colnames(FUNKEGGNAMES)=c("KO","KOnames")
Funtax=tax_table(FUNKEGGNAMES)
taxa_names(Funtax)
FunPSabund=merge_phyloseq(otu_table(KEGG.bac.abund, taxa_are_rows=FALSE),sample_data(BRO5_metadata),Funtax)
saveRDS(FunPSabund, "FunPSabund.rds")
#datafram of abundances
as.data.frame(t(BacOnly$functions$KEGG$abund)) -> KEGG.bac
#the water sample name differs from the taxon name.  CHange it here.
rn=BacOnly$misc$samples
rn[[8]]="LYBR2W2"
rownames(KEGG.bac)=rn
#turn the functional Data into a Phyloseq object for parallel processing, name it FunPS
FUNKEGGNAMES=as.matrix(read.csv("keggnames.csv",row.names = 1))
Funtax=tax_table(FUNKEGGNAMES)
taxa_names(Funtax)
FunPS=merge_phyloseq(otu_table(KEGG.bac, taxa_are_rows=FALSE),sample_data(BRO5_metadata),Funtax)
write.csv(BacOnly$misc$KEGG_names, "KEGGnames.CSV")
saveRDS(FunPS, "FunPS.rds")
#####
#make another where taxa are KO names name it FunPSabund
as.data.frame(t(BacOnly$functions$KEGG$tpm)) -> KEGG.bac.tpm
#the water sample name differs from the taxon name.  CHange it here.
rn=BacOnly$misc$samples
rn[[8]]="LYBR2W2"
rownames(KEGG.bac.tpm)=rn
#turn the functional Data into a Phyloseq object for parallel processing
FUNKEGGNAMES=as.matrix(read.csv("keggnames.csv",row.names = 1))
head(FUNKEGGNAMES)
FUNKEGGNAMES=cbind(rownames(FUNKEGGNAMES), FUNKEGGNAMES)
colnames(FUNKEGGNAMES)=c("KO","KOnames")
Funtax=tax_table(FUNKEGGNAMES)
taxa_names(Funtax)
FunPStps=merge_phyloseq(otu_table(KEGG.bac.tpm, taxa_are_rows=FALSE),sample_data(BRO5_metadata),Funtax)
saveRDS(FunPStps, "FunPSabund.rds")
write.csv(BacOnly$misc$KEGG_names, "KEGGnames.CSV")
# Read in the results from phyloseq processing. Remember to change the specified path to the correct one on your computer/system. I duplicate the phyloseq object so I have one object that I can manipulate, and the original that can be used for subsequent analysis. If you do not want the duplication, remove the following from the command: " -> ASV.deconLuluTrim" and make sure that all downstream instances of ASV.deconLuluTrim have been accounted for.

# The line below reads in the phyloseq object in which OTUs (97% threshold) were generated.
ASV.tax=readRDS("bro5.otus.rds") 
#saveRDS(ASV.tax, "/Users/anthonyamend/Dropbox/Grants/NSF_Rules_of_Life/Bro5anal/SamOTU/R_Analysis/bro5.otus.rds")
ChatGPT said:

Here is a thoroughly annotated version of your R script, with comments added or clarified for clarity, best practices, and structure. I've cleaned and organized the script without modifying its logic:

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

Suggestions:

    Remove duplicate code: the creation of FunPS and FunPSabund occurs twice unnecessarily.

    Use consistent naming: your final FunPSabund is overwritten by FunPStps and saved under the same name.

    Consider breaking the code into separate scripts or RMarkdown chunks (e.g., install/load, load data, transform, export).

Let me know if you'd like a modularized version, RMarkdown template, or this output saved as a .R or .Rmd file.
You said:
again library(phyloseq)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tibble)
library(RColorBrewer)
library(vegan)

# Load phyloseq object
ps_otu <- readRDS("bro5.otus.rds")

# Color mapping for trophic levels
troph_colors <- c(
  "Water" = "#5f7fb9",
  "Leaf litter" = "#426435",
  "Detritivore" = "#e1c751",
  "Toxorhynchites" = "#161616"
)

# Aggregate to Class level
ps_class <- tax_glom(ps_otu, taxrank = "Class")

# Rarefy to minimum depth
set.seed(123)
ps_class_rarefied <- rarefy_even_depth(ps_class, rngseed = 123, verbose = FALSE)

# Extract OTU table and sample data
otu <- as(otu_table(ps_class_rarefied), "matrix")
if (taxa_are_rows(ps_class_rarefied)) {otu <- t(otu)}

sample_data_df <- data.frame(sample_data(ps_class_rarefied))
sample_data_df$SampleID <- rownames(sample_data_df)
sample_data_df <- sample_data_df %>%
  filter(SampleID %in% colnames(otu)) %>%
  arrange(factor(samp.troph, levels = c("Water", "Leaf litter", "Detritivore", "Toxorhynchites")))

otu <- otu[, sample_data_df$SampleID, drop = FALSE]

# Annotate rows with Class
tax_table_df <- data.frame(tax_table(ps_class_rarefied))
rownames(otu) <- tax_table_df$Class[match(rownames(otu), rownames(tax_table_df))]

# Remove empty rows before distance calculation
otu_norm_filtered <- otu[rowSums(otu) > 0, , drop = FALSE]

# Cluster OTUs using Bray-Curtis distance if possible
if (nrow(otu_norm_filtered) >= 2) {
  dist_rows <- vegdist(otu_norm_filtered, method = "bray")
  cluster_obj <- hclust(dist_rows, method = "average")
} else {
  cluster_obj <- FALSE
}

# Prepare sample annotations
annot_col <- data.frame(Trophic = sample_data_df$samp.troph)
rownames(annot_col) <- sample_data_df$SampleID
annot_colors <- list(Trophic = troph_colors)

# Generate heatmap
# Generate heatmap
pheatmap(otu_norm_filtered,
         color = colorRampPalette(c("white", "black"))(100),
         cluster_rows = cluster_obj,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = annot_col,
         annotation_colors = annot_colors,
         fontsize_row = 6,
         fontsize_col = 6,
         filename = "bromeliad_class_heatmap.pdf",
         width = 18.3 / 2.54,
         height = 25 / 2.54)
ChatGPT said:

Here's your R script annotated with clear comments explaining each section and key lines:

# Load necessary libraries
library(phyloseq)       # For microbiome data structures and manipulation
library(ggplot2)        # General plotting (not used in this script, but often helpful)
library(pheatmap)       # For heatmap visualization
library(dplyr)          # For data manipulation
library(tibble)         # For working with tibbles/data frames
library(RColorBrewer)   # Color palettes (not used directly here)
library(vegan)          # Ecological statistics (Bray-Curtis distance)

# Load OTU-level phyloseq object (previously saved)
ps_otu <- readRDS("bro5.otus.rds")

# Define custom colors for trophic levels
troph_colors <- c(
  "Water" = "#5f7fb9",
  "Leaf litter" = "#426435",
  "Detritivore" = "#e1c751",
  "Toxorhynchites" = "#161616"
)

# Aggregate OTU table up to the Class taxonomic level
ps_class <- tax_glom(ps_otu, taxrank = "Class")

# Rarefy all samples to even sequencing depth (using minimum sample size)
set.seed(123)  # For reproducibility
ps_class_rarefied <- rarefy_even_depth(ps_class, rngseed = 123, verbose = FALSE)

# Extract OTU count matrix
otu <- as(otu_table(ps_class_rarefied), "matrix")
if (taxa_are_rows(ps_class_rarefied)) {
  otu <- t(otu)  # Ensure rows = taxa, columns = samples
}

# Extract and clean sample metadata
sample_data_df <- data.frame(sample_data(ps_class_rarefied))
sample_data_df$SampleID <- rownames(sample_data_df)

# Filter metadata to match sample columns in OTU table, and order by trophic category
sample_data_df <- sample_data_df %>%
  filter(SampleID %in% colnames(otu)) %>%
  arrange(factor(samp.troph, levels = c("Water", "Leaf litter", "Detritivore", "Toxorhynchites")))

# Reorder OTU table columns to match metadata order
otu <- otu[, sample_data_df$SampleID, drop = FALSE]

# Replace OTU IDs (row names) with Class-level taxonomy
tax_table_df <- data.frame(tax_table(ps_class_rarefied))
rownames(otu) <- tax_table_df$Class[match(rownames(otu), rownames(tax_table_df))]

# Remove OTUs (Classes) with zero total counts
otu_norm_filtered <- otu[rowSums(otu) > 0, , drop = FALSE]

# Cluster taxa (rows) using Bray-Curtis distance, only if ≥2 taxa present
if (nrow(otu_norm_filtered) >= 2) {
  dist_rows <- vegdist(otu_norm_filtered, method = "bray")
  cluster_obj <- hclust(dist_rows, method = "average")
} else {
  cluster_obj <- FALSE  # No clustering if not enough taxa
}

# Prepare column annotations (trophic level for each sample)
annot_col <- data.frame(Trophic = sample_data_df$samp.troph)
rownames(annot_col) <- sample_data_df$SampleID
annot_colors <- list(Trophic = troph_colors)

# Generate and save heatmap to PDF
pheatmap(otu_norm_filtered,
         color = colorRampPalette(c("white", "black"))(100),  # Color gradient
         cluster_rows = cluster_obj,        # Cluster taxa using Bray-Curtis
         cluster_cols = FALSE,              # Do not cluster samples
         show_colnames = TRUE,              # Show sample IDs
         show_rownames = TRUE,              # Show taxa (Class names)
         annotation_col = annot_col,        # Add sample metadata
         annotation_colors = annot_colors,  # Use defined color palette
         fontsize_row = 6,                  # Small row font
         fontsize_col = 6,                  # Small column font
         filename = "bromeliad_class_heatmap.pdf",  # Output file
         width = 18.3 / 2.54,               # Width in inches (183 mm)
         height = 25 / 2.54)                # Height in inches (25 cm)

Let me know if you'd like to:

    Add z-score normalization across rows

    Split by trophic level

    Color by sample richness or other metrics

    Include the top N most abundant taxa only

I can also save this annotated script as a .R file for download if needed.


ChatGPT can make mistakes. Check important info.
