
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


ChatGPT can make mistakes. Check important info.
