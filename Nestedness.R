# Nestedness

## Preliminary analysis
require(vegan)
require(phyloseq)
require(tidyverse)
require(ggh4x)
require(patchwork)
ASV.tax=readRDS("bro5.otus.rds")
ASV.tax@sam_data[,c("samp.troph", "bro_number")] -> info.TropBro

# Coerce the object into dataframe object and order the rows by trophic compartment.
data.frame(info.TropBro) -> info.TropBro
arrange(info.TropBro, match(samp.troph, c("Leaf litter", "Toxorhynchites", "Detritivore", "Water"))) -> order.TropBro
rarefy_even_depth(ASV.tax, rngseed = TRUE, sample.size = min(sample_sums(ASV.tax)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE) -> ASV.rare
# Generate the nestedness values for taxonomy based on the 16S data and make a quick plot.
nest.asv=nestednodf(ASV.rare@otu_table)
# N columns  : 28.26645 
# N rows     : 33.6054 
# NODF       : 28.26753 
# Matrix fill: 0.2021898 
saveRDS(nest.asv, "nest.asv")
#calculate genus level nestedness
tax_glom(ASV.tax, taxrank = "Genus") -> ASVgenus
rarefy_even_depth(ASVgenus, rngseed = TRUE, sample.size = min(sample_sums(ASVgenus)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE) -> genus.rare
nestednodf(genus.rare@otu_table) -> nest.genus
# N columns  : 49.05493 
# N rows     : 54.95536 
# NODF       : 49.06238 
# Matrix fill: 0.3225865
saveRDS(nest.genus, "nest.genus")
#calculate family level nestedness
tax_glom(ASV.tax, taxrank = "Family") -> ASVfam
rarefy_even_depth(ASVfam, rngseed = TRUE, sample.size = min(sample_sums(ASVfam)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE) -> family.rare
nestednodf(family.rare@otu_table) -> nest.family
nest.family
# N columns  : 57.94213 
# N rows     : 64.58097 
# NODF       : 57.96977 
# Matrix fill: 0.3811258 
#calculate order level nestedness
tax_glom(ASV.tax, taxrank = "Order") -> ASVorder
rarefy_even_depth(ASVorder, rngseed = TRUE, sample.size = min(sample_sums(ASVorder)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE) -> order.rare
nestednodf(order.rare@otu_table) -> nest.order
nest.order
# N columns  : 64.46042 
# N rows     : 71.27819 
# NODF       : 64.54739 
# Matrix fill: 0.4319767
#calculate class level nestedness
tax_glom(ASV.tax, taxrank = "Class") -> ASVclass
rarefy_even_depth(ASVclass, rngseed = TRUE, sample.size = min(sample_sums(ASVclass)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE) -> class.rare
nestednodf(class.rare@otu_table) -> nest.class
nest.class
# N columns  : 70.75812 
# N rows     : 79.66508 
# NODF       : 71.3435 
# Matrix fill: 0.4560811
#calculate phylum level nestedness
tax_glom(ASV.tax, taxrank = "Phylum") -> ASVphylum
rarefy_even_depth(ASVphylum, rngseed = TRUE, sample.size = min(sample_sums(ASVphylum)), replace = TRUE, trimOTUs = TRUE, verbose = TRUE) -> phylum.rare
nestednodf(phylum.rare@otu_table) -> nest.phylum
nest.phylum
# N columns  : 79.16415 
# N rows     : 84.64756 
# NODF       : 80.91222 
# Matrix fill: 0.5982759 
# Extract the community dataframe from the taxonomic nestedness output above.
as.data.frame(nest.genus$comm) -> nesttax.df

# Reorder the rows in the extracted sample data to match that of the nested pattern.
info.TropBro[rownames(nesttax.df),] -> taxinfo.TropBro

# Append the sample data to the nestedness dataframe.
nesttax.df$TrophicLevel <- taxinfo.TropBro$samp.troph

# Create a new column that matches the row names such that the dataframe can be melted according to two ID variables.
nesttax.df$SampNames <- row.names(nesttax.df)

# Melt the dataframe such that a custom ggplot can be successfuly made.
reshape2::melt(nesttax.df, id.vars = c("SampNames", "TrophicLevel")) -> nesttax.melt

# Define colors for facet wrapping.
stripcolors <- strip_themed(background_y = elem_list_rect(fill = c("#5f7fb9", "#426435", "#e1c751", "#161616")), text_y=elem_list_text(colour = "white"))

# Make a nice figure from the 16S taxonomic nestedness output.
nested_tax=ggplot(data = nesttax.melt, mapping = aes(x = variable, y = factor(SampNames), fill = value)) +
  geom_tile() +
  scale_fill_gradient(high = "gray35", low = "white") +
  facet_wrap2(~factor(TrophicLevel,levels = c("Water", "Leaf litter", "Detritivore","Toxorhynchites"), labels = c("Water", "Leaf", "Detr", "Tox")), ncol = 1, scales = "free_y", strip.position = "left", strip = stripcolors) +
  labs(x = "Presence or absence of taxon", y = "Trophic level") +
  theme(axis.title = element_blank(),  # Remove axis titles
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.05, "lines"))
# Double check for colorblind accessibility.
cvdPlot()


## Functional nestedness
KEGG.bac.abund=readRDS("FunPS.rds")
#subsample to minimum depth, here 458562
#892890, and 1463994 are next lowest
min(sample_sums(KEGG.bac.abund))
rar.kegg.abund=rarefy_even_depth(KEGG.bac.abund, min(sample_sums(KEGG.bac.abund)))
sample_sums(rar.kegg.abund)
# Generate the nestedness value for function based on the KEGG output and make a quick plot.
nestednodf(rar.kegg.abund@otu_table) -> nest.KEGG; plot(nest.KEGG, names = TRUE, col = "thistle3", main = "Quick nestedness plot for function", xaxt = 'n')
# N columns  : 61.81421 
# N rows     : 84.89705 
# NODF       : 61.81429 
# Matrix fill: 0.4260923 
saveRDS(nest.KEGG, "nest.KEGG")
as.data.frame(nest.KEGG$comm) -> nestfn.df

# Reorder the rows in the extracted sample data to match that of the nested pattern. First, assign the same row names to the same matrices.
order.TropBro[rownames(nestfn.df),] -> fninfo.TropBro

# Append the sample data to the nestedness dataframe.
nestfn.df$TrophicLevel <- fninfo.TropBro$samp.troph

# Create a new column that matches the row names such that the dataframe can be melted according to two ID variables.
nestfn.df$SampNames <- row.names(nestfn.df)

# Melt the dataframe such that a custom ggplot can be successfuly made.
reshape2::melt(nestfn.df, id.vars = c("SampNames", "TrophicLevel")) -> nestfn.melt

# Define colors for facet wrapping.
stripcolors <- strip_themed(background_y = elem_list_rect(fill = c("#5f7fb9", "#426435", "#e1c751", "#161616")), text_y=elem_list_text(colour = "white"))

# Make a nice figure from the functional nestedness output.
nested_function=ggplot(data = nestfn.melt, mapping = aes(x = variable, y = factor(SampNames), fill = value)) +
  geom_tile() +
  scale_fill_gradient(high = "gray35", low = "white") +
  facet_wrap2(~factor(TrophicLevel,levels = c("Water", "Leaf litter", "Detritivore","Toxorhynchites"), labels = c("Water", "Leaf", "Detr", "Tox")), ncol = 1, scales = "free_y", strip.position = "left", strip = stripcolors) +
  labs(x = "Presence or absence of taxon", y = "Trophic level") +
  theme(axis.title = element_blank(),  # Remove axis titles
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0.05, "lines"))

# Double check for colorblind accessibility.
cvdPlot()

nested_taxonomy=as.data.frame(cbind(nest.asv$statistic[3], nest.genus$statistic[3],nest.family$statistic[3],nest.order$statistic[3], nest.class$statistic[3], nest.phylum$statistic[3]))
colnames(nested_taxonomy)=c("ASV","Genus","Family","Order","Class","Phylum")
nested_taxonomy <- data.frame(Metric = factor(c("ASV", "Genus", "Family", "Order", "Class", "Phylum"),
                  levels = c("ASV", "Genus", "Family", "Order", "Class", "Phylum")),
  Value = c(28.26753, 49.06238, 57.96977, 64.54739, 71.3435, 80.91222))

line_plot <- ggplot(nested_taxonomy, aes(x = Metric, y = Value, group = 1)) +
  geom_line(color = "gray35", size = 1) +
  geom_point(color = "gray35", size = 2) +
  labs(x = "", y = "NODF") +  # set x-axis label to blank
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#######
#######
#Combine
########
  # Top row: nested_tax (left) and line_plot (right) with line_plot set to 2 cm wide.
  # Assuming the total figure width is 8.9 cm, nested_tax gets ~6.6 cm.
  top_row <- nested_tax + line_plot + plot_layout(widths = c(6.5, 2.4))
  
  # Bottom row: nested_function
  combined_plot <- top_row / nested_function + 
    plot_annotation(tag_levels = 'A')
  
  # --- 5. Export the figure ---
  # Export as a PDF 8.9 cm wide (8.9 cm ≈ 3.5 inches)
  ggsave("nestedness_figure.pdf", combined_plot,
         width = 8.9/2.54,  # convert cm to inches
         height = 5,        # adjust height as needed
         units = "in")
### Null model simulation for nestedness

#Below, I compute the null model simulations to make sure that the nestedness observed in the above section is significantly different from a rondom null community. The output is suppressed (i.e., "eval = FALSE") for this markdown because each line takes approximately 5-6 hours to run. You can run these from the original code.

# To run the simulation as a background job, remember to load the "vegan" package first. Uncomment the library first.
# library(vegan)

# Run null model simulations for the taxonomic data.
oecosimu(comm = clean.tax.df, nestfun = nestednodf, method = "r0", nsimul = 999) -> null.tax; null.tax

# Run null model simulations for the functional data.
oecosimu(comm = KEGG.bac, nestfun = nestednodf, method = "r0", nsimul = 999) -> null.fn; null.fn

```

## Analysis by bromeliad block

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
