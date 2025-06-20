# Load vegan for Mantel tests
require(vegan)
library(phyloseq)

## ---- LOAD DATA ----

FunPS <- readRDS("FunPSabund.rds")      # Functional KEGG abundance phyloseq object
ASV.tax <- readRDS("bro5.otus.rds")     # ASV-level taxonomic phyloseq object

## ---- GENERATE HIGHER TAXONOMIC LEVEL PHYLOSEQ OBJECTS ----

# Aggregate ASVs to higher taxonomic levels using tax_glom()
ASVgenus  <- tax_glom(ASV.tax, taxrank = "Genus")
ASVfam    <- tax_glom(ASV.tax, taxrank = "Family")
ASVorder  <- tax_glom(ASV.tax, taxrank = "Order")
ASVphylum <- tax_glom(ASV.tax, taxrank = "Phylum")

## ---- CREATE BRAY-CURTIS DISTANCE MATRICES ----

# Functional (KEGG) distance matrix
kegg.braydist <- vegdist(decostand(otu_table(FunPS), "hellinger"), method = "bray")

# Taxonomic distance matrices at different levels
braydist.asv    <- vegdist(decostand(otu_table(ASV.tax), "hellinger"), method = "bray")
braydist.genus  <- vegdist(decostand(otu_table(ASVgenus), "hellinger"), method = "bray")
braydist.fam    <- vegdist(decostand(otu_table(ASVfam), "hellinger"), method = "bray")
braydist.order  <- vegdist(decostand(otu_table(ASVorder), "hellinger"), method = "bray")
braydist.phylum <- vegdist(decostand(otu_table(ASVphylum), "hellinger"), method = "bray")

## ---- RUN MANTEL TESTS ----

mant.asv    <- mantel(kegg.braydist, braydist.asv, method = "pearson"); mant.asv
mant.genus  <- mantel(kegg.braydist, braydist.genus, method = "pearson"); mant.genus
mant.fam    <- mantel(kegg.braydist, braydist.fam, method = "pearson"); mant.fam
mant.order  <- mantel(kegg.braydist, braydist.order, method = "pearson"); mant.order
mant.phylum <- mantel(kegg.braydist, braydist.phylum, method = "pearson"); mant.phylum

## ---- PLOT FUNCTIONAL VS ASV DISTANCE ----

# Convert to vectors for scatterplot
Taxvector <- as.vector(braydist.asv)
Keggvector <- as.vector(kegg.braydist)

# Plot to PDF
pdf("distance_matrix_correlation.pdf", width = 7, height = 7)
par(mar = c(5.1, 5.5, 4.1, 2.1))  # Margins: bottom, left, top, right

plot(Taxvector, Keggvector,
     xlab = "Taxonomic Distance",
     ylab = "Functional Distance",
     pch = 16,
     col = rgb(0, 0, 1, 0.5),
     cex.lab = 2)

cor_test <- cor.test(Taxvector, Keggvector, method = "pearson")

legend("topright",
       legend = c(
         paste("r =", round(cor_test$estimate, 3)),
