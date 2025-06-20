
# Load required libraries
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)
library(multcompView)
library(tibble)

# Load phyloseq objects
ASV.tax <- readRDS("bro5.otus.rds")
KOASV   <- readRDS("FunPStpm.rds")

# NMDS for ASVs
nmds.tax <- metaMDS(decostand(otu_table(ASV.tax), "hellinger"))
nmds.taxscores <- as.data.frame(scores(nmds.tax)$sites)
nmds.taxscores <- cbind(nmds.taxscores,
                        Bromeliad = ASV.tax@sam_data$bro_number,
                        TrophicLevel = ASV.tax@sam_data$samp.troph)
nmds.taxscores$Bromeliad <- as.factor(nmds.taxscores$Bromeliad)
nmds.taxscores$label_color <- ifelse(as.character(nmds.taxscores$TrophicLevel) == "Toxorhynchites", "white", "black")

taxplot <- ggplot(nmds.taxscores, aes(NMDS1, NMDS2)) +
  stat_ellipse(aes(fill = TrophicLevel), geom = "polygon", alpha = 0.4, level = 0.8, color = NA) +
  geom_line(aes(group = Bromeliad), color = "black") +
  geom_point(size = 6, shape = 21, aes(fill = TrophicLevel), color = "black") +
  geom_text(aes(label = Bromeliad, color = label_color), fontface = "bold", size = 3, show.legend = FALSE) +
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) +
  scale_fill_manual(values = c("#e1c751", "#426435", "#161616", "#5f7fb9")) +
  scale_color_identity() +
  ggtitle("ASV")

# NMDS for KEGG
KOnmds <- metaMDS(decostand(otu_table(KOASV), "hellinger"))
KOnmdsscores <- as.data.frame(scores(KOnmds)$sites)
KOnmdsscores <- cbind(KOnmdsscores,
                      KOBromeliad = KOASV@sam_data$bro_number,
                      KOTrophicLevel = KOASV@sam_data$samp.troph)
KOnmdsscores$KOBromeliad <- as.factor(KOnmdsscores$KOBromeliad)
KOnmdsscores$KOlabel_color <- ifelse(as.character(KOnmdsscores$KOTrophicLevel) == "Toxorhynchites", "white", "black")

KOplot <- ggplot(KOnmdsscores, aes(NMDS1, NMDS2)) +
  stat_ellipse(aes(fill = KOTrophicLevel), geom = "polygon", alpha = 0.4, level = 0.8, color = NA) +
  geom_line(aes(group = KOBromeliad), color = "black") +
  geom_point(size = 6, shape = 21, aes(fill = KOTrophicLevel), color = "black") +
  geom_text(aes(label = KOBromeliad, color = KOlabel_color), fontface = "bold", size = 3, show.legend = FALSE) +
  theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA)) +
  scale_fill_manual(values = c("#e1c751", "#426435", "#161616", "#5f7fb9")) +
  scale_color_identity() +
  ggtitle("KEGG")

# Combine NMDS plots
combined_plot <- taxplot + KOplot +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", plot.margin = margin(2, 2, 2, 2))

ggsave("Combined_ASV_KEGG.pdf", combined_plot, width = 183, height = 74, units = "mm")

# Procrustes analysis
protest.out <- protest(nmds.tax, KOnmds)
summary(protest.out)
prores.df <- data.frame(
  samp.troph = KOASV@sam_data$samp.troph,
  bro_number = KOASV@sam_data$bro_number,
  protest.residual = residuals(protest.out)
)

# ANOVA and Tukey HSD on Procrustes residuals
fml <- aov(protest.residual ~ samp.troph, data = prores.df)
tukey <- TukeyHSD(fml, conf.level = 0.95)
tuk.letters <- multcompLetters4(fml, tukey)

tktbl <- prores.df %>%
  group_by(samp.troph) %>%
  summarise(resmean = mean(protest.residual), sd = sd(protest.residual)) %>%
  arrange(samp.troph)
tuk.letters <- as.data.frame.list(tuk.letters$samp.troph)
tktbl$Letters <- tuk.letters$Letters
df.boxplot <- left_join(prores.df, tktbl, by = "samp.troph")

ggplot(data = df.boxplot, aes(x = samp.troph, y = protest.residual)) +
  geom_boxplot(aes(fill = samp.troph), alpha = 0.75, outlier.shape = 21) +
  scale_fill_manual(values = c("#e1c751", "#426435", "#161616", "#5f7fb9")) +
  theme(panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA)) +
  xlab("Trophic level") + ylab("Procrustes Residual")

# PERMANOVA
taxpart <- adonis2(decostand(ASV.tax@otu_table, "hellinger") ~
                   as.factor(ASV.tax@sam_data$bro_number):ASV.tax@sam_data$samp.troph +
                   as.factor(ASV.tax@sam_data$bro_number))

funpart <- adonis2(decostand(KOASV@otu_table, "hellinger") ~
                   KOASV@sam_data$samp.troph:KOASV@sam_data$bro_number +
                   KOASV@sam_data$bro_number)
