## **Figure 1a-b:** Cladogram of fungal taxa contributing significantly to community composition variation (a) and abundance of top 10 genera for each of the top 3 phyla across genotypes (B73, NIL 1, and NIL 2) at V8 and VT growth stages (b), irrespective of management.

# Make your phyloseq object
mydata <- import_biom(BIOMfilename = "taxonomy.feature-table.biom", treefilename = "tree.nwk", parseFunction = parse_taxonomy_default)
mapfile <- import_qiime_sample_data("metadata.txt")
physeqITS <- merge_phyloseq(mydata, mapfile)

# Subset the levels of the factors
VT_N0 <- subset_samples(physeqITS, Growth.stage == "VT" | N.dosages == "N:0" | Inoculant == "none")

## a)
mm_lefse <- run_lefse(
  VT_N0,
  wilcoxon_cutoff = 0.01,
  group = "Genotype",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4
)
plot_cladogram(mm_lefse, color = MG)

## b)

library(ggplot2)

# Subset taxa based on specified Phyla
top <- subset_taxa(physeqITS, Phylum %in% c("Ascomycota", "Basidiomycota", "Mortierellomycota"))
top <- transform_sample_counts(top, function(OTU) OTU / sum(OTU))

# Agglomerate taxa at the Genus level
ph.glom <- tax_glom(top, taxrank = "Genus")
# Prune taxa with no abundance
ph.glom_1 <- prune_taxa(taxa_sums(ph.glom) > 0, ph.glom)
ph.glom_C <- tax_glom(ph.glom_1, taxrank = "Genus")
Genus10 <- names(sort(taxa_sums(ph.glom_C), decreasing = TRUE)[1:20])
top <- prune_taxa(Genus10, ph.glom_C)

# Rename columns in the tax_table
colnames(tax_table(top)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Plotting
o <- plot_bar(top, x = "Genotype", fill = "Phylum") +
  scale_fill_manual(values = MG) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        strip.text.x = element_text(family = "Arial", size = 12, color = "black", face = "bold"),
        strip.text.y = element_text(family = "Arial", size = 12, color = "black", face = "bold"),
        strip.background = element_rect(colour = "black", fill = "gray98")) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "", y = "Relative abundance (%)") +
  facet_grid(cols = vars(paste(Inoculant, N.dosages, sep = "_")), 
             rows = vars(Phylum), 
             scales = "free_y", switch = 'y') +
  ggtitle("Dominant Genera Across Genotypes") +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"))

o


