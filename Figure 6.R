# Tidy phyloseq object# you can find my tidy_phyloseq object in the_little_things branch
glom22_vt_net <- tidy_phyloseq(glom22_vt_net)

# Subset taxa for network
NETWORKglom22_vt_net <- subset_taxa(glom22_vt_net, Species %in% c("uncultured_Gemmatimonas", "Conexibacter_sp.", "Vitiosangium_sp.", "Trichoderma_viride", "Mortierella_elongata", "uncultured_crenarchaeote", "SBR1031"))

# Save tidy phyloseq object
saveRDS(NETWORKglom22_vt_net, "NETWORKglom22_vt_net.rds")

# Define species colors
species_colors <- c(
  "uncultured_Gemmatimonas" = "#66a182",
  "Conexibacter_sp." = "#2e4057",
  "Vitiosangium_sp." = "#8d96a0",
  "Trichoderma_viride" = "olivedrab3",
  "Mortierella_elongata" = "grey",
  "uncultured_crenarchaeote" = "#0e669b",
  "uncultured_Gemmatimonas" = "dodgerblue")

# Subset samples for high N level
NETWORKglom22_vt_Nhigh <- subset_samples(NETWORKglom22_vt_net, N.dosages == "N:67" & Inoculant == "Proven")

# Perform taxonomic aggregation
ph.glom21_V5 <- tax_glom(NETWORKglom22_vt_Nhigh, taxrank = "Species")
ph.glom21_V5_1 <- prune_taxa(taxa_sums(ph.glom21_V5) > 0, ph.glom21_V5)
ph.glom21_V5_C <- tax_glom(ph.glom21_V5_1, taxrank = "Genus")
Genus10 <- names(sort(taxa_sums(ph.glom21_V5_C), decreasing = TRUE)[1:10])
gen10 <- prune_taxa(Genus10, ph.glom21_V5_C)
top_v5 <- prune_taxa(taxa_sums(gen10) > 0, gen10)
colnames(tax_table(top_v5)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Factorize variables
sample_data(top_v5)$Growth.stage <- factor(sample_data(physeq_relsum22)$Growth.stage, levels = c("V8", "VT"))
sample_data(top_v5)$N.dosage <- factor(sample_data(physeq_relsum22)$N.dosage, levels = c("N:0", "N:67"))
sample_data(top_v5)$Genotype <- factor(sample_data(physeq_relsum22)$Genotype, levels = c("B73", "NIL 1", "NIL 2"))
sample_data(top_v5)$Inoculant <- factor(sample_data(physeq_relsum22)$Inoculant, levels = c("None", "Proven"))

# Plot bar chart
o <- plot_bar(top_v5, x = "Genotype", fill = "Phylum") +
  scale_fill_manual(values = c(
    "Gemmatimonadota" = "#66a182",
    "Actinobacteriota" = "#2e4057",
    "Myxococcota" = "#8d96a0",
    "Ascomycota" = "olivedrab3",
    "Mortierellomycota" = "grey",
    "Crenarchaeota" = "#0e669b",
    "Gemmatimonadota" = "dodgerblue"  )) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2)) +
  labs(x = "", y = "Relative abundance %") +
  facet_grid(cols = vars(paste(Inoculant, N.dosages, sep = "_")), scales = "free_y", switch = 'x') +
  mytheme +  theme(strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"))

# Additional theme adjustments
ooo <- o +   theme(strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold")) +
  theme(strip.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"))

v <- ooo +  theme(strip.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"))

k1 <- v + theme(strip.background = element_rect(colour = "black", fill = "gray98")) +
  theme(strip.background = element_rect(colour = "white", fill = "gray97"),
        strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
        legend.position = "right",
        legend.key = element_blank())

k1

# Additional processing for another set of data

NETWORKglom22_vt_Non <- subset_samples(NETWORKglom22_vt_net, N.dosages == "N:0" & Inoculant == "None")

ph.glom21_V5 <- tax_glom(NETWORKglom22_vt_Non, taxrank = "Genus")
ph.glom21_V5_1 <- prune_taxa(taxa_sums(ph.glom21_V5) > 0, ph.glom21_V5)
ph.glom21_V5_C <- tax_glom(ph.glom21_V5_1, taxrank = "Genus")
Genus10 <- names(sort(taxa_sums(ph.glom21_V5_C), decreasing = TRUE)[1:10])
gen10 <- prune_taxa(Genus10, ph.glom21_V5_C)
top_v5 <- prune_taxa(taxa_sums(gen10) > 0, gen10)

colnames(tax_table(top_v5)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

sample_data(top_v5)$Growth.stage <- factor(sample_data(physeq_relsum22)$Growth.stage, levels = c("V8", "VT"))
sample_data(top_v5)$N.dosage <- factor(sample_data(physeq_relsum22)$N.dosage, levels = c("N:0", "N:67"))
sample_data(top_v5)$Genotype <- factor(sample_data(physeq_relsum22)$Genotype, levels = c("B73", "NIL 1", "NIL 2"))
sample_data(top_v5)$Inoculant <- factor(sample_data(physeq_relsum22)$Inoculant, levels = c("None", "Proven"))

top_v5 <- subset_taxa(top_v5, Genus != "uncultured")

o <- plot_bar(top_v5, x = "Genotype", fill = "Genus") +
  scale_fill_manual(values = c(
    "Gemmatimonas" = "#66a182",
    "Conexibacter" = "#2e4057",
    "Archangium" = "#8d96a0",
    "Trichoderma" = "olivedrab3",
    "Mortierella" = "tomato3",
    "Nitrososphaeraceae" = "#0e669b",
    "SBR1031" = "#F77F00"  )) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3)) +
  labs(x = "", y = "Relative abundance %") +
  facet_grid(cols = vars(paste(Inoculant, N.dosages, sep = "_")), scales = "free_y", switch = 'x') +
  mytheme +  theme(strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"))
ooo <- o +   theme(strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"))
v <- ooo +   theme(strip.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"))
k2 <- v +   theme(strip.background = element_rect(colour = "black", fill = "gray98")) +
  theme(strip.background = element_rect(colour = "white", fill = "gray97"),
        strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
        legend.position = "right",
        legend.key = element_blank())

k2
