# Re-build the phyloseq object with the sum of the factors (sum up replications)
#introduce the Rep as a factor
Rep<-as.factor(physeq_rel22)$Rep
#summerize metadata based on replications
sum.net2=physeq_rel22 %>% 
  group_by(Rep) %>% 
  select_if(is.numeric) %>%
  summarise_all(list(mean = ~mean(.)))


#alternatively create a new column and merge it with your phyloseq obj
vec <- interaction(Genotype, Inoculant, N.dosages, Gorwth.stages)
Meta22 %>%
  mutate(newcol = vec)



#now re-merge metadata with biome file

# Compositional data
physeq_rel22 <- transform_sample_counts(physeq22, function(x) 100 * x / sum(x))
saveRDS(physeq_rel22, file = "physeq_rel22.rds")
physeq_rel22 <- readRDS(file = "physeq_rel22.rds")

# Extract Nitrifiers
Nitrifiers.tax22 <- subset_taxa(physeq_rel22, Genus %in% c("Nitrosospira", "Candidatus_Nitrososphaera", "Nitrososphaeraceae", "Nitrospira_japonica", "Nitrobacter", "Nitrosomonas", "Nitrosospira", "Nitrospira"))

# Subset V8 and VT
nit_V8 <- subset_samples(Nitrifiers.tax22, Growth.stages == "V8")
nit_VT <- subset_samples(Nitrifiers.tax22, Growth.stages == "VT")

# Using glom function for clean-up
nit.glom22_V8 <- tax_glom(nit_V8, taxrank = "Genus")
nit.glom22_V8_1 <- prune_taxa(taxa_sums(nit.glom22_V8) > 0, nit.glom22_V8)
# Get top 20 Genus
Genus10 <- names(sort(taxa_sums(nit.glom22_V8_1), TRUE)[1:20])
top_v5 <- prune_taxa(Genus10, nit.glom22_V8_C)

# reintroduce factros
sample_data(top_v5)$Growth.stage <- factor(sample_data(physeq_rel22)$Growth.stage, levels = c("V8", "VT"))
sample_data(top_v5)$N.dosage <- factor(sample_data(physeq_rel22)$N.dosage, levels = c("N:0", "N:67"))
sample_data(top_v5)$Genotype <- factor(sample_data(physeq_rel22)$Genotype, levels = c("B73", "NIL 1", "NIL 2"))
sample_data(top_v5)$Inoculant <- factor(sample_data(physeq_rel22)$Inoculant, levels = c("None", "Proven"))

                                        
# Function to plot bar charts
plot_bar_chart <- function(data, stage) {
  p <- plot_bar(data, x = "Genotype", fill = "Genus") +
    scale_fill_manual(values = MG) +
    theme(legend.position = "bottom",
          strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold")) +
    labs(x = "", y = "Relative abundance") +
    facet_grid(cols = vars(paste(Inoculant, N.dosages, sep = "_")), scales = "free_y", switch = 'x') +
    ggtitle(paste("Fungi abundance at", stage, "growth stage")) +
    mytheme
  return(p)
}


# Plotting
plot_v8 <- plot_bar_chart(top_v5_v8, stage = "V8")
plot_vt <- plot_bar_chart(top_v5_vt, stage = "VT")

# Arrange plots
ggarrange(plot_v8, plot_vt, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")

