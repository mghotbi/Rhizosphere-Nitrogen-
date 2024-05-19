#a
# Re-build the phyloseq object with the sum of the factors (sum up replications)
#introduce the Rep as a factor
Rep<-as.factor(physeq_rel22)$Rep
#summerize metadata based on replications
sum.net2=physeq_rel22 %>% 
  group_by(Rep) %>% 
  select_if(is.numeric) %>%
  summarise_all(list(mean = ~mean(.)))


#alternatively create a new column, remove duplicates and merge it with your phyloseq obj
vec <- interaction(Genotype, Inoculant, N.dosages, Gorwth.stages)
new_df <- meta22 %>%
  mutate(newcol = vec) %>%
  distinct(newcol, .keep_all = TRUE) %>%
  na.omit()


#now re-merge new col of metadata with biome file

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

### b
##individual fig for each genotype

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(gridExtra)
nit.glom@sam_data
nit_melt<-psmelt(nit.glom)
write.csv(nit_melt,"nit_melt.csv")
PotentialN22<-read.csv("PotentialN22_N67_4.csv")
nit_melt<-read.csv("nit_melt.csv")

Final3 <-  dplyr::left_join(PotentialN22,nit_melt, by="sample.id")%>%
  group_by(sample.id) %>%
  arrange(desc(Potential_N_micro)) 
write.csv(Final3,"Final3.csv")
                                        
#created using left_join by sample.id/common in both files
dt <- read.csv("Final3.csv")

# Ensure necessary columns are factors
dt$Genotype <- as.factor(dt$Genotype)
dt$N.dosage <- as.factor(dt$N.dosage)
dt$Inoculant <- as.factor(dt$Inoculant)
dt$management <- as.factor(dt$management)
dt$Growth.stage <- as.factor(dt$Growth.stage)
print(levels(dt$management))

# Define the custom theme
mytheme <- theme_bw() + 
  theme( panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = 'white', color = "#e1deda"),
    panel.border = element_rect(linetype = "solid", fill = NA),
    axis.line.x = element_line(colour = 'black', size = 0.6),
    axis.line.y = element_line(colour = 'black', size = 0.6),
    axis.ticks = element_line(colour = 'black', size = 0.35),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12), 
    legend.key.size = unit(0.8, 'cm'),
    axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(family = "Times New Roman", size = 11, angle = 0, color = "black", face = "bold"),
    axis.text.y = element_text(family = "Times New Roman", size = 11, color = "black", face = "bold"),
    plot.title = element_text(color = "black", size = 12, face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
    strip.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"),
    legend.position = "top",
    legend.key = element_blank()
  )

create_plot <- function(data, genotype, management_status) {
  # Filter data based on genotype and management status
  filtered_data <- data %>%
    filter(Genotype == genotype & management == management_status)
  
  if (nrow(filtered_data) == 0) {
    return(NULL)
  }
  
  # maximum y-value for combination of Growth.stage and management
  max_y_values <- filtered_data %>%
    group_by(Growth.stage) %>%
    summarise(max_y = max(log1p(Abundance), na.rm = TRUE)) %>%
    ungroup()
  
  # Create df_label DataFrame and join with max_y_values to get the max_y for each facet
  df_label <- filtered_data %>%
    group_by(Growth.stage) %>%
    summarise(
      Inter = if(n_distinct(Potential_N_micro) > 1) lm(log1p(Abundance) ~ log1p(Potential_N_micro))$coefficients[1] else NA,
      Coeff = if(n_distinct(Potential_N_micro) > 1) lm(log1p(Abundance) ~ log1p(Potential_N_micro))$coefficients[2] else NA,
      pval = if(n_distinct(Potential_N_micro) > 1) summary(lm(log1p(Abundance) ~ log1p(Potential_N_micro)))$coefficients[2, 4] else NA,
      r2 = if(n_distinct(Potential_N_micro) > 1) summary(lm(log1p(Abundance) ~ log1p(Potential_N_micro)))$r.squared else NA
    ) %>%
    ungroup() %>%
    mutate(
      Label = paste(
        "bold(italic(y))==bold(", round(Inter, 3), ")", 
        ifelse(Coeff < 0, " - ", " + "), 
        "bold(", round(abs(Coeff), 3), ") * bold(italic(x))",
        "~~~~bold(italic(R^2))==bold(", round(r2, 3), ")",
        "~~bold(italic(p))==bold(", round(pval, 3), ")", sep = ""
      )
    ) %>%
    left_join(max_y_values, by = "Growth.stage")
  
  df_label <- df_label %>%
    group_by(Growth.stage) %>%
    mutate(
      y_pos = max_y - seq(0, length.out = n()) * (max_y / (n() + 1))  # Adjust spacing by adding +1 to n()
    ) %>%
    ungroup()
  
  # if you need to auto-adjust the x position for labels/*** edit it as you need
  x_pos <- max(log1p(filtered_data$Potential_N_micro), na.rm = TRUE) * 0.65  # Move left by setting factor < 0.75
  
  # Plotting
  p <- ggplot(filtered_data, aes(x = log1p(Potential_N_micro), y = log1p(Abundance), color = Growth.stage, group = Growth.stage)) +
    geom_smooth(aes(color = Growth.stage, fill = Growth.stage), method = "lm", se = TRUE) +  # Change method to lm for linear model
    geom_point(size = 2) +  # Reduce the size of points
    geom_text_repel(data = df_label,
                    aes(x = x_pos, y = y_pos,  # Adjust x and y positions
                        label = Label, color = Growth.stage), 
                    show.legend = FALSE, parse = TRUE, size = 5,  # Increase the size of the equations
                    box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, 
                    nudge_y = -0.15) +  # Nudge labels down slightly
    scale_color_manual(values = MG) +
    scale_fill_manual(values = MG) +
    ggtitle(paste("Genotype:", genotype, "Management:", management_status)) +
    labs(x = "ln(Potential Nitrification)", y = "ln(Nitrifier Abundance)") +  # Add "ln" to axis labels
    mytheme +  # Apply the custom theme
    guides(color = guide_legend(override.aes = list(size = 2)))  # Adjust the size of points in the legend
  
  return(p)}
#please recall I have already introduced mytheme and MG (color pallet)
# Genotypes 
genotypes <- c("NIL 1", "NIL 2", "B73")

# save plots 
for (genotype in genotypes) {
  dt_genotype <- filter(dt, Genotype == genotype)
  management_levels <- unique(dt_genotype$management)
  
  plots <- list()
  for (management_status in management_levels) {
    plot_name <- paste(genotype, management_status, sep = "_")
    print(paste("Creating plot for:", plot_name))  # Debug print statement
    plot <- create_plot(dt_genotype, genotype, management_status)
    if (!is.null(plot)) {
      plots[[plot_name]] <- plot
    }
  }
   if (length(plots) > 0) {
    # Arrange plots in a grid
    grid_plot <- marrangeGrob(plots, nrow = 2, ncol = 2, top = paste("Genotype:", genotype))
    
    # Save the grid plot to a file
    print(paste("Saving plots to file for:", genotype))  # Debug print statement
    ggsave(filename = paste0(genotype, "_management_plots.tiff"), plot = grid_plot, width = 16, height = 12, units = "in", dpi = 600)
  } else {
    print(paste("No plots were created for:", genotype))
  }}


  ## all together 
  #############The best one :)
nit.glom@sam_data
nit_melt<-psmelt(nit.glom)
write.csv(nit_melt,"nit_mel1t.csv")
PotentialN22<-read.csv("PotentialN22_N67_4.csv")
nit_melt<-read.csv("nit_mel1t.csv")
dt<-read.csv("dt.csv")

library(tidyr)
dto <- dt %>%
  group_by(sample.id) %>%
  arrange(desc(Potential_N_micro)) 


#############
MG <- c("#66a182","#2e4057","#8d96a0","#0e669b","#00798c","dodgerblue4", "steelblue2","lightskyblue4","#82cfd0","#b2e0e4","honeydew3","mintcream","#8d96a3","lavender","#CC6686","lavenderblush2","mistyrose3","#e1deda","darkgoldenrod","burlywood","papayawhip","wheat4","cornsilk3","khaki2","beige","gray60","gray80","gray96")

# Define your custom theme
mytheme <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(), # gets rid of grey and lines in the middle
    panel.grid.major = element_blank(), # gets rid of grey and lines in the middle
    panel.background = element_rect(fill = "white"), # makes entire background white
    plot.background = element_rect(fill = 'white', color = "#e1deda"),
    panel.border = element_rect(linetype = "solid", fill = NA), # gets rid of square going around the entire graph
    axis.line.x = element_line(colour = 'black', size = 0.6), # sets the axis line size
    axis.line.y = element_line(colour = 'black', size = 0.6), # sets the axis line size
    axis.ticks = element_line(colour = 'black', size = 0.35),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12), 
    legend.key.size = unit(0.8, 'cm'), # sets the tick lines
    axis.title.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), # size of x-axis title
    axis.title.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), # size of y-axis title
    axis.text.x = element_text(family = "Times New Roman", size = 11, angle = 0, color = "black", face = "bold"), # size of x-axis text
    axis.text.y = element_text(family = "Times New Roman", size = 11, color = "black", face = "bold"), # size of y-axis text
    plot.title = element_text(color = "black", size = 12, face = "bold"),
    plot.subtitle = element_text(size = 11),
    strip.background = element_rect(colour = "white", fill = "white"), # set strip background color to white
    strip.text.x = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), # set strip text size to 12 for x
    strip.text.y = element_text(family = "Times New Roman", size = 12, color = "black", face = "bold"), # set strip text size to 12 for y
    legend.position = "top", # set legend position to top
    legend.key = element_blank() # remove legend key
  )

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# Read data
dt <- read.csv("Final.csv")

# Ensure necessary columns are factors
dt$Genotype <- as.factor(dt$Genotype)
dt$N.dosage <- as.factor(dt$N.dosage)
dt$Inoculant <- as.factor(dt$Inoculant)
dt$management <- as.factor(dt$management)
dt$Growth.stage <- as.factor(dt$Growth.stage)

max_y_values <- dt %>%
  group_by(Growth.stage, management) %>%
  summarise(max_y = max(sqrt(Abundance), na.rm = TRUE)) %>%
  ungroup()

df_label <- dt %>%
  group_by(Genotype, management, Growth.stage) %>%
  summarise(
    Inter = if(n_distinct(Potential_N_micro) > 1) lm(Abundance ~ Potential_N_micro)$coefficients[1] else NA,
    Coeff = if(n_distinct(Potential_N_micro) > 1) lm(Abundance ~ Potential_N_micro)$coefficients[2] else NA,
    pval = if(n_distinct(Potential_N_micro) > 1) summary(lm(Abundance ~ Potential_N_micro))$coefficients[2, 4] else NA,
    r2 = if(n_distinct(Potential_N_micro) > 1) summary(lm(Abundance ~ Potential_N_micro))$r.squared else NA
  ) %>%
  ungroup() %>%
  mutate(
    Label = paste(
      "bold(italic(y))==bold(", round(Inter, 3), ")", 
      ifelse(Coeff < 0, " - ", " + "), 
      "bold(", round(abs(Coeff), 3), ") * bold(italic(x))",
      "~~~~bold(italic(R^2))==bold(", round(r2, 3), ")",
      "~~bold(italic(p))==bold(", round(pval, 3), ")", sep = ""
    )  ) %>%
  left_join(max_y_values, by = c("Growth.stage", "management"))

df_label <- df_label %>%
  group_by(Growth.stage, management) %>%
  mutate(
    y_pos = max_y - seq(0, length.out = n()) * (max_y / (n() + 1))  # Adjust spacing by adding +1 to n()
  ) %>%
  ungroup()

x_pos <- max(log1p(dt$Potential_N_micro), na.rm = TRUE) * 0.65  # Move left by setting factor < 0.75

# Plotting
NILs_B73 <- ggplot(dt, aes(x = log1p(Potential_N_micro), y = log1p(Abundance), color = Genotype, group = Genotype)) +
  geom_smooth(aes(color = Genotype, fill = Genotype), method = "lm", se = TRUE) +  # Change method to lm for linear model
  geom_point(size = 2) +  # Reduce the size of points
  geom_text_repel(data = df_label,
                  aes(x = x_pos, y = y_pos,  # Adjust x and y positions
                      label = Label, color = Genotype), 
                  show.legend = FALSE, parse = TRUE, size = 3,  # Make equations smaller
                  box.padding = 0.3, point.padding = 0.3, max.overlaps = 10, 
                  nudge_y = -0.15) +  # Nudge labels down slightly
  scale_color_manual(values = MG) +
  scale_fill_manual(values = MG) +
  ggtitle("") +
  facet_grid(rows = vars(Growth.stage), cols = vars(management), scales = "free_y", switch = 'x') +
  labs(x = "Potential Nitrification", y = "Nitrifier Abundance") +
  mytheme +  # Apply the custom theme
  guides(color = guide_legend(override.aes = list(size = 2)))  # Adjust the size of points in the legend
print(NILs_B73)
