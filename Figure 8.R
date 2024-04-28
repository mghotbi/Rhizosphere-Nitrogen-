#calculate R2 and p
lm_fit <- lm(log(Grain.N.concentration) ~ log(VT.Plant.N.Content), data = Plantdata22NIL1)
lm_fit$model
# Extract R-squared and p-value
r_squared <- summary(lm_fit)$r.squared
p_value <- summary(lm_fit)$coefficients[8]
# Create the R2 AND p-value label by paste
p_label <- paste("P value =", formatC(p_value, digits = 2))
R_label <- paste("R² =", formatC(r_squared, digits = 2))


# Create ggplots
NIL2 <- ggplot(Plantdata22NIL2, aes(x = log(VT.Plant.N.Content), y = log(Grain.N.concentration), color = management)) +
  geom_point(size = 2) +  
  geom_smooth(method = lm, se = TRUE, fill = "gray93") + 
  scale_size_manual(values = c(1.5, 1.5, 1.5)) +
  scale_shape_manual(values = c(3, 16, 17)) +
  scale_color_manual(values = MG) +
  labs(y = NULL, x = NULL) + 
  theme(axis.text.y = element_blank(), legend.position = "bottom") +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "top",
               formula = y ~ x, parse = TRUE) +
  ggtitle("NIL 2") + mytheme

# Create ggplot for B73
B73 <- ggplot(Plantdata22B73, aes(x = log(VT.Plant.N.Content), y = log(Grain.N.concentration), color = management)) +
  geom_point(size = 2) +  
  geom_smooth(method = lm, se = TRUE, fill = "gray93") + 
  scale_size_manual(values = c(1.5, 1.5, 1.5)) +
  scale_shape_manual(values = c(19, 19, 19)) +
  scale_color_manual(values = MG) +
  labs(y = NULL, x = NULL) + 
  theme(axis.text.y = element_blank(), legend.position = "bottom") +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "top",
               formula = y ~ x, parse = TRUE) +
  ggtitle("B73") + mytheme

# Create ggplot for NIL1
NIL1 <- ggplot(Plantdata22NIL1, aes(x = log(VT.Plant.N.Content), y = log(Grain.N.concentration), color = management)) +
  geom_point(size = 2) +  
  geom_smooth(method = glm, se = TRUE, fill = "gray93") + 
  scale_size_manual(values = c(1.5, 1.5, 1.5)) +
  scale_shape_manual(values = c(19, 19, 19)) +
  scale_color_manual(values = MG) +
  labs(y = NULL, x = NULL) + 
  theme(axis.text.y = element_blank(), legend.position = "bottom") +
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.x.npc = "right", label.y.npc = "top",
               formula = y ~ x, parse = TRUE) +
  ggtitle("NIL 1") + mytheme

# Arrange ggplots
figure <- ggarrange(B73, NIL1, NIL2, ncol = 3, nrow = 1, common.legend = TRUE, align = "hv")
