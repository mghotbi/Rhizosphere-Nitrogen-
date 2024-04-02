#Figure 9
cors <- function(df) {
  M <- Hmisc::rcorr(as.matrix(df)) 
  Mdf <- map(M, ~data.frame(.x)) return(Mdf) }
formatted_cors <- function(df){
  cors(df) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, "measure2")) %>% 
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA)) }


#  libraries
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
library(RColorBrewer)

# code
mat <- as.matrix(cor(meta, ASV))
wrote
pheatmap::pheatmap(mat, 
                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
)
