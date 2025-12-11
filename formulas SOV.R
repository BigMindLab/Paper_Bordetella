

setwd("C:/Users/jalme/OneDrive/Escritorio/Immunoinformatics")
install.packages("ggsignif")
library(ggsignif)
library(dplyr)
library(purrr)

calc_sov99_zemla <- function(annotated, predicted) {
  N_total <- 0
  sov_sum_99 <- 0
  
  if (nrow(annotated) == 0) return(NA)
  
  for (i in seq_len(nrow(annotated))) {
    seg1 <- annotated[i, ]
    
    # Buscar segmentos predichos que se superponen
    overlaps <- predicted %>% filter(Start <= seg1$End & End >= seg1$Start)
    
    len1 <- seg1$End - seg1$Start + 1
    N_total <- N_total + len1
    
    if (nrow(overlaps) == 0) next  # No hay superposición
    
    for (j in seq_len(nrow(overlaps))) {
      seg2 <- overlaps[j, ]
      
      # calcular intersección
      minov <- min(seg1$End, seg2$End) - max(seg1$Start, seg2$Start) + 1
      if (minov <= 0) next
      
      maxov <- max(seg1$End, seg2$End) - min(seg1$Start, seg2$Start) + 1
      len2 <- seg2$End - seg2$Start + 1
      
      # calcular delta según Zemla (1999)
      delta <- min(
        c(
          minov,
          maxov - minov,
          floor(len1 / 2),
          floor(len2 / 2)
        )
      )
      
      # SOV'99
      sov_sum_99 <- sov_sum_99 + ((minov + delta) * len1) / maxov
    }
  }
  
  if (N_total == 0) return(NA)
  
  return(sov_sum_99 / N_total)
}

# Asegurarse que IDs coincidan
df_tm$ID <- gsub("^#\\s*sp\\_|\\_.*$", "", df_tm$ID)
protein_ids <- union(unique(Uniprot_domains$ID), unique(topcons2$ID))

# Aplicar por proteína #TOPCONS2
sov99_results <- map_dfr(protein_ids, function(pid) {
  ann <- Uniprot_domains %>% filter(ID == pid)
  pred <- topcons2 %>% filter(ID == pid)
  score <- calc_sov99_zemla(ann, pred)
  tibble(ID = pid, SOV99 = score)
})



# protein_ids <- intersect(unique(Uniprot_domains$ID), unique(df_tm$ID))


print(protein_ids)


# Aplicar por proteína #DeepTMHMM
protein_ids_1 <- union(unique(Uniprot_domains$ID), unique(df_tm$ID))
sov99_results_1 <- map_dfr(protein_ids_1, function(pid) {
  ann <- Uniprot_domains %>% filter(ID == pid)
  pred <- df_tm %>% filter(ID == pid)
  score <- calc_sov99_zemla(ann, pred)
  tibble(ID = pid, SOV99 = score)
})


sov99_results$Software <- "TOPCONS2"
sov99_results_1$Software <- "DeepTMHMM"

sov99_merge <- rbind(sov99_results, sov99_results_1)

sov99_merge$SOV99 <- ifelse(sov99_merge$SOV99 > 1, 1, sov99_merge$SOV99)



#PLOT

plot <- ggplot(data = sov99_merge, aes(x = Software, y = SOV99)) +
  theme_bw() +
  geom_violin(scale = "width", alpha = 0.15, aes(fill = Software),
              show.legend = FALSE, trim = FALSE, linewidth = 0.4) +
  geom_point(aes(fill = Software), shape = 21, size = 2, alpha = 0.8,
             show.legend = FALSE, position = position_jitter(width = 0.05)) +
  geom_boxplot(width = 0.05, fill = "grey90", outlier.shape = NA) +
  geom_errorbar(stat = "boxplot", width = 0.2, linewidth = 1, color = "red",
                aes(ymin = after_stat(middle), ymax = after_stat(middle))) +
  scale_fill_manual(values = c("DeepTMHMM" = "#b9ffffff", "TOPCONS2" = "#fab8c2ff")) +
  scale_y_continuous(breaks = seq(0, 1.13, by = 0.1), limits = c(0.4, 1.13)) +
  labs(y = expression("Score")) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, vjust = 2)) +
  geom_signif(comparisons = list(c("DeepTMHMM", "TOPCONS2")),
              map_signif_level = TRUE,
              textsize = 5,
              y_position = 1.05)  # Ajusta la altura de la barra
plot       


score_final <- data.frame(ID = merged_SOV99$ID, Score = merged_SOV99$SOV99_res1, 
                          Software = merged_SOV99$Software_res1)

score_final$Method <- "SOV99"
score_final_a <- data.frame(ID = merged_SOV99$ID, Score = merged_SOV99$SOV99_res2, 
                            Software = merged_SOV99$Software_res2)

score_final_a$Method <- "SOV99"

score_test_1 <- rbind(score_final, score_final_a) #usar este para gráfico




# Diferencias entre DeepTMHMM y TOPCONS2
diferencias <- sov99_results_1$SOV99 - sov99_results$SOV99

# Prueba de normalidad
shapiro.test(diferencias)

qqnorm(diferencias)
qqline(diferencias, col = "red")


##TEST DE WILCOXON
# Filtrar por software
deep <- sov99_results_1 %>% filter(Software == "DeepTMHMM")
topcons <- sov99_results %>% filter(Software == "TOPCONS2")

# Asegurar orden
deep <- deep[order(deep$ID), ]
topcons <- topcons[order(topcons$ID), ]

# Comparación pareada
wilcox.test(deep$SOV99, topcons$SOV99, paired = TRUE)


# Convertir a data.table
dt1_SOV99 <- as.data.table(sov99_results)
dt2_SOV99 <- as.data.table(sov99_results_1)

# Suponiendo que ambas tablas tienen una columna "ID"
# Hacer merge por ID, conservando todos los IDs de ambos
merged_SOV99 <- merge(dt1_SOV99, dt2_SOV99, by = "ID", all = TRUE,
                      suffixes = c("_res1", "_res2"))

merged_SOV99$SOV99_res1[is.na(merged_SOV99$SOV99_res1)] <- 0
merged_SOV99$SOV99_res2[is.na(merged_SOV99$SOV99_res2)] <- 0
merged_SOV99$Software_res2 <- "DeepTMHMM"
merged_SOV99$Software_res1 <- "TOPCONS2"


diferencias_SOV99 <- merged_SOV99$SOV99_res1 - merged_SOV99$SOV99_res2

# Prueba de normalidad
shapiro.test(diferencias_SOV99)

qqnorm(diferencias_SOV99)
qqline(diferencias_SOV99, col = "red")


##TEST DE WILCOXON
# Filtrar por software
deep_SOV99_2 <- merged_SOV99 %>% filter(Software_res2 == "DeepTMHMM")
topcons_SOV99_2 <- merged_SOV99 %>% filter(Software_res1 == "TOPCONS2")

# Asegurar orden
deep_SOV99_2 <- deep_SOV99_2[order(deep_SOV99_2$ID), ]
topcons_SOV99_2 <- topcons_SOV99_2[order(topcons_SOV99_2$ID), ]

# Comparación pareada
wilcox.test(deep_SOV99_2$SOV99_res2, topcons_SOV99_2$SOV99_res1, paired = TRUE)







