
library(dplyr)
library(purrr)
library(tidyverse)

topcons2 <- read_tsv("TOPCONS2_length_2.tsv")

# Función para calcular índice de Jaccard por proteína
calc_jaccard <- function(annotated, predicted) {
  total_len_annotated <- sum(annotated$End - annotated$Start + 1)
  total_len_predicted <- sum(predicted$End - predicted$Start + 1)
  
  # Caso 2: no hay anotaciones, pero sí predicciones → score = 0
  if (nrow(annotated) == 0) return(0)
  
  # Caso 3: no hay predicciones, pero sí anotaciones → score = 0
  if (nrow(predicted) == 0) return(0)
  
  # Caso 4: hay ambos → calcular intersección
  total_intersection <- 0
  for (i in seq_len(nrow(annotated))) {
    seg1 <- annotated[i, ]
    overlaps <- predicted %>%
      filter(Start <= seg1$End & End >= seg1$Start)
    
    for (j in seq_len(nrow(overlaps))) {
      seg2 <- overlaps[j, ]
      overlap_len <- min(seg1$End, seg2$End) - max(seg1$Start, seg2$Start) + 1
      if (overlap_len > 0) {
        total_intersection <- total_intersection + overlap_len
      }
    }
  }
  
  total_union <- total_len_annotated + total_len_predicted - total_intersection
  return(total_intersection / total_union)
}

protein_ids <- union(unique(Uniprot_domains$ID), unique(topcons2$ID))
protein_ids_1 <- union(unique(Uniprot_domains$ID), unique(df_tm$ID))

# Calcular Jaccard por proteína
jaccard_results <- map_dfr(protein_ids, function(pid) {
  ann <- Uniprot_domains %>% filter(ID == pid)
  pred <- topcons2 %>% filter(ID == pid)
  score <- calc_jaccard(ann, pred)
  tibble(ID = pid, Jaccard = score)
})


jaccard_results_1 <- map_dfr(protein_ids_1, function(pid) {
  ann <- Uniprot_domains %>% filter(ID == pid)
  pred <- df_tm %>% filter(ID == pid)
  score <- calc_jaccard(ann, pred)
  tibble(ID = pid, Jaccard = score)
})

jaccard_results$Software <- "TOPCONS2"
jaccard_results_1$Software <- "DeepTMHMM"


library(data.table)

# Convertir a data.table
dt1 <- as.data.table(jaccard_results)
dt2 <- as.data.table(jaccard_results_1)

# Suponiendo que ambas tablas tienen una columna "ID"
# Hacer merge por ID, conservando todos los IDs de ambos
merged <- merge(dt1, dt2, by = "ID", all = TRUE, suffixes = c("_res1", "_res2"))

merged$Software_res2 <- "DeepTMHMM"

# Reemplazar NA por 1 en los Jaccard faltantes
merged[is.na(Jaccard_res1), Jaccard_res1 := 1]
merged[is.na(Jaccard_res2), Jaccard_res2 := 1]


merged[, diferencia := Jaccard_res1 - Jaccard_res2]


diferencias_2.1 <- merged$Jaccard_res2 - merged$Jaccard_res1 
# Prueba de normalidad
shapiro.test(diferencias_2.1)

qqnorm(diferencias_2.1)
qqline(diferencias_2.1, col = "red")


##TEST DE WILCOXON
# Filtrar por software
deep_jaccard_2 <- merged %>% filter(Software_res2 == "DeepTMHMM")
topcons_jaccard_2 <- merged %>% filter(Software_res1 == "TOPCONS2")

# Asegurar orden
deep_jaccard_2 <- deep_jaccard_2[order(deep_jaccard_2$ID), ]
topcons_jaccard_2 <- topcons_jaccard_2[order(topcons_jaccard_2$ID), ]

# Comparación pareada
wilcox.test(deep_jaccard_2$Jaccard_res2, topcons_jaccard_2$Jaccard_res1, paired = TRUE)


###

score_final_b <- data.frame(ID = merged$ID, Score = merged$Jaccard_res1, 
                          Software = merged$Software_res1)

score_final_b$Method <- "Jaccard"
score_final_c <- data.frame(ID = merged$ID, Score = merged$Jaccard_res2, 
                            Software = merged$Software_res2)

score_final_c$Method <- "Jaccard"

score_test_2 <- rbind(score_final_b, score_final_c) #usar este para gráfico





# Diferencias entre DeepTMHMM y TOPCONS2
diferencias_2 <- jaccard_results$Jaccard - jaccard_results_1$Jaccard

# Prueba de normalidad
shapiro.test(diferencias_2)

qqnorm(diferencias)
qqline(diferencias, col = "red")


##TEST DE WILCOXON
# Filtrar por software
deep_jaccard <- jaccard_results_1 %>% filter(Software == "DeepTMHMM")
topcons_jaccard <- jaccard_results %>% filter(Software == "TOPCONS2")

# Asegurar orden
deep_jaccard <- deep_jaccard[order(deep_jaccard$ID), ]
topcons_jaccard <- topcons_jaccard[order(topcons_jaccard$ID), ]

# Comparación pareada
wilcox.test(deep_jaccard$Jaccard, topcons_jaccard$Jaccard, paired = TRUE)





