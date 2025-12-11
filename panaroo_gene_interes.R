
setwd("C:/Users/jalme/OneDrive/Escritorio/Immunoinformatics")
library(readr)
library(dplyr)

# Leer archivo de presencia/ausencia
pan <- read_csv("gene_presence_absence.csv")


# Contar cuántos genomas tienen cada gen
pan_categories <- pan %>%
  mutate(No_isolates = rowSums(!is.na(.[ , 4:ncol(.)])))
# Categorizar (suponiendo que ya sabes cuántos genomas hay)
n_genomes <- length(colnames(pan)) - 3  # Resta columnas no-genómicas

pan_categories_group <- pan_categories %>%
  mutate(Category = case_when(
    No_isolates >= 0.99 * n_genomes            ~ "Core",
    No_isolates >= 0.95 * n_genomes     ~ "Soft-core",
    No_isolates >= 0.15 * n_genomes     ~ "Shell",
    TRUE                                ~ "Cloud"
  ))

pan_categories_group %>% count(Category)

pan_categories_group$Category

genes_interes <- pan_categories_group %>%
  filter(Category %in% c("Core", "Soft-core"))

# Extraer solo genes Core y Soft-core gene_list
shell_genes <- pan_categories_group %>%
  filter(Category %in% c("Shell")) %>%
  pull(Gene)

# Guardarlos en un archivo
writeLines(shell_genes, "gene_shell_list.txt")




##ROARY
# Leer archivo de presencia/ausencia
pan_1 <- read_csv("gene_presence_absence_roary.csv")
n_genomes_1 <- length(grep("^GCF_|^GCA_", colnames(pan_1)))  # número total de genomas

pan_categories_1 <- pan_1 %>%
  mutate(Category = case_when(
    `No. isolates` == n_genomes_1            ~ "Core",
    `No. isolates` >= 0.95 * n_genomes_1     ~ "Soft-core",
    `No. isolates` >= 0.15 * n_genomes_1     ~ "Shell",
    TRUE                                   ~ "Cloud"
  ))

pan_categories_1 <- pan_1 %>%
  mutate(Category = case_when(
    `No. isolates` >= 0.99 * n_genomes_1            ~ "Core",
    `No. isolates` >= 0.95 * n_genomes_1     ~ "Soft-core",
    `No. isolates` >= 0.15 * n_genomes_1     ~ "Shell",
    TRUE                                   ~ "Cloud"
  ))

genes_interes_1 <- pan_categories_1 %>%
  filter(Category %in% c("Core", "Soft-core"))





# Resumen por tamaño del clúster (número de genomas donde está)
cluster_hist <- pan_categories_1 %>%
  count(`No. isolates`, Category)

ggplot(cluster_hist, aes(x = `No. isolates`, y = n, fill = Category)) +
  geom_col(width = 1, color = "black", linewidth = 0.1) +
  scale_fill_manual(values = c("Core" = "white", "Soft-core" = "orange", "Shell" = "gold", "Cloud" = "red")) +
  labs(x = "Number of genomes in cluster", y = "Number of gene clusters") +
  theme_bw()

# Crear los bins (puedes ajustar el tamaño con binwidth = 10, etc.)
cluster_hist_binned <- cluster_hist %>%
  mutate(bin = cut(`No. isolates`, breaks = seq(0, max(`No. isolates`)+51, by = 51), right = FALSE)) %>%
  group_by(bin, Category) %>%
  summarise(total_genes = sum(n), .groups = "drop")

bin_width <- 50 
cluster_hist_binned <- cluster_hist_binned %>%
  mutate(bin_mid = as.numeric(sub("\\[(\\d+),.*", "\\1", bin)) + bin_width / 2)
library(ggplot2)


# Graficar por bins
ggplot(cluster_hist_binned, aes(x = bin, y = total_genes, fill = Category)) +
  theme_bw() +
  geom_col(color = "black", linewidth = 0.1, position = "stack") +
  scale_fill_manual(values = c("Cloud" = "red", "Shell" = "orange", "Soft-core" = "yellow", "Core" = "white")) +
  scale_x_discrete(breaks = cluster_hist_binned$bin[seq(1, length(cluster_hist_binned$bin), by = 5)]) +
  labs(x = "Number of genomes per cluster (binned)", y = "Number of gene clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_blank())

ggplot(cluster_hist_binned, aes(x = bin, y = total_genes, fill = Category)) +
  theme_bw() +
  geom_col(color = "black", linewidth = 0.1, position = "stack") +
  scale_fill_manual(values = c("Cloud" = "red", 
                               "Shell" = "orange", 
                               "Soft-core" = "yellow", 
                               "Core" = "white")) +
  # Aquí cambiamos la etiqueta: tomamos solo el número inicial
  scale_x_discrete(labels = function(x) sub("\\[|\\(|,.*", "", x)) +
  labs(x = "Number of genomes per cluster (binned)", 
       y = "Number of gene clusters") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.title = element_blank())


## ojo
n_genomes <- max(cluster_hist$`No. isolates`)
bin_width <- 50

# Separar core manualmente
cluster_hist_binned <- cluster_hist %>%
  mutate(
    bin = case_when(
      `No. isolates` == n_genomes ~ "Core_bin",  # bin exclusivo para Core
      TRUE ~ as.character(cut(`No. isolates`, breaks = seq(0, n_genomes - 2, by = bin_width), right = FALSE))
    )
  ) %>%
  mutate(bin = cut(`No. isolates`, breaks = seq(0, max(`No. isolates`)+10, by = 10), right = FALSE)) %>%
  group_by(bin, Category) %>%
  summarise(total_genes = sum(n), .groups = "drop")
  
# Crear bin y agrupar
cluster_hist_binned <- cluster_hist %>%
  mutate(
    bin = case_when(
      `No. isolates` == n_genomes ~ "Core_bin",  # Bin separado para Core
      TRUE ~ as.character(cut(`No. isolates`, breaks = seq(0, n_genomes, by = bin_width), right = FALSE))
    )
  ) %>%
  group_by(bin, Category) %>%
  summarise(total_genes = sum(n), .groups = "drop") %>%
  # Reordenar niveles del bin para que Core_bin vaya al final
  mutate(bin = factor(bin, levels = c(
    sort(unique(bin[bin != "Core_bin"])),
    "Core_bin"
  )))





# Extraer solo genes Core y Soft-core gene_list
core_soft_genes_1 <- pan_categories_1 %>%
  filter(Category %in% c("Core", "Soft-core")) %>%
  pull(Gene)


##heatmap
# Usando data.table (más rápido y recomendado)
library(data.table)
library(pheatmap)
rtab_data <- fread("gene_presence_absence.Rtab")

# Convertir primera columna a rownames
mat <- as.matrix(rtab_data[, -1, with = FALSE]) 
rownames(mat) <- rtab_data[[1]]

# Hacer heatmap
pheatmap(mat, 
         color = pal(),   # 0 = blanco, 1 = negro
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         show_rownames = FALSE,         # oculta nombres de genes (suelen ser miles)
         fontsize_col = 6)   


###
pan_2 <- read_csv("gene_presence_absence_moderate.csv")



# Contar cuántos genomas tienen cada gen
pan_categories_2 <- pan_2 %>%
  mutate(No_isolates = rowSums(!is.na(.[ , 4:ncol(.)])))
# Categorizar (suponiendo que ya sabes cuántos genomas hay)
n_genomes <- length(colnames(pan_2)) - 3  # Resta columnas no-genómicas

pan_categories_group_2 <- pan_categories_2 %>%
  mutate(Category = case_when(
    No_isolates >= 0.99 * n_genomes     ~ "Core",
    No_isolates >= 0.95 * n_genomes     ~ "Soft-core",
    No_isolates >= 0.15 * n_genomes     ~ "Shell",
    TRUE                                ~ "Cloud"
  ))

pan_categories_group_2 %>% count(Category)

ids <- colnames(genes_interes_2)[4:ncol(genes_interes_2)]
writeClipboard(ids)                         # uno por línea
# o separados por coma:
writeClipboard(paste(ids, collapse = ","))


genes_interes_2 <- pan_categories_group_2 %>%
  filter(Category %in% c("Core", "Soft-core"))

row85 <- genes_interes_2[85, ]
sum(!is.na(unlist(row85)))

GEN <- "group_738"   # cambia por el gene/cluster que te interesa
row <- genes_interes[Gene==GEN, ..idx]
faltan <- names(row)[vapply(row, function(x) is_blank(x), logical(1))]
faltan  # genomas donde está vacío



# Extraer solo genes Core y Soft-core gene_list
genes_interes_2_final <- pan_categories_group_2 %>%
  filter(Category %in% c("Core", "Soft-core")) %>%
  pull(Gene)
  
writeLines(genes_interes_2_final, "gene_list_moderate.txt")


gpa <- read.csv("gene_presence_absence.csv", check.names = FALSE)
ann <- which(names(gpa) == "Annotation")
idx <- (ann+1):ncol(gpa)

GEN <- "group_2034"

row <- gpa[gpa$Gene == GEN, idx, drop = FALSE]

is_blank <- function(x) is.na(x) || trimws(x) == ""

library(data.table)

# 1) leer y preparar
gpa <- fread("gene_presence_absence.csv", check.names = FALSE)
setDT(gpa)                               # asegura data.table
ann <- match("Annotation", names(gpa))   # posición de la columna "Annotation"
stopifnot(!is.na(ann))                   # si falla, no existe esa columna
idx <- (ann + 1):ncol(gpa)               # columnas de genomas

# 2) fila del gen y columnas de genomas




gpa <- read.csv("gene_presence_absence.csv", check.names = FALSE)
ann <- match("Annotation", names(gpa)); stopifnot(!is.na(ann))
idx <- (ann + 1):ncol(gpa)

GEN <- "group_2034"
row <- gpa[gpa$Gene == GEN, idx, drop = FALSE]

# ---- GUARDAR A CSV
fwrite(data.table(gene=GEN, present_pct, n_present=length(present),
                  n_total=length(idx), missing_genomes=paste(faltan, collapse=";")),
       "resultado_un_gene.csv")











  