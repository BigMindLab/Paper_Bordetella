
setwd("C:/Users/jalme/OneDrive/Escritorio/Immunoinformatics")

install.packages("tidyverse")
library(openxlsx)
library(tidyverse)
library(dplyr)
library(tidyr)

TP <- read.xlsx("transmembrane_bordetella_TP.xlsx")
TN <- read.xlsx("transmembrane_bortedella_TN.xlsx")
transmembrane_benchmark <- read.xlsx("Benchmark_transmembrane.xlsx")

TN$UniProt <- NA

bordetella_transmembrane <- rbind(TP, TN)

identical(bordetella_transmembrane$Entry, transmembrane_benchmark$Entry) #Ver si todos están en el mismo orden y si son idénticos

##
transmembrane_benchmark <- transmembrane_benchmark %>%
  filter(!(UniProt == 0 & DeepTMHMM == 0 & TOPCONS == 0))



##Opción1
# lógica condicional para cada fila
transmembrane_benchmark$score_deep <- mapply(function(a, b) {
  if (is.na(a) || is.na(b)) {
    return(NA)  # Devuelve NA si falta algún dato
  } else if (a == b) {
    return(a / b)
  } else if (a > b) {
    return((b) / a)
  } else {
    return(((b - a) * -1 + a) / a)
  }
}, transmembrane_benchmark$UniProt, transmembrane_benchmark$DeepTMHMM)


##Opción2
transmembrane_benchmark$score_deep <- mapply(function(a, b) {
  if (is.na(a) || is.na(b)) {
    return(NA)  # Devuelve NA si falta algún dato
  } else if (a == 0 && b == 0) {
    return(1) 
  } else if (a == b) {
    return(a / b)
  } else if (a > b) {
    return((b) / a)
  } else {
    return(((b - a) * -1 + a) / a)
  }
}, transmembrane_benchmark$UniProt, transmembrane_benchmark$DeepTMHMM)


# Mostrar resultados
print(transmembrane_benchmark)

##Opción 1

transmembrane_benchmark$score_topcons2 <- mapply(function(a, b) {
  if (is.na(a) || is.na(b)) {
    return(NA)  # Devuelve NA si falta algún dato
  } else if (a == b) {
    return(a / b)
  } else if (a > b) {
    return((b) / a)
  } else {
    return(((b - a) * -1 + a) / a)
  }
}, transmembrane_benchmark$UniProt, transmembrane_benchmark$TOPCONS)



##Opción 2

transmembrane_benchmark$score_topcons2 <- mapply(function(a, b) {
  if (is.na(a) || is.na(b)) {
    return(NA)  # Devuelve NA si falta algún dato
  } else if (a == 0 && b == 0) {
    return(1)  # Caso especial: ambos son cero
  } else if (a == b) {
    return(a / b)  # Esto devolvería 1, salvo que a = b = 0, que ya está cubierto arriba
  } else if (a > b) {
    return(b / a)
  } else {
    return(((b - a) * -1 + a) / a)
  }
}, transmembrane_benchmark$UniProt, transmembrane_benchmark$TOPCONS)


##

dataset_1 <- data.frame(
  Entry = transmembrane_benchmark$Entry,
  Score = transmembrane_benchmark$score_deep,
  Software = "DeepTMHMM"
)

dataset_2 <- data.frame(
  Entry = transmembrane_benchmark$Entry,
  Score = transmembrane_benchmark$score_topcons2,
  Software = "TOPCONS2"
)

dataset_new <- rbind(dataset_1, dataset_2)
colnames(dataset_new)[colnames(dataset_new) == "Entry"] <- "ID"
dataset_new$Method <- "Presence/Absence"


length(unique(na.omit(diferencias_score)))
any(is.infinite(diferencias_score))
length(na.omit(diferencias_score))

# Diferencias entre DeepTMHMM y TOPCONS2
diferencias_score <- dataset_1$Score - dataset_2$Score

# Reemplazar -Inf con NA
diferencias_score[diferencias_score == -Inf] <- NA

# O, mejor, limpiar cualquier valor no numérico válido (Inf, -Inf, NaN)
diferencias_score[!is.finite(diferencias_score)] <- NA

# Prueba de normalidad
shapiro.test(diferencias_score)

qqnorm(diferencias_score)
qqline(diferencias_score, col = "red")


##TEST DE WILCOXON
# Filtrar por software
deep_score <- dataset_1 %>% filter(Software == "DeepTMHMM")
topcons_score <- dataset_2 %>% filter(Software == "TOPCONS2")

# Asegurar orden
deep_score <- deep_score[order(deep_score$Entry), ]
topcons_score <- topcons_score[order(topcons_score$Entry), ]

# Comparación pareada
wilcox.test(deep_score$Score, topcons_score$Score, paired = TRUE)





write.csv(dataset, "dataset.csv")


lines <- readLines("TMRs_1.gff3")

# Extraer el ID de la proteína cuando aparezca la línea "Length:"
protein_ids <- ifelse(grepl("^#\\s*sp_", lines), 
                      sub(".*#\\s*(\\S+)\\s+Length.*", "\\1", lines), 
                      NA)

# Rellenar hacia adelante para que cada línea tenga su correspondiente ID
protein_ids <- zoo::na.locf(protein_ids)

# Extraer solo las líneas TMhelix
tm_lines <- grep("TMhelix|Beta sheet", lines, value = TRUE)

# Asociar el ID correspondiente a cada línea TMhelix
tm_ids <- protein_ids[grep("TMhelix|Beta sheet", lines)]

# Crear data frame con inicio y fin
df_tm <- tibble(line = tm_lines, ID = tm_ids) %>%
  separate(line, into = c("ID2", "Feature", "Start", "End"), sep = "\\s+", extra = "drop") %>%
  mutate(Start = as.numeric(Start),
         End = as.numeric(End),
         Length = End - Start + 1) %>%
  select(ID, Start, End, Length)

df_tm <- df_tm %>%
  mutate(ID = gsub("^#\\s*sp_|_.*$", "", ID))

# Crear un vector con el ID actual para cada línea
current_id <- NA
id_vector <- character(length(lines))

for (i in seq_along(lines)) {
  line <- lines[i]
  if (grepl("^#\\s*sp_", line)) {
    current_id <- sub(".*#\\s*(\\S+)\\s+Length.*", "\\1", line)
  }
  id_vector[i] <- current_id
}

# Extraer solo las líneas con TMhelix
tm_lines <- lines[grepl("TMhelix", lines)]
tm_ids <- id_vector[grepl("TMhelix", lines)]

# Armar el dataframe
df_tm <- tibble(line = tm_lines, ID = tm_ids) %>%
  separate(line, into = c("ID2", "Feature", "Start", "End"), sep = "\\s+", extra = "drop") %>%
  mutate(Start = as.numeric(Start),
         End = as.numeric(End),
         Length = End - Start + 1) %>%
  select(ID, Start, End, Length)



root_dir <- "rst_7ngr9w3l/rst_7ngr9w3l"

# Encuentra todas las carpetas seq_*
seq_dirs <- list.dirs(root_dir, recursive = FALSE, full.names = TRUE)
seq_dirs <- seq_dirs[grepl("seq_\\d+$", seq_dirs)]

# Inicializa resultado
all_tm_results <- list()

for (dir in seq_dirs) {
  file_path <- file.path(dir, "query.result.txt")
  if (!file.exists(file_path)) next
  
  lines <- readLines(file_path)
  top_line_index <- grep("TOPCONS predicted topology:", lines)
  if (length(top_line_index) == 0) next
  
  topology <- lines[top_line_index + 1]  # línea después de "predicted topology"
  topology_chars <- unlist(strsplit(topology, split = ""))
  
  id <- basename(dir)
  tm_positions <- which(topology_chars == "M")
  
  if (length(tm_positions) == 0) next
  
  # Detectar cambios de bloques
  diffs <- c(Inf, diff(tm_positions))
  start_idx <- which(diffs != 1)
  
  starts <- tm_positions[start_idx]
  ends <- c(tm_positions[start_idx[-1] - 1], tm_positions[length(tm_positions)])
  lengths <- ends - starts + 1
  
  df <- tibble(ID = id, Start = starts, End = ends, Length = lengths)
  all_tm_results[[id]] <- df
}

# Unir todo en un solo data frame
final_df <- bind_rows(all_tm_results)




file.exists("rst_7ngr9w3l/rst_7ngr9w3l/seq_0/query.result.txt")

list.dirs("rst_7ngr9w3l/rst_7ngr9w3l", recursive = FALSE)




plot <- ggplot(data = dataset, aes(x = dataset$Software, y = dataset$Score)) +
  theme_bw() + facet_wrap(~dataset, nrow = 1, scales = "free_y") +
  geom_violin(scale = "count", alpha = 0.05, fill = "yellow",
              show.legend = FALSE, trim = TRUE, linewidth = 0.4) +
  geom_boxplot(width = 0.15, fill = "grey90", outlier.shape = NA) +
  geom_errorbar(stat = "boxplot", width = 0.4, linewidth = 1, color = "red",
                aes(ymin = after_stat(middle), ymax = after_stat(middle))) +
  geom_point(aes(fill = type), shape = 21, size = 3, alpha = 1,
             show.legend = FALSE, position = position_jitter(width = 0.12)) +
  scale_fill_manual(values = c("Control" = "#b9ffffff", "TB" = "#fab8c2ff")) +
  scale_y_continuous(breaks = function(x) seq(floor(min(x, na.rm = TRUE)),
                                              ceiling(max(x, na.rm = TRUE)), by = 1)) +
  labs(title = "Benchmark transmembrane topology", y = expression("Value range") +

  theme(plot.title = element_text(size = 22),
        strip.text = element_text(size = 17),
        strip.background = element_blank(),
        axis.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18))

  df.signif$ypos <- NA

  for (batch in unique_batches) {
    batch_rows <- df.signif$dataset == batch
    max_value <- max(df.counts$log2_gc[is.finite(df.counts$log2_gc) & df.counts$dataset == batch], na.rm = TRUE)
    prop <- 0.08*diff(range(df.counts$log2_gc[is.finite(df.counts$log2_gc) & df.counts$dataset == batch], na.rm = TRUE))

    df.signif$ypos[batch_rows] <- max_value + prop
  }

  df.signif <- df.signif[df.signif$labels != "ns", ]
  plot <- plot + geom_signif(data = df.signif, manual = TRUE, size = 1, colour = "black",
                             aes(xmin = group1, xmax = group2, annotations = labels, y_position = ypos),
                             textsize = 8, vjust = 0.5, tip_length = 0.02)

  ggsave(paste0("VSB/", genes[i], ".jpg"), plot = plot, width = 6, height = 7, dpi = 600)
}


signalpeptide_benchmark <- read.xlsx("Signal_Peptide.xlsx")



confusion_matrix <- function(df, col_a = "UniProt", col_b = "TOPCONS") {
  # Inicializar contadores
  true_positive <- 0
  false_negative <- 0
  false_positive <- 0
  true_negative <- 0

  # Recorrer filas
  for (i in 1:nrow(df)) {
    a <- as.numeric(df[[col_a]][i])
    b <- as.numeric(df[[col_b]][i])

    if (!is.na(a) && !is.na(b)) {
      if (a == 1 && b == 1) {
        true_positive <- true_positive + 1
      } else if (a == 1 && b == 0) {
        false_negative <- false_negative + 1
      } else if (a == 0 && b == 1) {
        false_positive <- false_positive + 1
      } else if (a == 0 && b == 0) {
        true_negative <- true_negative + 1
      }
    }
  }

  return(data.frame(
    True_Positive = true_positive,
    False_Negative = false_negative,
    False_Positive = false_positive,
    True_Negative = true_negative
  ))
}


resultado <- confusion_matrix(df = signalpeptide_benchmark, col_a = "UniProt", col_b = "TOPCONS")

print(resultado)






