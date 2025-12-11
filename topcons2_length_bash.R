
setwd("C:/Users/jalme/OneDrive/Escritorio/Immunoinformatics")
library(openxlsx)
library(tidyverse)
library(dplyr)
library(tidyr)

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

write_tsv(final_df, "TOPCONS2_length.tsv")

topcons <- read_tsv("TOPCONS2_length.tsv")


