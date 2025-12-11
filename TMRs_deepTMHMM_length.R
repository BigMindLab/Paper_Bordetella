
library(dplyr)
library(tidyr)
library(zoo)

lines <- readLines("TMRs_1.gff3")

# ID de proteína
protein_ids <- ifelse(grepl("^#\\s*sp_", lines), 
                      sub(".*#\\s*(\\S+)\\s+Length.*", "\\1", lines), 
                      NA)
protein_ids <- zoo::na.locf(protein_ids)

# Filtrar TMhelix o Beta sheet
tm_idx <- grep("TMhelix|Beta sheet", lines)
tm_lines <- lines[tm_idx]
tm_ids <- protein_ids[tm_idx]

# Data frame
df_tm <- tibble(line = tm_lines, ID = tm_ids) %>%
  mutate(
    # separa en tokens y toma los dos últimos como Start y End
    Start = as.numeric(sapply(strsplit(line, "\\s+"), function(x) x[length(x)-1])),
    End   = as.numeric(sapply(strsplit(line, "\\s+"), function(x) x[length(x)])),
    Length = End - Start + 1
  ) %>%
  select(ID, Start, End, Length)

df_tm
df_tm <- df_tm %>%
  mutate(ID = gsub("^#\\s*sp_|_.*$", "", ID))
