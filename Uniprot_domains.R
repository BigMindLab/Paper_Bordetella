
library(readxl)
library(tidyverse)

# Cargar archivo Excel (ajusta el nombre si es necesario)
df <- read_xlsx("transmembrane_bordetella_TP.xlsx")

# Asegúrate de que las columnas relevantes estén correctamente nombradas
# Ejemplo: Entry = ID, Transmembrane = la columna con los segmentos

# Extraer segmentos transmembrana por fila
Uniprot_domains <- df %>%
  filter(!is.na(Transmembrane)) %>%
  rowwise() %>%
  mutate(matches = str_extract_all(Transmembrane, "TRANSMEM\\s+(\\d+)\\.\\.(\\d+)")) %>%
  unnest(matches) %>%
  mutate(
    Start = as.integer(str_extract(matches, "\\d+")),
    End = as.integer(str_extract(matches, "(?<=\\.\\.)(\\d+)")),
    Length = End - Start + 1
  ) %>%
  select(ID = Entry, Start, End, Length)

# Mostrar resultados
print(Uniprot_domains)

# Guardar como TSV si deseas
# write_tsv(result, "uniprot_tm_segments.tsv")

ids <- c(
  "A0A0H3LKL4","J7QLC0","Q07703","Q9NPF7","P04978","Q45340","P0A3I5","P35077",
  "Q7WF17","Q7WFT9","Q7WD07","Q7WGU8","Q70I53","P0DP27","Q775D6","Q45374",
  "G3XD94","G3XD23","Q9HZ76","P04981","Q7VYQ1","P04979","Q7VSX7","Q7VTJ9",
  "A9I0E6","P0A3R5","P0A4H2","P81549","Q04064","Q7VSX9","Q7W0D3","Q7VSX3",
  "Q7VSX5","P0A323","Q7VSX6","A0A0H3LQK8","Q7WLV0","Q7WLU9","J7RC67","Q48252",
  "Q8L1U3","Q8FXK7","Q9RPX5","Q7CEG3","Q7W980","Q7W2Q0","P67675","Q7W0K2",
  "Q7VV85","P67673","P12255","Q7VU58","A0A0H3LX82","P0A321","A0A0H3LKV3",
  "Q06064","Q9RPX9","Q05115","Q7VU61","Q7WKU6","B3EWP1"
)

# Supongamos que tu data frame se llama df y la columna con IDs es "Accession"
Uniprot_filtrado <- Uniprot_domains[Uniprot_domains$ID %in% ids, ]



