

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ---- las 22 proteínas que quieres graficar ----
keep <- c("bcrD","BrkA","BscC","BscJ","CycB","FhaC","Fim3","FimD","FimX",
          "FlgA","FlgD","FlgE","FlgI","FlgL","FliD","FliF","MotB","PtlE",
          "PtlF","Vag8","WlbG","WlbL")

# -------- 1) Mapa de longitudes por proteína --------
# Usa S2 y, si no existe 'Protein_length', toma la mayor coordenada 'End' como proxy.
# Calcula longitud de proteína como el máximo 'End' observado
make_len_map <- function(df) {
  df %>%
    mutate(
      End = suppressWarnings(as.numeric(End)),
      Start = suppressWarnings(as.numeric(Start))
    ) %>%
    group_by(Protein) %>%
    summarise(prot_len = max(End, na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(prot_len))
}
len_map <- make_len_map(S3)

# (opcional) comprueba si falta alguna de las 22 en el mapa de longitudes
setdiff(keep, len_map$Protein)  # debería devolver character(0)

# -------- 2) Función para armar el scatter --------
plot_len_vs_epitopes <- function(S3, class_to_plot = "HLA-II", keep, len_map) {
  S3 %>%
    mutate(
      Type2 = toupper(gsub("_", "-", as.character(Type))), # normaliza "HLA_II" -> "HLA-II"
      Start = suppressWarnings(as.numeric(Start)),
      End   = suppressWarnings(as.numeric(End)),
      uid   = paste(Peptide, Start, End, sep = "|"), # ID único de epítopo
      Protein = str_to_title(Protein),
    ) %>%
    filter(Type2 == class_to_plot, Protein %in% str_to_title(keep)) %>%
    count(Protein, name = "n_epitopes") %>%
    inner_join(len_map %>% mutate(Protein = str_to_title(Protein)), by = "Protein") %>%
    ggplot(aes(x = prot_len, y = n_epitopes)) +
    geom_point(aes(fill = n_epitopes), shape = 21, color = "black", size = 4, stroke = 1) +
    geom_text_repel(aes(label = Protein), size = 3,
                    max.overlaps = Inf, box.padding = .35, point.padding = .5,
                    min.segment.length = 0,
                    vjust = -0.5,         # 🔹 Ajusta posición vertical
                    hjust = 0.5, 
                    seed = 1, nudge_y = 0.9) +
    labs(x = "Protein Length",
         y = paste("Number of", class_to_plot, "epitopes")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16))
}

# -------- 3) Llamadas (elige fuente de epítopos) --------
# A) Como en tu figura “general”: usa S2 (post filtro humano)
p_HLAII_S2 <- plot_len_vs_epitopes(S2, "HLA-II", keep, len_map)
# B) Si quieres solo CS≥7: usa S3
p_HLAII_S3 <- plot_len_vs_epitopes(S3, "HLA-II", keep, len_map)
p_HLAII_S3

# Panel comparativo A/B con S2 vs S3
p_HLAII_S2 + p_HLAII_S3 + plot_layout(ncol = 2)




#####EXTRA
library(dplyr)
library(stringr)
library(ggplot2)
library(purrr)

# 22 proteínas
keep <- c("bcrD","brkA","bscC","bscJ","cycB","fhaC","fim3","fimD","fimX",
          "flgA","flgD","flgE","flgI","flgL","fliD","fliF","motB","ptlE",
          "ptlF","vag8","wlbG","wlbL")

# --- 0) Chequeo: ahora sí deben calzar en S3$Protein ---
length(intersect(unique(S4$Protein), keep))
# debería ser > 0

# --- helpers para “empaquetar” segmentos y marcar solape ---
pack_lanes <- function(df){
  df <- df[order(df$Start, df$End), , drop = FALSE]
  ends <- numeric(0); lane <- integer(nrow(df))
  for(i in seq_len(nrow(df))){
    placed <- FALSE
    for(k in seq_along(ends)){
      if(df$Start[i] > ends[k]) { lane[i] <- k; ends[k] <- df$End[i]; placed <- TRUE; break }
    }
    if(!placed){ ends <- c(ends, df$End[i]); lane[i] <- length(ends) }
  }
  df$lane <- lane; df
}
mark_overlaps <- function(df){
  n <- nrow(df); ov <- rep(FALSE, n)
  if(n >= 2){
    for(i in 1:(n-1)) for(j in (i+1):n)
      if(df$Start[j] <= df$End[i] && df$End[j] >= df$Start[i]) { ov[i] <- TRUE; ov[j] <- TRUE }
  }
  df$overlap <- ov; df
}

# --- 1) Construir dataset de posiciones (solo HLA-II y tus 22 proteínas) ---
epi_pos <- S4 %>%
  mutate(
    Type2 = toupper(gsub("_","-", as.character(Type))),  # normaliza por si acaso
    Start = suppressWarnings(as.numeric(Start)),
    End   = suppressWarnings(as.numeric(End))
  ) %>%
  filter(Type2 == "HLA-II", Protein %in% keep, is.finite(Start), is.finite(End)) %>%
  arrange(Protein, Start, End) %>%
  group_by(Protein) %>%
  group_modify(~ .x %>% pack_lanes() %>% mark_overlaps()) %>%
  ungroup()

# --- (opcional) longitudes desde S2 para escalar ejes ---
len_map <- S3 %>%
  mutate(End = suppressWarnings(as.numeric(End))) %>%
  group_by(Protein) %>% summarise(prot_len = max(End, na.rm = TRUE), .groups="drop")

epi_pos <- epi_pos %>% left_join(len_map, by = "Protein")

# --- 2) Plot: Start–End por proteína; rojo = solapados ---
ggplot(epi_pos, aes(x = Start, xend = End, y = lane, yend = lane)) +
  geom_segment(aes(color = overlap), linewidth = 2, lineend = "round") +
  scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#E63946"),
                     labels = c("Unique", "Overlapping"),
                     name   = "Epitope") +
  facet_wrap(~ Protein, scales = "free_x") +
  labs(x = "Position (aa)", y = "Lane",
       subtitle = "HLA-II epitopes; red segments indicate overlaps") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10))

##
OJO
# ---- Lista de 22 proteínas de interés ----
keep <- c("BcrD","BrkA","BscC","BscJ","CycB","FhaC","Fim3","FimD","FimX",
          "FlgA","FlgD","FlgE","FlgI","FlgL","FliD","FliF","MotB","PtlE",
          "PtlF","Vag8","WlbG","WlbL")

# ---- Conteo de epítopos únicos por proteína ----
summary_S4 <- S4 %>%
  filter(Protein %in% keep) %>%
  distinct(Protein, Peptide) %>%       # evita contar duplicados
  count(Protein, name = "n_epitopes") %>%
  right_join(tibble(Protein = keep), by = "Protein") %>% # incluye las que no tengan epitopos
  replace_na(list(n_epitopes = 0)) %>%
  arrange(desc(n_epitopes))

# ---- Resultado ----
print(summary_S4)




plot_len_vs_epitopes <- function(S3, class_to_plot = "HLA-II", keep, len_map) {
  S3 %>%
    mutate(
      Type2 = toupper(gsub("_", "-", as.character(Type))), # Normaliza HLA_II -> HLA-II
      Start = suppressWarnings(as.numeric(Start)),
      End   = suppressWarnings(as.numeric(End)),
      uid   = paste(Peptide, Start, End, sep = "|"),       # ID único de epítopo
    ) %>%
    filter(Type2 == class_to_plot, Protein %in% keep) %>%
    count(Protein, name = "n_epitopes") %>%
    inner_join(len_map, by = "Protein") %>%
    ggplot(aes(x = prot_len, y = n_epitopes)) +
    geom_point(shape = 21, fill = "skyblue", color = "black", size = 3, stroke = 1) + # Un solo color
    geom_text_repel(
      aes(label = Protein), size = 4,
      max.overlaps = Inf, box.padding = .35, point.padding = .5,
      min.segment.length = 0, 
      direction = "y",   # 🔹 mueve solo en vertical
      force_pull = 0, nudge_y = 0.3, seed = 1
    ) +
    labs(x = "Protein Length",
         y = paste("Number of", class_to_plot, "epitopes")) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(size = 16),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16))
}


S3$Protein <- paste0(
  toupper(substr(S3$Protein, 1, 1)),   # primera letra en mayúscula
  substr(S3$Protein, 2, nchar(S3$Protein))  # resto igualito
)

p <- plot_len_vs_epitopes(S3, class_to_plot = "HLA-II", keep, len_map)
p

unique(df$Protein)
unique(len_map$Protein)
keep






