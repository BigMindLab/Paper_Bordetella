##Gráficos adicionales epitope prediction 

# ====== Librerías ======
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

# ====== Archivos (ajusta si hace falta) ======
f_s1 <- "Supplementary_1_All_epitopes_from_49_proteins.xlsx"
f_s2 <- "Suplementary_2_ALL_epitopes_conserved_without_100hit_humans_from_49_protteins.xlsx"
f_s3 <- "Suplementary_3_Epitopes_score_morethan7.xlsx"
f_s4 <- "Supplementary_4_I_INSIDE_IIMAYOR7.xlsx"
f_s5 <- "Suplementary_5_Epitope_with_some_reference_iedb_from_49.xlsx"

# ====== Lee primera hoja de cada archivo (sin renombrar columnas) ======
S1 <- read_excel(f_s1, sheet = 1)
S2 <- read_excel(f_s2, sheet = 1)
S3 <- read_excel(f_s3, sheet = 1)
S4 <- read_excel(f_s4, sheet = 1)
S5 <- read_excel(f_s5, sheet = 1)

# Países (tal como están nombrados en tus Excels)
country_cols <- c("ARG","BOL","BRA","CHI","COL","ECU","PAR","PER","VEN")

# Asegurar orden de clases
lvl_class <- c("HLA-I","HLA-II")
norm_type <- function(x) factor(as.character(x), levels = lvl_class)

# =========================================================
# FIGURA 1 (A+B en una sola): Conteos globales + Histograma CS
# =========================================================

# 1A) Conteos globales por clase: crudo (S1) vs filtrado humano (S2)
# -------- Fig 1A --------
counts_raw <- S1 |> count(Type, name = "n") |> mutate(Etapa = "Raw")
counts_flt <- S2 |> count(Type, name = "n") |> mutate(Etapa = "Filtered")
counts <- bind_rows(counts_raw, counts_flt) |>
  mutate(
    Type  = norm_type(Type),
    Etapa = factor(Etapa, levels = c("Raw", "Filtered"))  # Raw primero
  )

pal_etapa <- c(          # Para p1 (Etapa)
  "Raw"      = "#5FD9C3",  # menta pastel viva
  "Filtered" = "#FF76B4"   # rosa pastel encendido
)

pal_type <- c(            # Para p2 (Type)
  "HLA-I"  = "#79ADFF",   # azul pastel vivo
  "HLA-II" = "#FFBE78"    # durazno pastel vivo
)

p1 <- ggplot(counts, aes(Type, n, fill = Etapa)) +
  geom_col(position = position_dodge(width = .7), width = .65) +
  geom_text(aes(label = comma(n)),
            position = position_dodge(width = .7), vjust = -0.25, size = 3.5) +
  scale_fill_manual(values = pal_etapa, name = "Type") +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  scale_x_discrete(expand = expansion(mult = c(0.4, 0.4))) +
  labs(x = NULL, y = "Number of epitopes") +
  theme_bw() +
  theme(legend.position = "top",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

p1

# 1B) Distribución de Score (CS) con corte CS≥7 (desde S1)
p2 <- S1 |>
  mutate(Type = norm_type(Type)) |>
  ggplot(aes(x = Score, fill = Type)) +
  geom_histogram(binwidth = 0.5, alpha = .7, position = "identity") +
  geom_vline(xintercept = 7, linetype = 2) +
  scale_fill_manual(values = pal_type, name = " ") +
  scale_x_continuous(breaks = seq(0, 8, 2), limits = c(0, 8)) +
  labs(x = "Coverage Score (CS)", y = "Frequency") +
  theme_bw() +
  theme(legend.position = "top",
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

Fig1 <- p1 + p2 + plot_layout(ncol = 2)
print(Fig1)

# =========================================================
# FIGURA 2: Barras por proteína (post filtro humano; S2)
# =========================================================
S2_prot <- S2 |>
  count(Protein, Type, name = "n") |>
  group_by(Protein) |>
  mutate(total = sum(n)) |>
  ungroup()

Fig2 <- ggplot(S2_prot, aes(x = reorder(Protein, total), y = n, fill = Type)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Proteína", y = "Epítopos predichos (post filtro humano)") +
  theme_bw() +
  theme(legend.position = "top")
print(Fig2)

# =========================================================
# FIGURA 3: Heatmap CS≥7 por proteína × país (S3)
#        (usa columnas ARG..VEN con 0/1; suma por celda)
# =========================================================
stopifnot(all(country_cols %in% names(S3)))
S3_long <- S3 |>
  pivot_longer(cols = all_of(country_cols), names_to = "Country", values_to = "Presence") |>
  filter(Presence > 0) |>
  count(Protein, Country, name = "n")

Fig3 <- ggplot(S3_long, aes(Country, Protein, fill = n)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Conteo") +
  labs(x = "País", y = "Proteína") +
  theme_bw()
print(Fig3)

# =========================================================
# FIGURA 4: HLA-I anidados por epítopo HLA-II (S4)
#        (cuenta filas por HLA_II_NAME)
# =========================================================
stopifnot("HLA_II_NAME" %in% names(S4))
nested <- S4 |> count(HLA_II_NAME, name = "nested_count")

Fig4 <- ggplot(nested, aes(nested_count)) +
  geom_histogram(binwidth = 1) +
  labs(x = "HLA-I anidados por epítopo HLA-II", y = "Epítopos HLA-II (n)") +
  theme_bw()
print(Fig4)

# =========================================================
# FIGURA 5: Evidencia IEDB por proteína (S5)
#        (separa por CS≥7 vs CS<7 usando Score)
# =========================================================
stopifnot("iedb" %in% names(S5))
S5_bar <- S5 |>
  mutate(CS_tier = ifelse(Score >= 7, "CS≥7", "CS<7")) |>
  count(Protein, CS_tier, name = "n")

Fig5 <- ggplot(S5_bar, aes(x = reorder(Protein, n), y = n, fill = CS_tier)) +
  geom_col(width=.7) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Proteína", y = "Epítopos con referencia IEDB", fill = "CS") +
  theme_bw()
print(Fig5)

# =========================================================
# FIGURA 6: Funnel del pipeline (valores de tus tablas)
# =========================================================
n_total <- nrow(S1)      # 20,159
n_post  <- nrow(S2)      # 11,489 + 6,522 = 18,006 (desde S2)
n_cs7   <- nrow(S3)      # 166 (S3)
n_sel   <- 88            # seleccionados por criterios estructurales

flow <- tibble::tibble(
  step = c("Predichos (I+II)", "Filtro identidad humana", "CS ≥ 7 (HLA-II)", "Seleccionados (estructural)"),
  n    = c(n_total, n_post, n_cs7, n_sel)
)

Fig6 <- ggplot(flow, aes(x = reorder(step, -n), y = n)) +
  geom_col(width=.6, fill="grey60") +
  geom_text(aes(label = comma(n)), vjust = -0.3) +
  labs(x = NULL, y = "Conteo") +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0,.1))) +
  coord_flip()
print(Fig6)

# ====== Guardado opcional ======
ggsave("Fig1_counts_CS.png", Fig1, width = 11, height = 5, dpi = 300)
ggsave("Fig2_by_protein_posthuman.png", Fig2, width = 8, height = 10, dpi = 300)
ggsave("Fig3_CS7_heatmap_country.png", Fig3, width = 10, height = 10, dpi = 300)
ggsave("Fig4_nested_hist.png", Fig4, width = 7, height = 5, dpi = 300)
ggsave("Fig5_IEDB_bars.png", Fig5, width = 7, height = 5, dpi = 300)
ggsave("Fig6_funnel.png", Fig6, width = 6, height = 4, dpi = 300)




###extra
# --- Conteos (exactos con tus columnas) ---
n_total <- nrow(S1)
n_post  <- nrow(S2)

# Validados: si el archivo S5 son solo los validados, usa su número de filas.
# (Opcional: descomenta una de estas dos líneas si prefieres deduplicar)
# n_valid <- n_distinct(S5$Peptide)                      # por péptido único
# n_valid <- n_distinct(S5[, c("Protein","Peptide")])    # combinación proteína-péptido
n_valid <- nrow(S5)

summary_tbl <- tibble(
  Etapa = factor(c("Total predichos (49 proteínas)","Post filtro humano","Con evidencia IEDB"),
                 levels = c("Total predichos (49 proteínas)","Post filtro humano","Con evidencia IEDB")),
  n     = c(n_total, n_post, n_valid)
) %>%
  mutate(pct_total = n / n_total,
         label = paste0(comma(n), " (", percent(pct_total, accuracy = 0.1), ")"))

# --- Opción A: Funnel simple (3 barras) ---
p_funnel <- ggplot(summary_tbl, aes(x = Etapa, y = n)) +
  geom_col(width = .6, fill = "grey60") +
  geom_text(aes(label = label), vjust = -0.3, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0,.08))) +
  labs(x = NULL, y = "Número de epítopos") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

# --- Opción B: Dumbbell/slope (misma info, más compacta) ---
p_dumb <- ggplot(summary_tbl, aes(x = Etapa, y = n, group = 1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_text(aes(label = label), nudge_y = max(summary_tbl$n)*0.05, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0,.15))) +
  labs(x = NULL, y = "Número de epítopos") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

# Muestra ambos uno al lado del otro (elige el que prefieras)
p_funnel + p_dumb + plot_layout(ncol = 2)
