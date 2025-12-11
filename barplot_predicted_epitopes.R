
setwd("C:/Users/jalme/OneDrive/Escritorio/Immunoinformatics")

# Paquetes
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)
library(RColorBrewer)

# 1) Leer y preparar datos
path <- "Supplementary_1_All_epitopes_from_49_proteins.xlsx"
df <- read_excel(path)
df_HLAI <- df[df$Type == "HLA-I", ]
df <- df_HLAI
df_HLAII <- df[df$Type == "HLA-II", ]
df <- df_HLAII

# Asegura que la columna Protein exista y las columnas de países sean numéricas
paises <- c("ARG","BOL","BRA","CHI","COL","ECU","PAR","PER","VEN")

df2 <- df %>%
  mutate(across(all_of(paises), as.numeric)) %>%
  group_by(Protein) %>%
  summarise(across(all_of(paises), ~sum(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(Total = rowSums(across(all_of(paises))))


# Formato largo (Protein, Pais, Frecuencia)
df_long <- df2 %>%
  select(-Total) %>%
  pivot_longer(cols = all_of(paises), names_to = "Country",
               values_to = "Frequency")

# Ordenar proteínas por total (útil para todos los gráficos)
df_long <- df_long %>%
  left_join(df2 |> select(Protein, Total), by = "Protein") %>%
  mutate(Protein = fct_reorder(Protein, Total))

df_long$Protein <- str_replace(df_long$Protein, "^(.)", toupper)

# Reordenar por frecuencia total
df_long <- df_long %>%
  group_by(Protein) %>%
  mutate(Total = sum(Frequency, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Protein = fct_reorder(Protein, Total, .fun = sum, .desc = FALSE))

ggplot(df_long, aes(x = Protein, y = Frequency, fill = Country)) +
  geom_col() +
  coord_flip() +
  labs(x = "Proteins", y = "Epitope Frequency",
       title = "") +
  theme_minimal() +
  theme(legend.position = "right") + theme_bw()

ggplot(df_long, aes(x = Protein, y = Frequency, fill = Country)) +
  geom_col() +
  coord_flip() +
  labs(x = "Proteins", y = "Epitope Frequency") +
  theme_bw() +
  scale_fill_manual(values = c(
    "ARG" = "#08306B",
    "BOL" = "#08519C",
    "BRA" = "#2171B5",
    "CHI" = "#4292C6",
    "COL" = "#6BAED6",
    "ECU" = "#9ECAE1",
    "PAR" = "#C6DBEF",
    "PER" = "#DEEBF7",
    "VEN" = "#045A8D"
  ))

ggplot(df_long, aes(x = Protein, y = Frequency, fill = Country)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  scale_fill_brewer(option = "Set2")

###
# Cargar RColorBrewer
library(RColorBrewer)

# Orden sugerido por total (opcional)
df_long <- df_long %>%
  group_by(Protein) %>% mutate(Total = sum(Frequency, na.rm=TRUE)) %>% ungroup() %>%
  mutate(Protein = fct_reorder(Protein, Total))

# Paleta pastel con buen contraste (9 países)
pal <- c("ARG"="lightblue","BOL"="#FDBF6F","BRA"="#B2DF8A","CHI"="#CAB2D6",
         "COL"="#FFDD89","ECU"="#CCEBC5","PAR"="#BC80BD","PER"="#FB9A99","VEN"="#8DD3C7")
pal <- c(
  ARG = "#B7D3F2",  # cielo suave
  BOL = "#FFC8A2",  # melocotón
  BRA = "#C6E5B3",  # verde oliva claro
  CHI = "#E2C5F1",  # lila cálido
  COL = "#FFE7A1",  # mantequilla cálida
  ECU = "#D6EEEB",  # aqua tibio (muy suave)
  PAR = "#F0C3CF",  # rosa-lila cálido
  PER = "#FFB6C1",  # rosa bebé
  VEN = "#E1F0C4"   # azul glaciar   # pistacho claro   # verde agua brillante
)



ggplot(df_long, aes(x = Protein, y = Frequency, fill = Country)) +
  geom_col(width = 0.85, color = "white") +   # sin alpha; borde suave para separar capas
  coord_flip() +
  scale_fill_manual(values = pal_named, name = "Country", drop = FALSE) +  # LEYENDA
  labs(x = "Protein", y = "Epitope Frequency", title = NULL) +
  theme_bw( ) +
  theme(
    legend.position = "right",
    legend.key.height = unit(5, "mm"),
    panel.grid.major.y = element_blank() # limpia las líneas detrás de las barras
  )

###
# Defino un set base de colores suaves
base_cols <- c("#A1C9F4", "#FFB6C1", "#FFE29A", "#BDE7E1")  # celeste, rosa bebé, amarillo suave, aqua

# Genero la función de paleta
pal_fun <- colorRampPalette(base_cols)

# Saco 9 colores (uno por país)
pal_crp <- pal_fun(9)

# Asigno a los países
pal_crp_named <- setNames(pal_crp, c("ARG","BOL","BRA","CHI","COL","ECU","PAR","PER","VEN"))

df_long_I <- df_long

# Uso en ggplot
ggplot(df_long, aes(x = Protein, y = Frequency, fill = Country)) +
  geom_col(width = 0.85, color = "white") +
  coord_flip() +
  scale_fill_manual(values = pal_named, name = "Country") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Protein", y = "HLA-I Epitope Frequency ") +
  theme_bw(base_size = 14 ) + 
  theme(
    legend.position = c(0.9, 0.48), # coordenadas relativas (x, y)
    legend.background = element_rect(fill = alpha("white", 1), color = "black"), # fondo semitransparente
    legend.key = element_rect(fill = NA),
    axis.title.y = element_text(size = 16, margin = margin(r = 1)),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    strip.text = element_text(size = 14)# sin fondo en cada clave
  )


#pasteles un poquito más saturados
# colores base pastel pero más saturados
base_cols <- c("#8DBDF0", "#FF9EB5", "#FFD480", "#9EE0D5")

# generar la paleta con más vida
pal_fun <- colorRampPalette(base_cols)
pal_stronger <- pal_fun(9)  # o 9 según países

# asignar a tus países
pal_named <- setNames(pal_stronger, c("ARG","BOL","BRA","CHI","COL","ECU","PAR","PER","VEN"))


####
ggplot(df_long, aes(x = Pais, y = Protein, fill = Frecuencia)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "País", y = "Proteína",
       title = "Mapa de calor: frecuencia de epítopos por proteína y país",
       fill = "Frecuencia") +
  theme_minimal() +
  theme(panel.grid = element_blank())

###
ggplot(df_long, aes(x = Protein, y = Frecuencia)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  facet_wrap(~ Pais, ncol = 3, scales = "free_y") +
  labs(x = "Proteína", y = "Frecuencia de epítopos",
       title = "Epítopos por proteína, facetado por país") +
  theme_minimal()



##heatmap
ggplot(df_long, aes(x = Pais, y = Protein, fill = Frecuencia)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "País", y = "Proteína",
       title = "Mapa de calor: frecuencia de epítopos por proteína y país",
       fill = "Frecuencia") +
  theme_minimal() +
  theme(panel.grid = element_blank())


##facetas
ggplot(df_long, aes(x = Protein, y = Frecuencia)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  facet_wrap(~ Pais, ncol = 3, scales = "free_y") +
  labs(x = "Proteína", y = "Frecuencia de epítopos",
       title = "Epítopos por proteína, facetado por país") +
  theme_minimal()



##100%proporciones 
df_prop <- df_long %>%
  group_by(Protein) %>%
  mutate(Prop = Frecuencia / sum(Frecuencia, na.rm = TRUE)) %>%
  ungroup()

ggplot(df_prop, aes(x = Protein, y = Prop, fill = Pais)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Proteína", y = "Proporción de epítopos",
       title = "Distribución proporcional de epítopos por proteína y país") +
  theme_minimal()
