
filter_pangenome <- read_excel("pangenome.xlsx")
df_long_pangenome <- pivot_longer(filter_pangenome, cols = c(Total_Genes, Pangenome),
                        names_to = "Condition", values_to = "Value")

ggplot(df_long_pangenome, aes(x = Category, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Value),
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3) +
  scale_y_continuous(limits = c(0, 2700), expand = c(0,0)) +
  theme_bw(base_size = 14) +
  labs(y = "Number of Genes", x = " ") +
  scale_color_manual(values =  c("Total_Genes" = "skyblue", 
                                "Pangenome" = "navy")) +
  theme(
    legend.position = c(0.85, 0.8), # coordenadas relativas (x, y)
    legend.background = element_rect(fill = alpha("white", 1), color = "black"), # fondo semitransparente
    legend.key = element_rect(fill = NA) # sin fondo en cada clave
  )

ggplot(df_long_pangenome, aes(x = Category, y = Value, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Value),
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3) +
  scale_y_continuous(limits = c(0, 2700), expand = c(0,0)) +
  theme_bw(base_size = 14) +
  labs(y = "Number of Genes", x = " ") +
  scale_fill_manual(values =  c("Total_Genes" = "#A7D8F4", 
                                "Pangenome" = "#FFD6A5"), 
                    labels = c("Total_Genes" = "Pre-filtered Genes",
                               "Pangenome"   = "Pangenome Genes")) +
  theme(
    legend.position = c(0.78, 0.8), # coordenadas relativas (x, y)
    legend.background = element_rect(fill = alpha("white", 1), color = "black"), 
    legend.key = element_rect(fill = NA),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

#####
df_long_pangenome$Condition <- factor(
  df_long_pangenome$Condition,
  levels = c("Total_Genes", "Pangenome")  # primero Pre-filtered, luego Pangenome
)


#Gráfico de barras apiladas
ggplot(df_long_pangenome, aes(x = Condition, y = Value, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Value),
            position = position_stack(vjust = 0.5), size = 3, color = "black") +
  scale_fill_manual(values = c("Cloud"     = "#A7D8F4",  # azul
                               "Shell"     = "#FFD6A5",  # naranja
                               "Soft-core" = "#B2DF8A",  # verde
                               "Core"      = "#D7B9F7")) + # rojo
  scale_x_discrete(limits = c("Total_Genes", "Pangenome"),
                   labels = c("Total_Genes" = "Pre-filtered Genes",
                              "Pangenome" = "Pangenome Genes")) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
              axis.title.y = element_text(size = 16),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              axis.title.x = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)) +
  labs(y = "Number of genes", x = "",
       title = " ") 
