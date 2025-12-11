
##PLOT 
install.packages("ggpubr")
library("ggpubr")

final_version <- rbind(dataset, dataset_new)
final_test_varsion_1 <- rbind(score_test_1, score_test_2)
final_test_version_2 <- rbind(final_version, final_test_varsion_1)
final_test <- rbind(final_test_version_2, score_test_3)
final_test$Method <- factor(final_test$Method, levels = c("Score_1", "Score_2", "Jaccard", "SOV99", "SOVrefine"), 
                            ordered = TRUE)
final_test_2 <- final_test
final_test_2$Score <- ifelse(final_test_2$Score > 1 &
                               final_test_2$Method %in% c("SOV99", "SOVrefine"),
                             1, final_test_2$Score)

score_test_1$Score <- ifelse(score_test_1$Score > 1 &
                               score_test_1$Method %in% c("SOV99"),
                             1,score_test_1$Score)

any(score_test_1$Score > 1)

#renombrar
score_test_2$Method <- "Length accuracy"

final_test_js <- rbind(final_test_true, score_test_1)
final_test_true <- rbind(dataset_new, score_test_2)

final_test_js <- final_test_js %>%
  mutate(
    Method   = str_squish(Method),
    Method   = str_remove(Method, "^'+"),                           # quita comilla inicial
    Method   = recode(Method,
                      "Presence/Absence" = "Presence/Absence"), # corrige typo
    Method   = fct_relevel(Method,
                           "Presence/Absence",
                           "Length accuracy",
                           "SOV99"),
    Software = fct_relevel(Software, "DeepTMHMM", "TOPCONS2")       # opcional: orden en X
  )

plot <- ggplot(data = final_test_js, aes(x = Software, y = Score)) +
  theme_bw() +
  geom_violin(scale = "width", alpha = 0.15, aes(fill = Software),
              show.legend = FALSE, trim = FALSE, linewidth = 0.4) +
  geom_point(aes(fill = Software), shape = 21, size = 2, alpha = 0.8,
             show.legend = FALSE, position = position_jitter(width = 0.05)) +
  geom_boxplot(width = 0.05, fill = "grey90", outlier.shape = NA) +
  geom_errorbar(stat = "boxplot", width = 0.5, linewidth = 1, color = "red",
                aes(ymin = after_stat(middle), ymax = after_stat(middle))) +
  scale_fill_manual(values = c("DeepTMHMM" = "#b9ffffff", "TOPCONS2" = "#fab8c2ff")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1.3)) +
  facet_wrap(~ Method, nrow= 1) +
  labs(y = expression("Score")) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, vjust = 2),
        strip.background = element_blank(),
        strip.text = element_text(size = 14)) +
  #geom_signif(comparisons = list(c("DeepTMHMM", "TOPCONS2")),
              #map_signif_level = TRUE,
              #textsize = 7,
              #y_position = 1.2) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("DeepTMHMM", "TOPCONS2")),
                     label = "p.signif", size = 5, paired = TRUE, label.y = 1.2)

plot 

?geom_signif

ggsave(filename = "scores_graph.jpg", plot = plot, width = 10, height = 7)
