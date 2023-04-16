library(tidyverse)
library(Seurat)
library(ggbeeswarm)

dat <- readRDS("../TIL_mono_macro.RDS")

t <- DimPlot(dat)


df <- read.csv("newMeta.csv")

dfNew <- df %>% group_by(cell_clusters) %>% summarise(mean(latent_time))

ggplot(df, aes(x = latent_time, 
               y = as.character(cell_clusters),
               col = as.character(cell_clusters))) +
  geom_boxplot(size = 1) +
  geom_quasirandom(groupOnX = FALSE, alpha = 0.2)+
  scale_color_brewer(palette="Set2") +
  ylab("Clusters") +
  xlab("latent-time") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size = 15, face = 2))
ggsave("../../3.Result/3.monomac_velocity_latent_time_with_box.png",
       width = 5, height = 3, units = "in", dpi = 300, bg = "white")


dfSub <- df[c('CellID', 'latent_time', 'cell_clusters')] |> 
  pivot_wider(names_from = cell_clusters, values_from = latent_time)

write.csv(dfSub, "../../3.Result/SourceData/fig5c_latent-time.csv")

testDf <- data.frame()
for (i in 2:7) {
  testDf <- rbind(data.frame(value = na.omit(dfSub[,i]),
                             group = paste0("clst", names(dfSub[,i]))) %>% 
                    rename('value' = 1),
                  testDf)
}

table(testDf$group)
result = kruskal.test(value~ group,
                      data = testDf)

print(result)

library(ggpubr)

result <- data.frame()
for (i in unique(testDf$group)) {
  t = compare_means(value ~ group,  data = testDf,
                    ref.group = i, method = "wilcox.test")
  
  result <- rbind(t, result)
}

write.csv(result, "wilcoxTestAccrossAllClusters.csv")
