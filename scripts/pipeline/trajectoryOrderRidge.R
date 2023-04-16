library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)
library(ggbeeswarm)
library(ggridges)

dat <- readRDS("MonoMacro.RDS")

DimPlot(dat)

dat <- subset(dat, idents = 5, invert = T)

pca_df <- data.frame(PC1 = Embeddings(dat)[,1],
                     PC2 = Embeddings(dat)[,2],
                     cell_type2 = dat@active.ident)

ggplot(pca_df, aes(x = PC1, group = cell_type2, 
                   fill = cell_type2)) +
  geom_density(position="fill") +
  xlab("Pseudo-time ->")

ggplot(pca_df, aes(x = PC1, group = cell_type2, 
                   fill = cell_type2,
                   col = cell_type2)) +
  geom_density(alpha = 0.5) +
  xlab("Pseudo-time →") +
  ylab("Cell density") +
  guides(color = guide_legend("Cluster"),
         fill = guide_legend("Cluster")) +
  theme_bw() +
  theme(text = element_text(size = 15, face = 2))
ggsave("../3.Result/pseudoTimeorder.png", width = 10, height = 5, bg = "White")
