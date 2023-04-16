library(Seurat)
library(tidyverse)

dat <- readRDS("MonoMacro.RDS")

dat <- subset(dat, samples %in% c("IDC", "ILC"))

dat <- FindNeighbors(dat)

dat <- FindClusters(dat, resolution = 0.4)

pct <- dat[["pca"]]@stdev / sum(dat[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

dat <- RunUMAP(dat, dims = 1:pcs)

DimPlot(dat)

dfPlot <- data.frame(Embeddings(dat, reduction = "umap"),
                     Idents(dat),
                     dat$samples)


# Downsize ----
table(dat$samples)

ratio = table(dat$samples)["ILC"] / table(dat$samples)["IDC"]

dfPlotIDC = data.frame()

for (i in unique(Idents(dat))) {
  dfSub = dfPlot %>% filter(dat.samples == "IDC", Idents.dat. == i)
  t <- dfSub[sample(1:nrow(dfSub), nrow(dfSub) * ratio),]
  dfPlotIDC = rbind(t, dfPlotIDC)
}

nrow(dfPlotIDC)

dfPlotNew <- rbind(dfPlotIDC, dfPlot %>% filter(dat.samples == "ILC"))

ggplot(dfPlotNew, aes(UMAP_1, UMAP_2, fill = Idents.dat.)) +
  geom_point(shape=21, color = "gray20", size = 3) +
  scale_fill_brewer(palette="Set2") +
  # facet_grid(~dat.samples) +
  coord_fixed(ratio = 1) +
  theme_void()

ggsave("../3.Result/mono_macor_re-cluster_down_zise.png", width = 5, height = 5, 
       units = "in", dpi = 300, bg = "white")

# Density plot ----
ggplot(dfPlot, aes(UMAP_1, UMAP_2)) +
  stat_density_2d(aes(fill = ..level..), 
                  geom = "polygon", bins = 500)+
  scale_fill_gradient(low = "white", high = "black") +
  facet_grid(~dat.samples) +
  coord_fixed(ratio = 1) +
  xlim(-8, 8) +
  ylim(-8, 8) +
  theme_void()

ggsave("../3.Result/mono_macor_re-cluster_density.png", width = 8, height = 3, 
       units = "in", dpi = 300, bg = "white")

# Hex plot ----
ggplot(dfPlot, aes(UMAP_1, UMAP_2, colour = ..count..)) +
  geom_hex(aes(colour = ..count..), bins = 30) +
  scale_fill_continuous(type = "viridis") +
  scale_color_continuous(type = "viridis") +
  guides(colour = 'none') +
  # facet_grid(~dat.samples) +
  coord_fixed(ratio = 1) +
  xlim(-8, 8) +
  ylim(-8, 8) +
  theme_void() 

ggsave("../3.Result/mono_macor_re-cluster_hex.png", width = 5, height = 5, 
       units = "in", dpi = 300, bg = "white")
