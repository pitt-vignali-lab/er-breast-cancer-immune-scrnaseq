library(slingshot)
library(tidyverse)
library(Matrix)
library(Seurat)
library(scales)

dat <- readRDS("MonoMacro.RDS")

dat <- subset(dat, idents = 5, invert = T)

DimPlot(dat)


dat$cells <- rownames(dat@meta.data)
dat$clust <- dat@active.ident

datSub1 <- subset(dat, clust %in% c(1, 3))
sce <- as.SingleCellExperiment(datSub1, assay = "RNA")
rd <- reducedDims(sce)$UMAP
cl <- sce$clust
lin1 <- getLineages(rd, cl, start.clus = '1', end.clus = '3')
t <- getCurves(SlingshotDataSet(getCurves(lin1)))
lin1 <- as.data.frame(t@metadata[["curves"]][["Lineage1"]][["s"]])

datSub2 <- subset(dat, clust %in% c(1, 2, 4))
sce <- as.SingleCellExperiment(datSub2, assay = "RNA")
rd <- reducedDims(sce)$UMAP
cl <- sce$clust
lin2 <- getLineages(rd, cl, start.clus = '1', end.clus = '4')
t <- getCurves(SlingshotDataSet(getCurves(lin2)))
lin2 <- as.data.frame(t@metadata[["curves"]][["Lineage1"]][["s"]])
# -----

plot(lin1)
plot(rbind(lin1, lin2))
points(lin2, col = "red")
abline(v = 0)
abline(h = c(0.78,3))

lin1Sub1 <- lin1[lin1$UMAP_1 <0 & lin1$UMAP_2 < 4 & lin1$UMAP_2 > 0.78,]
lin1Sub2 <- lin1[!lin1$UMAP_1 %in% lin1Sub1$UMAP_1 & !lin1$UMAP_2 %in% lin1Sub1$UMAP_2, ]
lin2Sub1 <- lin2[lin2$UMAP_1 <0 & lin2$UMAP_2 < 4 & lin2$UMAP_2 > 0.78,]

lin1new <- rbind(lin2Sub1, lin1Sub2)

#----

umapDf <- data.frame(Embeddings(dat, 'umap'),
                     cluster = dat@active.ident,
                     samples = dat$samples) 

ggplot(umapDf, aes(x = UMAP_1, y = UMAP_2, col = cluster)) +
  geom_point() +
  geom_path(data = lin1new,  aes(x = UMAP_1, y = UMAP_2), 
            col = "#001f54", size = 1.5, arrow = arrow())+
  geom_path(data = lin2,  aes(x = UMAP_1, y = UMAP_2), 
            col = "#c44900", size = 1.5, arrow = arrow()) +
  theme_classic() +
  theme(text = element_text(size = 10))

ggsave("../3.Result/3.3.diffusionMapGeneSelection/slingshotOverall.png", width = 6, height = 5)

ggplot(umapDf, aes(x = UMAP_1, y = UMAP_2, col = ifelse(cluster %in% c(1, 3), cluster, "gray"))) +
  geom_point() +
  guides(color = guide_legend("cluster")) +
  scale_color_manual(breaks = c(1, 3, "gray"), 
                      values=c(hue_pal()(4)[c(1, 3)], "gray")) +
  geom_path(data = lin1new,  aes(x = UMAP_1, y = UMAP_2), 
            col = "#001f54", size = 1.5, arrow = arrow()) +
  geom_path(data = lin2,  aes(x = UMAP_1, y = UMAP_2), 
            col = "gray50", size = 1.5, arrow = arrow())+
  theme_classic() +
  theme(text = element_text(size = 10))
ggsave("../3.Result/3.3.diffusionMapGeneSelection/slingshotMonoToMono.png", width = 6, height = 5)

ggplot(umapDf, aes(x = UMAP_1, y = UMAP_2, col = ifelse(cluster %in% c(1, 2, 4), cluster, "gray"))) +
  geom_point() +
  guides(color = guide_legend("cluster")) +
  scale_color_manual(breaks = c(1, 2, 4, "gray"), 
                     values=c(hue_pal()(4)[c(1, 2, 4)], "gray")) +
  geom_path(data = lin1new,  aes(x = UMAP_1, y = UMAP_2), 
            col = "gray50", size = 1.5, arrow = arrow()) +
  geom_path(data = lin2,  aes(x = UMAP_1, y = UMAP_2), 
            col = "#c44900", size = 1.5, arrow = arrow())+
  theme_classic() +
  theme(text = element_text(size = 10))
ggsave("../3.Result/3.3.diffusionMapGeneSelection/slingshotMonoToMacro.png", width = 6, height = 5)

ggplot(umapDf, aes(x = UMAP_1, y = UMAP_2, col = cluster)) +
  geom_point() +
  geom_path(data = lin1new,  aes(x = UMAP_1, y = UMAP_2), 
            col = "#001f54", size = 1.5, arrow = arrow())+
  geom_path(data = lin2,  aes(x = UMAP_1, y = UMAP_2), 
            col = "#c44900", size = 1.5, arrow = arrow()) +
  facet_wrap(~samples, ncol = 2) +
  theme_classic() +
  theme(text = element_text(size = 10),
        strip.background = element_rect(colour="white", fill="white"))

ggsave("../3.Result/3.3.diffusionMapGeneSelection/slingshotSplit.png", width = 6, height = 5)

# ----

library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)

dat <- readRDS("MonoMacro.RDS")

dat <- subset(dat, idents = 5, invert = T)

dat <- subset(dat, idents = c(1, 3))

dataset <- wrap_expression(
  counts = t(dat@assays$RNA@counts),
  expression = t(dat@assays$RNA@data)
)

dataset <- add_prior_information(
  dataset,
  start_id = "m107_AACTCTTCATCCGCGA_1"
)

dataset <- add_grouping(
  dataset,
  Idents(dat)
)

model <- infer_trajectory(dataset, ti_slingshot())

plot_dimred(
  model, 
  expression_source = dataset$expression, 
  grouping = dataset$grouping
)

phm <- plot_heatmap(
  model,
  expression_source = dataset$expression,
  grouping = dataset$grouping,
  features_oi = 50
)

png("../3.Result/3.3.diffusionMapGeneSelection/slingshotHeatmap1-3.png", width = 5, height = 8, units = "in", res = 300)
phm
dev.off()