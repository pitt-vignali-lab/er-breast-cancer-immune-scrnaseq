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

DimPlot(dat, split.by = "samples")

# saveRDS(dat, "TIL_mono_macro.RDS")

mks <- FindAllMarkers(dat)

mksSub <- mks %>% 
  filter(p_val_adj < 0.05) %>% 
  write.csv("../3.Result/mono_macor_DE_genes_in_all_cluster.csv", row.names = F)


# figure for manual ----

# UMAP

plotDf <- merge(Embeddings(dat, reduction = "umap"),
                dat@meta.data[c("samples", "integrated_snn_res.0.4")],
                by = 0)

names(plotDf)[5] <- "cluster"

ggplot(plotDf, aes(UMAP_1, UMAP_2, fill = cluster)) +
  geom_point(shape=21, color = "gray20", size = 3) +
  # scale_x_continuous(breaks=seq(1930,1940,2)) +
  scale_fill_brewer(palette="Set2") +
  theme_void() +
  coord_fixed() +
  theme()

ggsave("../3.Result/mono_macor_re-cluster.png", width = 6, height = 5, 
       units = "in", dpi = 300, bg = "white")

ggplot(plotDf, aes(UMAP_1, UMAP_2, fill = cluster)) +
  geom_point(shape=21, color = "gray20", size = 3) +
  # scale_x_continuous(breaks=seq(1930,1940,2)) +
  scale_fill_brewer(palette="Set2") +
  facet_grid(~samples) +
  theme_void() +
  coord_fixed() +
  theme()

ggsave("../3.Result/mono_macor_re-cluster_split.png", width = 10, height = 3, 
       units = "in", dpi = 300, bg = "white")

# marker genes 

mkgs <- strsplit("S100A4, S100A8, FCN1, CXCL8, IL1B, CD14, FCGR3A, CSF1R, MARCO, CX3CR1, CDKN1C, APOE, SPP1, MKI67, FCGR1A, CD68, C1QA, APOC1, C1QB, C1QC, FOLR2",
                  split = ", ")[[1]]

# marker genes ----
mkIdc <- DotPlot(subset(dat, samples == "IDC"), features = mkgs,
        assay = "RNA")

mkIlc <- DotPlot(subset(dat, samples == "ILC"), features = mkgs,
                 assay = "RNA")

rbind(mkIdc$data |> add_column(sample = "IDC"),
               mkIlc$data |> add_column(sample = "ILC")) |>
  mutate(id = paste0("Cluster ", id)) |>
  filter(pct.exp > 10) |>
  ggplot(aes(x = features.plot, y = id, size = pct.exp, col = avg.exp.scaled)) +
  geom_point(alpha = 0.8)+
  scale_colour_gradient2(name = "Average \nExpression", 
                         low="blue", mid="white", high="red",
                         midpoint = 0,
                         limits = c(-2.5, 2.5))+ 
  facet_wrap(~sample, ncol = 1) +
  guides(size = guide_legend(title="Persent \nExpression"),
         color = guide_colourbar(barwidth = 0.8)) +
  ylab("") +
  xlab("Genes") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust =  1),
        text = element_text(size = 12, face = 2),
        strip.background = element_blank())

ggsave("../3.Result/new_monomarc_markergenes_dotplot.png", width = 7, height = 5,
       units = "in", dpi = 300, bg = "white")

# test ----
DotPlot(subset(dat, samples == "ILC"), features = mkgs, 
        assay = "RNA", cols = c("white", "blue"), col.min = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust =  1))


VlnPlot(dat, features = mkgs, 
        assay = "RNA") &
  scale_fill_brewer(palette="Set2")

VlnPlot(subset(dat, samples == "IDC"), features = mkgs, 
        assay = "RNA") &
  scale_fill_brewer(palette="Set2")

ggsave("../3.Result/mono_macor_re-cluster_selected_genes_in_all.png", width = 7, height = 10, 
       units = "in", dpi = 300, bg = "white")

DotPlot(dat, features = mkgs,
        assay = "RNA", cols = c("white", "blue"), col.min = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust =  1))