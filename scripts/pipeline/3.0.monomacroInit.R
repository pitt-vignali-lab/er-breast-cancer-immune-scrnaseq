library(Seurat)
library(tidyverse)

Dat <- readRDS("1.initialDat(res0.5,0.6).rds")

DatSub <- subset(Dat, idents = c(5, 7, 10, 13))

DimPlot(DatSub)

#----

DefaultAssay(DatSub) <- "integrated"

DatSub <- ScaleData(DatSub)

DatSub <- RunPCA(object = DatSub)

ElbowPlot(DatSub, ndims = 50)

#----
# Determine percent of variation associated with each PC
pct <- DatSub[["pca"]]@stdev / sum(DatSub[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
# ----

DatSub <- RunUMAP(DatSub,
                  dims = 1:pcs,
                  reduction = "pca")

DatSub <- FindNeighbors(DatSub, dims = 1:pcs)

DefaultAssay(DatSub) <- "integrated"

DatSub <- FindClusters(DatSub, resolution = 0.2)

DefaultAssay(DatSub) <- "RNA"

DimPlot(DatSub)

rm(Dat, plot_df, co1, co2, cumu, pct)

# Add meta data ----
metaDf <- DatSub@meta.data["samples"]
metaDf1 <- data.frame(samples = c("IDC", "ILC", "PBL.IDC", "PBL.ILC"),
                      types = rep(c("TIL", "PBL"), each  = 2))
metaDf <- metaDf %>%
  rownames_to_column("ID") %>%
  merge(metaDf1, by = "samples") %>%
  column_to_rownames(var = "ID") %>%
  select(-"samples")

DatSub <- AddMetaData(DatSub, metaDf)

rm(metaDf, metaDf1)

metaDf <- DatSub@meta.data["samples"]
metaDf1 <- data.frame(samples = c("IDC", "ILC", "PBL.IDC", "PBL.ILC"),
                      types = rep(c("TIL", "PBL"), each  = 2))
metaDf <- metaDf %>%
  rownames_to_column("ID") %>%
  merge(metaDf1, by = "samples") %>%
  column_to_rownames(var = "ID") %>%
  select(-"samples")

DatSub <- AddMetaData(DatSub, metaDf)

rm(metaDf, metaDf1)

saveRDS(DatSub, "MonoMacro.RDS")

# # DE ----
# DatSub <- readRDS("MonoMacro.RDS")
# 
# DatSub$des <- paste0(DatSub$samples, DatSub$integrated_snn_res.0.2)
# 
# Idents(DatSub) <- "des"
# 
# DEgene <- FindMarkers(DatSub, ident.1 = "ILC2", ident.2 = "IDC2")
# 
# write.csv(DEgene, "../3.Result/DE_TIL_ILC2vsIDC2.csv")
# 
# 
# 
# for (i in 0:5) {
#   subSub <- subset(DatSub, idents = i)
#   Idents(subSub) <- "types"
#   de <- FindMarkers(subSub, ident.1 = "TIL", ident.2 = "PBL")
#   write.csv(de, paste0("../3.Result/3.monomac_DE_Clst_", i, "_TILvsPBL.csv"))
#   rm(i, de, subSub)
# }
# 
# pbmc <- FindVariableFeatures(subSub, selection.method = "disp", nfeatures = 2000)
# top10 <- head(VariableFeatures(pbmc), 10)
# plot1 <- VariableFeaturePlot(pbmc)
# plot2 <- LabelPoints(plot = plot1, points = top10_1$X, repel = TRUE)
# plot1 + plot2
# top10_1 <- read.csv("../3.Result/3.monomac_DE_Clst_0_TILvsPBL.csv")[1:10,]
# 
# top10_1$X[order(top10_1$X)]
# 
# top10[order(top10)]
# 
# # heat map validation ----
# library(ComplexHeatmap)
# library(circlize)
# library(BBmisc)
# 
# htfun <- function(x){
#   df.1 <- subset(DatSub, types == "TIL" &  integrated_snn_res.0.2 ==  x)
#   df.1 <- AverageExpression(df.1, features = gList, assays = "RNA", group.by = "reps")[["RNA"]]
#   
#   df.2 <- subset(DatSub, types == "PBL" &  integrated_snn_res.0.2 ==  x)
#   df.2 <- AverageExpression(df.2, features = gList, assays = "RNA", group.by = "reps")[["RNA"]]
#   
#   return(cbind(df.1,df.2))
#   
#   cat(ncol(df.1))
# }
# 
# gList <- read.csv("../3.Result/3.monomac_DE_Clst_0_TILvsPBL.csv")$X
# 
# gclt1 <- htfun(0)
# 
# htmp <- normalize(as.matrix(gclt1),
#                   range = c(-4, 4))
# colnames(htmp) <- colnames(gclt1)
# 
# col_fun = circlize::colorRamp2(c(-2, 0, 2), c("#4b5cc4", "white", "#ff4e20"))
# 
# # png("DEHeatmapKM.png", width = 6, height = 20, units = "in", res = 300)
# # pdf("DEHeatmap.pdf", width = 6, height = 40)
# Heatmap(htmp,
#         na_col = "black",
#         col = col_fun,
#         # row_km = 3, 
#         #　row_order = 1:length(coreDEg)*3,
#         row_names_gp = gpar(fontsize = 2),
#         # cluster_rows = F,
#         # show_row_dend = F,
#         # cluster_columns = F,
#         # show_column_dend = F,
#         # row_split = rep(c("Memory CD8", "Dysfunctional CD8"), c(nrow(gclt1), nrow(gclt2))),
#         # column_split = rep(1:2, c(13,15)),
#         column_names_side = "bottom",
#         column_names_rot = 90,
#         column_names_gp = gpar(fontsize = 8),
#         # column_title = NULL,
#         show_column_names = T,
#         name = "Normalized Gene Expression",
#         top_annotation = HeatmapAnnotation(Type = rep(c("TIL", "PBL"), c(13,15)),
#                                            show_annotation_name = F,
#                                            col = list(Type = c("TIL" = "#ffa631", "PBL" = "#70f3ff")),
#                                            annotation_legend_param = list(Type = list(title = "Cancer Type",
#                                                                                       title_position = "leftcenter-rot"))),
#         heatmap_legend_param = list(# direction = "horizontal",
#           title_position = "leftcenter-rot",
#           legend_width = unit(5, "in")))
# 
# # dev.off()    
# 
# library(pheatmap)
# 
# htmp1 <- htmp
# annotation_col = data.frame(type = rep(c("TIL", "PBL"), c(13,15)))
# colnames(htmp1) <- 1:28
# pheatmap(htmp1,show_rownames=FALSE,show_colnames=FALSE,
#          annotation_col=annotation_col,
#          scale = "none",clustering_method="centroid",
#          clustering_distance_cols="euclidean")
