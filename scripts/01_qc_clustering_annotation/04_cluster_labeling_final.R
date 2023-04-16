library(Seurat)
library(tidyverse)
require(scales)

Dat <- readRDS("1.initialDat(res0.5,0.6).rds")
DefaultAssay(Dat) <- "RNA"

ClusterId <- c("CD4 T cells", # 0
               "Naive T cells", # 1
               "CD8 T cells", # 2
               "CD8 T cells", # 3
               "NK cells", # 4
               "monocytes and macs", # 5
               "B cells", # 6
               "monocytes and macs", # 7
               "Effector T cells", # 8
               "Tregs", # 9
               "monocytes and macs", # 10
               "DCs", # 11
               "pDCs", # 12
               "monocytes and macs" # 13
)

names(ClusterId) <- levels(Dat)

DatLab <- RenameIdents(Dat, ClusterId)

# datTreg <- subset(DatLab, idents = 'Tregs')
# saveRDS(datTreg, '../1.Data/r35TregSubset.RDS')

png("../3.Result/2.clusterLabeling_final.png", width = 8, height = 8, units = "in", res = 300)
DimPlot(DatLab, label = T) + NoLegend()
dev.off()

identities <- levels(DatLab@active.ident)
my_color_palette <- hue_pal()(length(identities))
names(my_color_palette) <- levels(DatLab@active.ident)
col_df <- data.frame(DatLab@active.ident)
col_df$ident <- rownames(col_df)
col_df_1 <- data.frame(my_color_palette)
col_df_1$DatLab.active.ident <- rownames(col_df_1)
col_df <- merge(col_df, col_df_1, by = "DatLab.active.ident")
col_df <- col_df[order(match(col_df$ident, names(DatLab@active.ident))),]

png("../3.Result/2.clusterLabeling_final_buble.png", width = 8, height = 8, units = "in", res = 300)
DimPlot(DatLab, label = F) + 
  geom_point(aes(x = DatLab[["umap"]]@cell.embeddings[,1],
                 y = DatLab[["umap"]]@cell.embeddings[,2]),
             shape=21, color = "gray20", size = 2, fill = col_df$my_color_palette) +
  theme_void() + 
  theme(legend.position="none")
dev.off()

png("../3.Result/2.clusterLabeling_final_buble_SourrceWise.png", width = 16, height = 16, units = "in", res = 300)
DimPlot(DatLab, label = F, split.by = "samples", cols = my_color_palette, ncol = 2)+ 
  geom_point(aes(x = DatLab[["umap"]]@cell.embeddings[,1],
                 y = DatLab[["umap"]]@cell.embeddings[,2]),
             shape=21, color = 'black', size = 2, fill = col_df$my_color_palette) + 
  theme_void() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.position="none")

dev.off()

table(DatLab@active.ident, DatLab$reps, DatLab$samples)

table(DatLab$integrated_snn_res.0.5, DatLab$reps, DatLab$samples)
