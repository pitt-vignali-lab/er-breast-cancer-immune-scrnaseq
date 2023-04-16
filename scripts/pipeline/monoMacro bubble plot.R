library(Seurat)
library(tidyverse)
require(scales)

DatLab <- readRDS("MonoMacro.RDS")
# id <- c(1, 1, 2, 3, 4, 5)
# names(id) <- levels(DatLab)
# DatLab <- RenameIdents(DatLab, id)
# saveRDS(DatLab, "MonoMacro.RDS")

DimPlot(DatLab, label = T) + NoLegend()

identities <- levels(DatLab@active.ident)
my_color_palette <- c(hue_pal()(length(identities) - 1), "#A65E2E")
names(my_color_palette) <- levels(DatLab@active.ident)
col_df <- data.frame(DatLab@active.ident)
col_df$ident <- rownames(col_df)
col_df_1 <- data.frame(my_color_palette)
col_df_1$DatLab.active.ident <- rownames(col_df_1)
col_df <- merge(col_df, col_df_1, by = "DatLab.active.ident")
col_df <- col_df[order(match(col_df$ident, names(DatLab@active.ident))),]

png("../3.Result/MonoMacro_final_buble.png", width = 6, height = 6, units = "in", res = 300)

DimPlot(DatLab, label = F) + 
  geom_point(aes(x = DatLab[["umap"]]@cell.embeddings[,1],
                 y = DatLab[["umap"]]@cell.embeddings[,2]),
             shape=21, color = "gray20", size = 2, fill = col_df$my_color_palette) +
  theme_void() + 
  theme(legend.position="none")

dev.off()
