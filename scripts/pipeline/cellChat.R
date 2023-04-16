library(Seurat)
library(tidyverse)
require(scales)

Dat <- readRDS("1.initialDat(res0.5,0.6).rds")

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

DimPlot(DatLab, label = T) + NoLegend()