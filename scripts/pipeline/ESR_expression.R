library(Seurat)
library(tidyverse)
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


FeaturePlot(DatLab, features = c("ESR1", "ESR2"), label = T, split.by = "samples")

p1 = VlnPlot(DatLab, features = "ESR1", 
        split.by = "samples", 
        ncol = 1,
        pt.size = 0) + 
  geom_jitter(mapping = aes(color = split), 
                # data = p$data, 
                position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9))+
  NoLegend()

p2 = VlnPlot(DatLab, features = "ESR2", 
             split.by = "samples", 
             ncol = 1,
             pt.size = 0) + 
  geom_jitter(mapping = aes(color = split), 
              # data = p$data, 
              position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.9))+
  guides(fill = FALSE) +
  labs(col = "Samples") +
  theme(legend.position = 'right')
p1 / p2
