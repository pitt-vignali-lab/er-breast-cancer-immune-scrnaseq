library(Seurat)
library(tidyverse)
require(scales)
library(ggpubr)

Dat <- readRDS("1.initialDat(res0.5,0.6).rds")
DefaultAssay(Dat) <- "RNA"

ClusterId <- c("CD4 T cells", # 0
               "Naive T cells", # 1
               "CD8 T cells", # 2
               "CD8 T cells", # 3
               "NK cells", # 4
               "monocytes \nand macs", # 5
               "B cells", # 6
               "monocytes \nand macs", # 7
               "tissue resident \n T cells", # 8
               "Tregs", # 9
               "monocytes \nand macs", # 10
               "DCs", # 11
               "pDCs", # 12
               "monocytes \nand macs" # 13
)

names(ClusterId) <- levels(Dat)

DatLab <- RenameIdents(Dat, ClusterId)

features <- strsplit("CD3E
CD3D
CD8A
CD4
IL7R
CCR7
SELL
TCF7
GZMK
EOMES
IFNG
FOXP3
CD69
CCL5
CTLA4
GNLY
NKG7
NCAM1
MS4A1
CD79A
CD1C
FCER1A
HLA-DRA
IL3RA
CLEC4C
CD14
MS4A7
CSF1R
CD68
LYZ
", "\n")[[1]]

p1 <- DotPlot(DatLab, features = features[1:18], cols = c("blue", "red"))
p2 <- DotPlot(DatLab, features = features[19:30], cols = c("blue", "red"))
p1$data$type = "Lymphoid "
p2$data$type = "Myeloid"
p <- rbind(p1$data, p2$data)
p$pct.exp <- log1p(p$pct.exp)

png("../3.Result/2.clusterLabeling_marker_genes_splited.png", width = 8, height = 6, units = "in", res = 300)
ggplot(p) +
  geom_point(aes(x = features.plot, y = id, size = pct.exp, col = avg.exp.scaled))+
  ylab("Cell type") +
  xlab("Marker gene") +
  labs(color = "Scaled average expression",
       size = "Average expression level (%)") +
  scale_size(range = c(0.5, 6), breaks = c(1, 50, 99), labels = c(1, 50, 100))+
  scale_color_gradient2(low="blue4", mid="#F1FAEE", high="#E63946", 
                        midpoint=0, limits=range(-2.5, 2.5))+# range(p$data$avg.exp.scaled)) +
  theme(panel.background = element_blank(),
        panel.border=element_rect(fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background=element_blank(),
        axis.text.x=element_text(colour="black", angle = 45, hjust = 1),
        axis.text.y=element_text(colour="black"),
        axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line"),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key = element_blank())+
  guides(size = guide_legend(title.position="top", title.hjust = 0.5, order = 1),
         color = guide_colourbar(title.position="top", title.hjust = 0.5, order = 2)) +
  facet_grid(~type, scales = "free_x", space = "free_x")
dev.off()

## table output
pSub = p[c('id', 'features.plot',  'pct.exp', 'avg.exp.scaled')]
names(pSub) <- c("Cell type", "Genes", "Average expression level",
                 "Scaled average expression")
write.csv(pSub, "../3.Result/SourceData/fig4b_markerGene.csv")