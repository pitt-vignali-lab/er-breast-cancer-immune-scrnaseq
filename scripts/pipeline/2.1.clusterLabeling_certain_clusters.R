library(Seurat)
library(tidyverse)

Dat <- readRDS("1.initialDat(res0.5,0.6).rds")
DefaultAssay(Dat) <- "RNA"

# markers for monocyte ----
png("../3.Result/2.clusterLabeling_Monocyte.png", width = 20, height = 30, units = "in", res = 300)
FeaturePlot(Dat, features = strsplit("CSF1R, MST1R, CD14, FCGR3A, MS4A7, TNFSF13", ", ")[[1]], label = T, split.by = "samples")
dev.off()

# markers for cluster 8 ----
p1 <- FeaturePlot(Dat, features = strsplit("SELL
CD69
TCF7
CCL5
CCR7
IL4R
CXCR6
KLF2
LEF1
", "\n")[[1]], label = T)

p2 <- VlnPlot(Dat, features = strsplit("SELL
CD69
TCF7
CCL5
CCR7
IL4R
CXCR6
KLF2
LEF1
", "\n")[[1]])


p1/p2

png("../3.Result/2.clusterLabeling_for_cluster_8.png", width = 10, height = 20, units = "in", res = 300)
p1/p2
dev.off()

# SingleR for all clusters ----
library(SingleR)
library(scater)

hpc <-  MonacoImmuneData()

anno.single <- SingleR(sumCountsAcrossCells(as.SingleCellExperiment(Dat), id = Dat@active.ident), 
                       ref = hpc, labels = hpc$label.fine, 
                       assay.type.test="sum")


annoDf <- cbind(annoDf, anno.single$labels)
annoDf <- data.frame(annoDf)
annoDf <- annoDf[-1]
# names(new.cluster.ids) <- levels(Dat)
# new.cluster.ids

annoDf <- cbind(levels(Dat), new.cluster.ids)
