library(Seurat)
library(tidyverse)

dat <- readRDS("MonoMacro.RDS")

dat <- subset(dat, samples %in% c("IDC", "ILC")) #IDC, ILC, PBL.IDC, PBL.ILC

dat <- FindNeighbors(dat)

dat <- FindClusters(dat, resolution = 0.4)

DimPlot(dat)

VlnPlot(dat, 
        assay = "RNA",
        features = c('CD14', 'FCGR3A', 'CSF1R', 'MARCO', 'CX3CR1', 'CDKN1C', 'APOE', 'SPP1', 'MKI67', 'FCGR1A', 'CD68'),
        ncol = 6)


DotPlot(dat, 
        assay = "RNA",
        features = c('CD14', 'FCGR3A', 'CSF1R', 'MARCO', 'CX3CR1', 'CDKN1C', 'APOE', 'SPP1', 'MKI67', 'FCGR1A', 'CD68'),
        col.min = 0,
        cols = c("white", "blue")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  

# OLD ----
Dat <- readRDS("1.initialDat(res0.5,0.6).rds")
DefaultAssay(Dat) <- "RNA"

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

png("../3.Result/3.monomac_overall.png", width = 6, height = 6, units = "in", res = 300)
DimPlot(DatSub, reduction = "umap", label = T, pt.size = 1)+NoLegend()
dev.off()

DefaultAssay(DatSub) <- "RNA"

table(Idents(DatSub), DatSub$reps, DatSub$samples)

png("../3.Result/3.monomac_split.png", width = 20, height = 5, units = "in", res = 300)
DimPlot(DatSub, reduction = "umap", label = T, pt.size = 1, split.by = "samples")+NoLegend()
dev.off()

png("../3.Result/3.monomac_reps_split.png", width = 20, height = 5, units = "in", res = 300)
DimPlot(DatSub, reduction = "umap", pt.size = 1, group.by = "reps", split.by = "samples")
dev.off()

# Overall cell distribution (sum)----
png("../3.Result/3.monomac_cell_number_sum_distribution.png", width = 10, height = 10, units = "in", res = 300)
data.frame(prop.table(table(Idents(DatSub), DatSub$samples),
                      margin = 2)) %>% 
  # filter(Freq != 0) %>%
  rename(CellType = Var1, Sample = Var2) %>%
  ggplot(aes(x = Sample, y = Freq*100, fill=CellType, label = paste0(round(Freq, 2) * 100,"%"))) +
  geom_bar(stat='identity', colour="black", alpha = 0.5) +
  ylab("Persentage (%)") +
  theme_bw()
dev.off()

# every cell distribution (mean)-----

data <- data.frame(prop.table(table(Idents(DatSub), 
                                    DatSub$samples, DatSub$reps),
                              margin = 2)) %>% 
  filter(Freq != 0)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(6)

# Calculates mean, sd, se and IC
for (i in 0:5) {
  assign(paste0("p", i),
         data %>%
           filter(Var1 == i) %>%
           group_by(Var2) %>%
           summarise( 
             n=n(),
             sum=sum(Freq),
             mean=mean(Freq),
             sd=sd(Freq)
           ) %>%
           mutate( se=sd/sqrt(n))  %>%
           mutate( ic=se * qt((1-0.05)/2 + .5, n-1)) %>%
           ggplot() +
           ylim(0, 0.125)+
           xlab("Sample type") +
           geom_bar( aes(x=Var2, y=mean), stat="identity", fill=cols[i+1], alpha=0.5) +
           geom_errorbar( aes(x=Var2, ymin=mean - se, ymax=mean + se), width=0.2, colour="black", alpha=0.9, size=0.5) +
           ggtitle(paste0("Cluster ", i)) +
           theme(plot.title = element_text(hjust = 0.5))+
           theme_bw())
}

png("../3.Result/3.monomac_cell_number_mean_distribution.png", width = 10, height = 10, units = "in", res = 300)
p0+p1+p2+p3+p4+p5
dev.off()

wilcoxRes <- c()
for (i in 0:5) {
  testDf <- data %>%
    filter(Var1 == i) %>%
    group_by(Var2)
  
  testRes <- rbind(
  wilcox.test(testDf[testDf$Var2 == "IDC",]$Freq,
              testDf[testDf$Var2 == "PBL.IDC",]$Freq)$p.value,
  wilcox.test(testDf[testDf$Var2 == "ILC",]$Freq,
              testDf[testDf$Var2 == "PBL.ILC",]$Freq)$p.value,
  wilcox.test(
    c(testDf[testDf$Var2 == "IDC",]$Freq,
      testDf[testDf$Var2 == "ILC",]$Freq),
    c(testDf[testDf$Var2 == "PBL.IDC",]$Freq,
     testDf[testDf$Var2 == "PBL.ILC",]$Freq)
  )$p.value)
  wilcoxRes <- cbind(wilcoxRes, testRes)
  
}

wilcoxRes <- data.frame(wilcoxRes)

rownames(wilcoxRes) <- c("IDCvsPBL.IDC", "ILCvsPBL.ILC", "TILvsPBL")

colnames(wilcoxRes) <- paste0("Clst", 0:5)

# Markers ----
# Classical monocytes

p1 <- FeaturePlot(DatSub, features = c("CD14", "FCGR3A", "FCN1"), split.by = "samples", label = T)

p2 <- VlnPlot(DatSub, features = c("CD14", "FCGR3A", "FCN1"), 
        group.by = "samples", 
        split.by = "integrated_snn_res.0.2", 
        ncol = 2) + theme(legend.position = 'right')

png("../3.Result/3.monomac_markers_for_Classical_monocytes.png", width = 12, height = 15, units = "in", res = 300)
p1/p2
dev.off()

# Intermediate monocytes

p1 <- FeaturePlot(DatSub, features = c("CD14", "FCGR3A", "ITGAL", "CDKN1C"), split.by = "samples", label = T)

p2 <- VlnPlot(DatSub, features = c("CD14", "FCGR3A", "ITGAL", "CDKN1C"), 
              group.by = "samples", 
              split.by = "integrated_snn_res.0.2", 
              ncol = 2) + theme(legend.position = 'right')

png("../3.Result/3.monomac_markers_for_Intermediate_monocytes.png", width = 12, height = 18, units = "in", res = 300)
p1/p2
dev.off()

# Non-classical monocytes

p1 <- FeaturePlot(DatSub, features = c("CD14", "FCGR3A"), split.by = "samples", label = T)

p2 <- VlnPlot(DatSub, features = c("CD14", "FCGR3A"), 
              group.by = "samples", 
              split.by = "integrated_snn_res.0.2", 
              ncol = 2) + theme(legend.position = 'right')

png("../3.Result/3.monomac_markers_for_Non-classical_monocytes.png", width = 12, height = 10, units = "in", res = 300)
p1/p2
dev.off()

# Up in macrophages

p1 <- FeaturePlot(DatSub, features = c("MRC1", "APOE", "C1QA", "C1QB", "CTSB", "CTSD"), split.by = "samples", label = T)

p2 <- VlnPlot(DatSub, features = c("MRC1", "APOE", "C1QA", "C1QB", "CTSB", "CTSD"), 
              group.by = "samples", 
              split.by = "integrated_snn_res.0.2", 
              ncol = 2) + theme(legend.position = 'right')

png("../3.Result/3.monomac_markers_for_Up_in_macrophages.png", width = 12, height = 20, units = "in", res = 300)
p1/p2
dev.off()

p1 <- FeaturePlot(DatSub, features = c("OLR1", "CDKN1C", "S100A12", "CXCL5", "SPP1", "APOE",
                                       "MMP7", "CHIT1", "CCL18", "FABP4", "CXCL9", "CXCL12",
                                       "CCL18", "MKI67"), split.by = "samples", label = T)

p2 <- VlnPlot(DatSub, features = c("OLR1", "CDKN1C", "S100A12", "CXCL5", "SPP1", "APOE",
                                  "MMP7", "CHIT1", "CCL18", "FABP4", "CXCL9", "CXCL12",
                                  "CCL18", "MK167"), 
              group.by = "samples", 
              split.by = "integrated_snn_res.0.2", 
              ncol = 2) + theme(legend.position = 'right')

png("../3.Result/3.monomac_markers_from_Zilionis.png", width = 12, height = 70, units = "in", res = 300)
p1/p2
dev.off()
