library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

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

datSub <- subset(dat, idents = 2, samples == "IDC")

datSub <- subset(dat, samples == "IDC")

DimPlot(dat)

matExp <- AverageExpression(dat, group.by = c('samples', 'seurat_clusters'),
                            assays = "integrated")[[1]]

matExp <- AverageExpression(datSub, group.by = c('samples', 'reps'),
                            assays = "RNA")[[1]]

matExp <- AverageExpression(datSub, group.by = c('seurat_clusters'),
                            assays = "RNA")[[1]]

matExp[1:5, 1:10]

# heatmap ----
matExpSub <- matExp[tarG[c(1:45, 86:125)],] # IDC
matExpSub <- apply(matExpSub, 2, function(x) (x - min(x))*100 / (max(x) - min(x)))

matExpSub <- matExp[tarG[c(46:85, 126:145)],]
colnames(matExpSub) <- c("IDC C2", "ILC C2", "IDC C3", "ILC C3")

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

htmp <- BBmisc::normalize(as.matrix(matExpSub), range = c(-2, 2))
# htmp <- matExpSub
colnames(htmp) <- colnames(matExpSub)

Heatmap(htmp,
        na_col = "black",
        rect_gp = gpar(col = "gray50", lwd = 0.8),
        row_dend_width = unit(2, "cm"),
        # col = col_fun,
        # right_annotation = ha,
        row_names_gp = gpar(fontsize = 7.5),
        cluster_rows = T,
        show_row_dend = T,
        cluster_columns = F)


# circos ----
matExpSub <- matExp[tarG,  c("IDC_2", "ILC_2", "IDC_3", "ILC_3")]
colnames(matExpSub) <- c("IDC C2", "ILC C2", "IDC C3", "ILC C3")

col_fun = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

htmp <- BBmisc::normalize(as.matrix(matExpSub), range = c(-2, 2))
htmp <- matExpSub
colnames(htmp) <- colnames(matExpSub)

png("../3.Result/MonoMacro_selected_genes_circos_heatmap.png",
    width = 5, height = 5, units = "in", res = 300)
circos.clear()
circos.par(gap.after = c(2, 2, 2, 12))
circos.heatmap(htmp, split = split, col = col_fun,
               track.height = 0.15,
               dend.side = "inside",
               rownames.side = "outside",
               show.sector.labels = TRUE,
               cluster = T)

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 4) { # the last sector
    cn = colnames(htmp)
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(0.1, "mm"), #x coordinate
                (1:n) *0.2 - 0.9,#y coordinate
                cn, #label
                cex = 0.5, adj = c(0, 1), facing = "inside")
  }
}, bg.border = NA)

lgd = Legend(title = "Normalized \nGene Expression", col_fun = col_fun, 
             direction = "horizontal",
             title_position = "topcenter")
grid.draw(lgd)

dev.off()


# genes ----
tarG <- strsplit("IRF7
MYOF
CEBPD
PRAM1
GZMB
CSNK2A1
NR4A2
SIRPG
ITK
TNF
KLRB1
DIP2A
TNFSF10
MS4A6A
CCL5
PPP2R2B
C2orf49
SLC11A1
SMIM14
CLEC10A
BCAS2
AKT2
ACAA1
AREG
MTCH2
RBP7
S100A9
PSAP
PLSCR1
LILRB2
DOCK7
RETN
LYZ
HMOX1
CD36
TGFBI
SERPINA1
LST1
FCER1G
HLA-DRB5
MILR1
CSTA
AIF1
KLF10
CD68
PRAM1
CCL5
RAB24
RBM43
IL10RA
CLEC2D
WDR54
TSEN54
DUSP11
TIGIT
VPS13A
KLF2
SVIP
IL32
RHOH
BCL2L13
TYMP
ACAP1
GATA3
GZMH
IFITM3
AIF1
HMOX1
HLA-DPA1
CD302
LGALS3
HLA-DRB5
CD14
GPBAR1
SAT1
HLA-DRA
TIMP1
JUN
FCN1
S100A8
LYZ
MS4A7
SERPINA1
C5AR1
S100A9
CTSB
ITM2C
NEAT1
GZMK
RETN
ARFIP1
FTH1
PIK3IP1
SH2D1A
CXCR3
JUN
TNRC6C
CD74
RBMS1
LLPH
CD69
UQCRC2
SF1
GOLGA8B
CS
MNDA
CSF3R
LGALS3
S100A10
MPEG1
NINJ1
PHLDA2
LGALS2
CDA
RGS18
GZMB
RBP7
RNASE2
FCER1G
S100A9
MS4A6A
CSTA
CD36
FCN1
SYNJ1
GABARAP
CXCR4
MALAT1
HLA-E
TSC22D3
LTB
BTG1
RPS5
S100A6
ISG20
STXBP2
PLBD1
VIM
FCN1
MNDA
LGALS2
KLF10
LYZ
VCAN
S100A9", split = "\n")[[1]]

split = c(rep("IDC C2 DE genes", 45),
          rep("ILC C2 DE genes", 40),
          rep("IDC C3 DE genes", 40),
          rep("ILC C3 DE genes", 20))