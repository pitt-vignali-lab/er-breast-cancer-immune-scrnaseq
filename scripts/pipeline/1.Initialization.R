library(Seurat)
library(tidyverse)

# Demutiplex function ----

readat <- function(raw.data, sample.type){
  
  reps = as.character(substitute(raw.data))
  
  colnames(raw.data) <- gsub("-1", "", colnames(raw.data))
  
  raw.data <- CreateSeuratObject(raw.data)
  
  meta.samp <- data.frame(samples=c(rep(sample.type,ncol(raw.data))),
                          reps = unlist(strsplit(reps,"\\."))[1])#V
  
  rownames(meta.samp) <- colnames(raw.data)
  
  dat <- AddMetaData(raw.data,meta.samp)
  
  return(dat)
}

unmt <- function(raw.data, raw.cite, cite.1.hto, cite.1.anno, cite.2.hto, cite.2.anno){
  
  reps = as.character(substitute(raw.data))
  
  raw.cite <- raw.cite[grep("^totalseq*", rownames(raw.cite)),]
  
  raw.cite <- raw.cite[c(cite.1.hto, cite.2.hto),]
  
  colnames(raw.data) <- gsub("-1", "", colnames(raw.data))
  
  cells.use <- intersect(colnames(raw.data), colnames(raw.cite))
  
  dat <- CreateSeuratObject(counts=raw.data[,cells.use])
  
  dat[["ADT"]] <- CreateAssayObject(counts=raw.cite[,cells.use])
  
  DefaultAssay(dat) <- "ADT"
  
  dat <- NormalizeData(dat,normalization.method="CLR")
  
  dat.norm <- GetAssayData(dat)
  
  dat.norm <- data.frame(t(as.matrix(dat.norm)))
  
  k.means.mean <- c()
  
  for (i in 1:2) {
    
    k.means.mean[i] <- mean(kmeans(dat.norm[,i],centers=2, nstart = 25)$centers)
    
  }
  
  dat.norm.tib <- as_tibble(dat.norm,rownames="cell_barcode")
  
  samp1 <- dat.norm.tib[dat.norm.tib[2] > k.means.mean[1] & #V
                          dat.norm.tib[3] < k.means.mean[2],]$cell_barcode
  
  samp2 <- dat.norm.tib[dat.norm.tib[2] < k.means.mean[1] & 
                          dat.norm.tib[3] > k.means.mean[2],]$cell_barcode
  
  meta.samp <- data.frame(samples=c(rep(cite.1.anno,length(samp1)), #V
                                    rep(cite.2.anno,length(samp2))),
                          reps = unlist(strsplit(reps, "\\."))[1])#V
  
  rownames(meta.samp) <- c(samp1,samp2)
  
  dat <- dat[,rownames(meta.samp)]
  
  dat <- AddMetaData(dat,meta.samp)
  
  return(dat)
}

clulab <- function(data, exist, subset = c()){
  if (length(subset) != 0){
    list.e <- list()
    for (i in 1:length(exist)) {
      list.e[[i]] <- data[data$gene == exist[i],]$cluster
    }
    cl.e <- Reduce(intersect, list.e)
    
    list.s <- list()
    for (i in 1:length(subset)) {
      list.s[[i]] <- data[data$gene == subset[i],]$cluster
    }
    cl.s <- Reduce(intersect, list.s)
    
    return(cl.e[! cl.e %in% cl.s])
  } else {
    list.e <- list()
    for (i in 1:length(exist)) {
      list.e[[i]] <- data[data$gene == exist[i],]$cluster
    }
    cl.e <- Reduce(intersect, list.e)
    return(cl.e)
  }
  
}

# Data read in ----
# M899 ---

m899.pbl <- Read10X("../1.Data/gene_expression/M899_PBMC/filtered_gene_bc_matrices/hg19/")
m899.pbl <- readat(m899.pbl, "PBL.ILC")

m899.til <- Read10X("../1.Data/gene_expression/M899_TIL/filtered_gene_bc_matrices/hg19/")
m899.til <- readat(m899.til, "ILC")

# M615 ---

m615.pbl <- Read10X("../1.Data/gene_expression/M615_PBL/filtered_gene_bc_matrices/hg19/")
m615.pbl <- readat(m615.pbl, "PBL.IDC")

m615.til <- Read10X("../1.Data/gene_expression/M615_TIL/filtered_gene_bc_matrices/hg19/")
m615.til <- readat(m615.til, "IDC")

# SO-BC (M546) ---

m546.pbl <- Read10X("../1.Data/gene_expression/BC1_PBMC/filtered_feature_bc_matrix/")
m546.pbl <- readat(m546.pbl, "PBL.IDC")

m546.til <- Read10X("../1.Data/gene_expression/BC1_TIL/filtered_feature_bc_matrix/")
m546.til <- readat(m546.til, "IDC")

# M888 ---

m888.til <- Read10X("../1.Data/gene_expression/M888_Breast_TIL/filtered_feature_bc_matrix/")
m888.til <- readat(m888.til, "IDC")

# M936 (cite, wrong) ---

# m936.data <- Read10X("../1.Data/gene_expression/TP19_M936_IDC/filtered_feature_bc_matrix/")
# 
# m936.cite <- Read10X("../1.Data/citeseq/TP19-M936_../1.Data/citeseq/umi_count/",gene.col=1)
# 
# m936.til.pbmc <- unmt(m936.data, m936.cite, 4, "IDC", 5, "PBMC.IDC")

# M951 (cite, wrong) ---

# m951.data <- Read10X("../1.Data/gene_expression/TP19_M951_IDC/filtered_feature_bc_matrix/")
# 
# m951.cite <- Read10X("../1.Data/citeseq/TP19-M951_../1.Data/citeseq/umi_count/", gene.column = 1)
# 
# m951.til.pbmc <- unmt(m951.data, m951.cite, 4, "IDC", 5, "PBMC.IDC")

# M108 (cite) ---

m108.dat <- Read10X("../1.Data/gene_expression/M108_IDC/filtered_feature_bc_matrix/")

m108.cite <- Read10X("../1.Data/citeseq/M108_IDC_CITEseq/umi_count/", gene.column = 1)

m108.til.pbl <- unmt(m108.dat, m108.cite, 1, "IDC", 2, "PBL.IDC")

# M194 (cite) ---

m194.data <- Read10X("../1.Data/gene_expression/M194_IDC/filtered_feature_bc_matrix/")

m194.cite <- Read10X("../1.Data/citeseq/M194_IDC_CITEseq/umi_count/", gene.column = 1)

m194.pbl.til <- unmt(m194.data, m194.cite, 3, "PBL.IDC", 2, "IDC")

# M226 (cite) ---

m226.data <- Read10X("../1.Data/gene_expression/TP20_M226_ILC/filtered_feature_bc_matrix/")

m226.cite <- Read10X("../1.Data/citeseq/TP20_M226_ILC_CITEseq/umi_count/", gene.column = 1)

m226.til.pbl <- unmt(m226.data, m226.cite, 1, "ILC", 2, "PBL.ILC")

# M234 (cite) ---

m234.data <- Read10X("../1.Data/gene_expression/TP20_M234_ILC/filtered_feature_bc_matrix/")

m234.cite <- Read10X("../1.Data/citeseq/TP20_M234_ILC_CITEseq/umi_count/", gene.column = 1)

m234.til.pbl <- unmt(m234.data, m234.cite, 1, "ILC", 2, "PBL.ILC")

# M268 (cite) ---

m268.data <- Read10X("../1.Data/gene_expression/TP20_M268_IDC_21.02.05/filtered_feature_bc_matrix/")

m268.cite <- Read10X("../1.Data/citeseq/TP20_M268_IDC_Cell_Hash/umi_count/", gene.column = 1)

m268.til.pbl <- unmt(m268.data, m268.cite, 1, "IDC", 4, "PBL.IDC")

# M302 (cite) ---

m302.data <- Read10X("../1.Data/gene_expression/TP20_M302_IDC/filtered_feature_bc_matrix/")

m302.cite <- Read10X("../1.Data/citeseq/TP20_M302_IDC_CITEseq/umi_count/", gene.column = 1)

m302.til.pbl <- unmt(m302.data, m302.cite, 1, "IDC", 2, "PBL.IDC")

# M365 (cite) ---

m365.data <- Read10X("../1.Data/gene_expression/TP20_M365_ILC/filtered_feature_bc_matrix/")

m365.cite <- Read10X("../1.Data/citeseq/TP20_M365_ILC_CITEseq/umi_count/", gene.column = 1)

m365.til.pbl <- unmt(m365.data, m365.cite, 2, "ILC", 3, "PBL.ILC")

# M272 (cite) ---

m272.data <- Read10X("../1.Data/gene_expression/TP20_M272_ILC/filtered_feature_bc_matrix/")

m272.cite <- Read10X("../1.Data/citeseq/TP20_M272_ILC_Cell_Hash/umi_count/", gene.column =  1)

m272.til.pbl <- unmt(m272.data, m272.cite, 1, "ILC", 4, "PBL.ILC")

# M349 (cite) ---

m349.data <- Read10X("../1.Data/gene_expression/TP20_M349_SO_BC2/filtered_feature_bc_matrix/")

m349.cite <- Read10X("../1.Data/citeseq/TP20_M349_SO_BC2_Cell_Hash/umi_count/", gene.column =  1)

m349.til.pbl <- unmt(m349.data, m349.cite, 1, "ILC", 2, "PBL.ILC")

# M571 ---

m571.pbl <- Read10X("../1.Data/gene_expression/TP20_M571_ILC/filtered_feature_bc_matrix/")
m571.pbl <- readat(m571.pbl, "PBL.ILC")

# M533 ---

m533.pbl <- Read10X("../1.Data/gene_expression/TP20_M533_ILC/filtered_feature_bc_matrix/")
m533.pbl <- readat(m533.pbl, "PBL.ILC")

# M107 ---

m107.data <- Read10X("../1.Data/gene_expression/TP20_M107_ILC/filtered_feature_bc_matrix/")

m107.cite <- Read10X("../1.Data/citeseq/TP20_M107_ILC_CITEseq/umi_count/", gene.column =  1)
m107.cite <-  m107.cite[grep("h_totalseq*", rownames(m107.cite)),]
m107.cite@Dimnames[[1]] <- gsub("h_", "", m107.cite@Dimnames[[1]])

m107.til.pbl <- unmt(m107.data, m107.cite, 1, "ILC", 2, "PBL.ILC")

# Integretion ----

rm(list = ls()[- grep("pbl|til",ls())])
projs <- ls()
id = unlist(strsplit(ls(), "\\."))

all.data <- merge(eval(parse(text = projs[1])), 
                  y=sapply(projs[-1], function(x) eval(parse(text = x))), 
                  add.cell.ids = id[grep("^m", id)],
                  project = "ALLCELLS")

rm(list = ls()[which(ls() != "all.data")])
nrow(all.data@meta.data)
table(all.data$reps)
table(all.data$samples)

DefaultAssay(all.data) <- "RNA"

all.data$log10GenesPerUMI <- log10(all.data$nFeature_RNA) / log10(all.data$nCount_RNA)

all.data$mitoRatio <- PercentageFeatureSet(object = all.data, pattern = "^MT-")
all.data$mitoRatio <- all.data@meta.data$mitoRatio / 100

# cell counts ----
png("../3.Result/QC_1_cell_number_distribution.png", width = 8, height = 5, units = "in", res = 300)
all.data@meta.data %>% 
  ggplot(aes(x=reps, fill=samples)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# UMI counts (transcripts) per cell ----
png("../3.Result/QC_2_UMI_counts_distribution.png", width = 8, height = 5, units = "in", res = 300)
all.data@meta.data %>% 
  ggplot(aes(color=samples, x=nCount_RNA, fill= samples)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 350)
dev.off()

# Genes detected per cell ----
png("../3.Result/QC_3_genes_per_cell.png", width = 8, height = 5, units = "in", res = 300)
all.data@meta.data %>% 
  ggplot(aes(color=samples, x=nFeature_RNA, fill= samples)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 180)
dev.off()

all.data@meta.data %>% 
  ggplot(aes(y=log10(nFeature_RNA), x= samples, fill= samples)) + 
  geom_boxplot(alpha = 0.2) + 
  # scale_x_log10() + 
  theme_classic() +
  ylab("lgo10(nGene)") +
  geom_vline(xintercept = 300)

# UMIs vs. genes detected ----
all.data@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~samples)

# Mitochondrial counts ratio ----
png("../3.Result/QC_5_Mitochondrial_gene_ratio.png", width = 8, height = 5, units = "in", res = 300)
all.data@meta.data %>% 
  ggplot(aes(color=samples, x=mitoRatio, fill=samples)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = c(0.003))
dev.off()

# all.data@meta.data ----
all.data@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = samples, fill=samples)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.79)

filtered.data <- subset(x = all.data, 
                        subset = nCount_RNA >= 350 & 
                          nFeature_RNA >= 180 & 
                          log10GenesPerUMI >= 0.79 & 
                          mitoRatio >= 0.0030 )

counts <- GetAssayData(object = filtered.data, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered.counts <- counts[keep_genes, ]
filtered.data <- CreateSeuratObject(filtered.counts, meta.data = filtered.data@meta.data)

rm(counts, filtered.counts, nonzero, keep_genes)

png("../3.Result/QC_6_cell_number_distribution_after_QC.png", width = 8, height = 5, units = "in", res = 300)
filtered.data@meta.data %>% 
  ggplot(aes(x=reps, fill=samples)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
dev.off()

# Batch removal ----

# Cell cycle checking ---
data.phase <- NormalizeData(filtered.data)

data.phase <- CellCycleScoring(data.phase, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)

data.phase <- FindVariableFeatures(data.phase, 
                                   selection.method = "vst",
                                   nfeatures = 2000, 
                                   verbose = FALSE)

data.phase <- ScaleData(data.phase)

data.phase <- RunPCA(data.phase)

png("../3.Result/QC_7_cell_cycle_PCA.png", width = 12, height = 5, units = "in", res = 300)
DimPlot(data.phase,
        reduction = "pca",
        group.by = "Phase",
        split.by = "Phase")
dev.off()

# Batch removal ---
split.data <- SplitObject(filtered.data, split.by = "reps")

for (i in 1:length(split.data)) {
  split.data[[i]] <- NormalizeData(split.data[[i]], verbose = TRUE)
  split.data[[i]] <- CellCycleScoring(split.data[[i]], g2m.features=cc.genes$g2m.genes, s.features=cc.genes$s.genes)
  split.data[[i]] <- SCTransform(split.data[[i]], vars.to.regress = c("mitoRatio", "S.Score", "G2M.Score"))
}


integ.features <- SelectIntegrationFeatures(object.list = split.data, 
                                            nfeatures = 3000) 

# options(future.globals.maxSize= 40*1024^3)
split.data <- PrepSCTIntegration(object.list = split.data, 
                                 anchor.features = integ.features)


integ.anchors <- FindIntegrationAnchors(object.list = split.data, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ.features)

gc()
# memory.limit(654800)

data.integrated <- IntegrateData(anchorset = integ.anchors, 
                                 normalization.method = "SCT")

# memory.limit(size = NA)

DefaultAssay(data.integrated) <- "integrated"

data.integrated <- ScaleData(data.integrated)

data.integrated <- RunPCA(object = data.integrated)

png("../3.Result/Init_1_PCA_Elbow_plot.png", width = 15, height = 5, units = "in", res = 300)
ElbowPlot(data.integrated, ndims = 50)
dev.off()

# Determine percent of variation associated with each PC
pct <- data.integrated[["pca"]]@stdev / sum(data.integrated[["pca"]]@stdev) * 100
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

png("../3.Result/Init_2_PCA_pcs_selection.png", width = 15, height = 5, units = "in", res = 300)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

rm(plot_df, co1, co2, cumu, pct)

png("../3.Result/Init_3_PCA_plot.png", width = 5, height = 5, units = "in", res = 300)
PCAPlot(data.integrated)# , split.by = "reps")
dev.off()

png("../3.Result/Init_4_PCA_plot_splited_by_replicates.png", width = 40, height = 10, units = "in", res = 300)
PCAPlot(data.integrated , split.by = "reps", ncol = 8)
dev.off()

png("../3.Result/Init_5_PCA_plot_splited_by_sample_types.png", width = 20, height = 5, units = "in", res = 300)
PCAPlot(data.integrated , split.by = "samples", ncol = 8)
dev.off()

data.integrated <- RunUMAP(data.integrated, 
                           dims = 1:pcs,
                           reduction = "pca")

DimPlot(data.integrated, split.by = "reps")

data.integrated <- FindNeighbors(data.integrated, dims = 1:pcs)

data.integrated <- FindClusters(data.integrated, resolution = 0.5)

png("../3.Result/Init_6_UMAP_plot(res0.5).png", width = 5, height = 5, units = "in", res = 300)
DimPlot(data.integrated, label = T)+NoLegend()
dev.off()

png("../3.Result/Init_7_UMAP_plot(res0.5)_splited_by_replicates.png", width = 40, height = 10, units = "in", res = 300)
DimPlot(data.integrated, label = T, split.by = "reps", ncol= 8)+NoLegend()
dev.off()

png("../3.Result/Init_8_UMAP_plot(res0.5)_splited_by_sample_types.png", width = 20, height = 5, units = "in", res = 300)
DimPlot(data.integrated, label = T, split.by = "samples", ncol= 8)+NoLegend()
dev.off()

saveRDS(data.integrated, "1.initialDat.rds")
save.image(file = "1.Initialization.RData")

gc()
