library(Seurat)
library(tidyverse)
library(msigdbr)
library(fgsea)
library(SPIA)

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

DimPlot(dat)

# DE ---- 
## Method 1
datSub <- subset(dat, samples == "ILC")

deG <- FindAllMarkers(datSub)

write.csv(deG %>% filter(cluster == 2), 
          "../3.Result/new_DeG_ILC_Clst2.csv")

## Method 2
datSub <- subset(dat, idents = 2)
Idents(datSub) <- "samples"
deG <- FindMarkers(datSub, ident.1 = "IDC", ident.2 = "ILC")
write.csv(deG, 
          "../3.Result/new_DeG_Clst2_IDC_vs_ILC.csv")

# GSEA ----
files <- list.files("../3.Result/")[grep("new_DeG", list.files("../3.Result/"))]

m_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

for(file in files){
  de <- read.csv(paste0("../3.Result/", file))
  
  if("gene" %in% colnames(de)){ 
         de = select(de, gene, avg_log2FC)
  }else{
           de = select(de, X, avg_log2FC)}
         
    de <- de %>%
      filter(avg_log2FC != 0) %>% 
    arrange(desc(avg_log2FC)) %>% 
    na.omit() %>% 
    deframe()
  
  fgseaRes <- fgsea(pathways = fgsea_sets,
                    stats = de,
                    minSize=15,
                    maxSize=500)
  
  data.table::fwrite(fgseaRes, paste0("../3.Result/", gsub("DeG", "GO_BP", file)))
}

# SPIA ----
refG <- unique(read.csv("../../RefData/hsaRefGeneGeneIDSymbol.csv"))
allGene <- unique(refG$GeneID)

for(file in files){
  de <- read.csv(paste0("../3.Result/", file))
  
  if("gene" %in% colnames(de)){ 
    de = select(de, gene, avg_log2FC)
  }else{
    de = select(de, X, avg_log2FC) %>% rename(gene = X)}
  
  de <- de %>%
    merge(refG, by.x = "gene", by.y = "Symbol") %>% 
    select(GeneID, avg_log2FC) %>% 
    filter(avg_log2FC != 0) %>% 
    unique() %>% 
    arrange(desc(avg_log2FC)) %>% 
    na.omit() %>% 
    deframe()
  
  spiaRes <- spia(de = de, all = allGene, organism = "hsa", 
                  data.dir = "../../RefData/", nB = 2000)
  
  data.table::fwrite(spiaRes, paste0("../3.Result/", gsub("DeG", "SPIA_KEGG", file)))
}
