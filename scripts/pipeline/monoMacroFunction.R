library(Seurat)
library(fgsea)
library(msigdbr)
library(tidyverse)

dat <- readRDS("MonoMacro.RDS")

DimPlot(dat)

dat <- subset(dat, idents = 5, invert = T)

datSub <- subset(dat, samples %in% c( "PBL.IDC", "PBL.ILC"))

sce <- as.SingleCellExperiment(DietSeurat(dat, dimreducs = "pca"))

DefaultAssay(dat) <- "RNA"

# Function analysis ----
m_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")

fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 

dfRes <- data.frame()

for (i in 1:4) {
  gFc <- FoldChange(dat, ident.1 = i)%>% 
    rownames_to_column(var = "X") %>% 
    dplyr::select(X, avg_log2FC) %>%
    arrange(desc(avg_log2FC)) %>% 
    na.omit() %>% 
    deframe()

  fgseaRes <- fgsea(pathways = fgsea_sets,
                    stats = gFc,
                    minSize=15,
                    maxSize=500)
  
  dfRes <- rbind(dfRes, 
                 fgseaRes %>% mutate(clust = i))
}

dfResSub <- dfRes[dfRes$padj < 0.05,]

sort(table(dfResSub$pathway))

data.table::fwrite(dfResSub, "../3.Result/monoClusterKEGG.csv")

