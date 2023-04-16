library(Seurat)
library(tidyverse)
require(scales)
library(philentropy)
library(fpc)

Dat <- readRDS("1.initialDat(res0.5,0.6).rds")
DefaultAssay(Dat) <- "RNA"

ClusterId <- c("CD4 T cells", # 0
               "Naive T cells", # 1
               "CD8 T cells", # 2
               "CD8 T cells", # 3
               "NK cells", # 4
               "monocytes and \nmacs", # 5
               "B cells", # 6
               "monocytes and \nmacs", # 7
               "tissue resident \n T cells", # 8
               "Tregs", # 9
               "monocytes and \nmacs", # 10
               "DCs", # 11
               "pDCs", # 12
               "monocytes and \nmacs" # 13
)

names(ClusterId) <- levels(Dat)

DatLab <- RenameIdents(Dat, ClusterId)

tempMeta <-   data.frame(samples = unique(DatLab$samples),
                         tissues = c("TIL", "PBL", "TIL", "PBL"))

newMeta <- DatLab@meta.data %>%
  select("samples") %>%
  rownames_to_column("cells") %>%
  merge(tempMeta, by = "samples") %>%
  column_to_rownames("cells") %>%
  select(-"samples")
  
DatLab <- AddMetaData(DatLab, newMeta)

rm(tempMeta, newMeta)
# ----

rndMatGen <- function(condi){
  
  sub1 <- subset(DatLab, tissues == condi)
  
  sub1 <- NormalizeData(sub1, verbose = F)
  
  sub1 <- FindVariableFeatures(sub1, verbose = F)
  
  sub1 <- ScaleData(sub1, verbose = F)
  
  sub1 <- RunPCA(sub1, features = VariableFeatures(object = sub1), verbose = F)
  
  df <- Embeddings(sub1, reduction = "pca")
  
  return(df)}

matGen <- function(condi, cellTp) {
  
  sub1 <- subset(DatLab, tissues == condi, idents = cellTp)
  
  sub1 <- NormalizeData(sub1, verbose = F)
  
  sub1 <- FindVariableFeatures(sub1, verbose = F)
  
  sub1 <- ScaleData(sub1, verbose = F)
  
  sub1 <- RunPCA(sub1, features = VariableFeatures(object = sub1), verbose = F)
  
  df <- Embeddings(sub1, reduction = "pca")
  
  return(df) }

distDf <- data.frame()

for (i in unique(ClusterId)) {
  mat1 <- matGen("TIL", i)
  mat2 <- matGen("PBL", i)
  mat3 <- rndMatGen("TIL")
  mat4 <- rndMatGen("PBL")
  
  for (j in 1:100) {
    tryCatch({
      sampleNo <- 500
      
      rm(.Random.seed, envir=globalenv())
      subM1 <- mat1[sample(rownames(mat1), sampleNo),]
      
      rm(.Random.seed, envir=globalenv())
      subM2 <- mat2[sample(rownames(mat2), sampleNo),]
      
      rm(.Random.seed, envir=globalenv())
      subM3 <- mat3[sample(rownames(mat3), sampleNo),]
      
      rm(.Random.seed, envir=globalenv())
      subM4 <- mat4[sample(rownames(mat4), sampleNo),]
      
      d1 <- bhattacharyya.dist(apply(subM1, 1, mean), apply(subM2, 1, mean),
                               cov(t(as.matrix(subM1))), cov(t(as.matrix(subM2))))
      
      d2 <- bhattacharyya.dist(apply(subM3, 1, mean), apply(subM4, 1, mean),
                               cov(t(as.matrix(subM3))), cov(t(as.matrix(subM4))))
      
      
      # cat(d1, d2, "\n")
      
      distDf <- rbind(distDf, 
                      data.frame(sample = "TIL vs PBL",
                                 Vlu = d1,
                                 type = "True",
                                 celltp = i,
                                 run = j),
                      data.frame(sample = "TIL vs PBL",
                                 Vlu = d2,
                                 type = "Random",
                                 celltp = i,
                                 run = j))
    }, error=function(e){})
  }
  
  cat(i, ": Done \n")
}
distDf$type <- factor(distDf$type, levels = c("True", "Random"))

png("../3.Result/2.clusterLabeling_bhattacharyya_dist_3.png", width = 10, height = 8, units = "in", res = 300)
ggplot(distDf) + 
  geom_boxplot(aes(x = celltp, y = Vlu, fill = type)) +
  facet_wrap(~sample, ncol = 1)+ 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()
