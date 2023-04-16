library(Seurat)
library(tidyverse)
require(scales)
library(ggsignif)
library(viridis)

# Function defination ----
Db <- function(pca.mat.x,pca.mat.y,pcs.use=1:50,num.sample=300,random.sample=FALSE) {
  
  if (random.sample==FALSE) {
    
    samp1 <- sample(1:nrow(pca.mat.x),num.sample, replace = T)
    samp2 <- sample(1:nrow(pca.mat.y),num.sample, replace = T)
    
    mat.x <- pca.mat.x[samp1,pcs.use]
    mat.y <- pca.mat.y[samp2,pcs.use]
    
    mean.mat1 <- apply(mat.x,2,mean)
    mean.mat2 <- apply(mat.y,2,mean)
    
    cov.mat1 <- cov(mat.x)
    cov.mat2 <- cov(mat.y)
    
    comb.cov <- (cov.mat1 + cov.mat2)/2
    
    (1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))
    
  } else {
    
    comb.mat <- rbind(pca.mat.x,pca.mat.y)
    
    samp1 <- sample(1:nrow(comb.mat),num.sample, replace = T)
    samp2 <- sample(1:nrow(comb.mat),num.sample, replace = T)
    
    mat.x <- comb.mat[samp1,pcs.use]
    mat.y <- comb.mat[samp2,pcs.use]
    
    mean.mat1 <- apply(mat.x,2,mean)
    mean.mat2 <- apply(mat.y,2,mean)
    
    cov.mat1 <- cov(mat.x[,pcs.use])
    cov.mat2 <- cov(mat.y[,pcs.use])
    
    comb.cov <- (cov.mat1 + cov.mat2)/2
    
    (1/8)*t(mean.mat1-mean.mat2)%*%solve(comb.cov)%*%(mean.mat1-mean.mat2)+0.5*log(det(comb.cov)/sqrt(det(cov.mat1)*det(cov.mat2)))
    
    
  }
  
}

# data readin ----
Dat <- readRDS("1.initialDat(res0.5,0.6).rds")
DefaultAssay(Dat) <- "RNA"

ClusterId <- c("CD4 T cells", # 0
               "Naive T cells", # 1
               "CD8 T cells", # 2
               "CD8 T cells", # 3
               "NK cells", # 4
               "Mono&macs", # 5
               "B cells", # 6
               "Mono&macs", # 7
               "Effector T cells", # 8
               "Tregs", # 9
               "Mono&macs", # 10
               "DCs", # 11
               "pDCs", # 12
               "Mono&macs" # 13
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

DatLab$cellType <- DatLab@active.ident
# PCs selecting ----
pct <- Dat[["pca"]]@stdev / sum(Dat[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
pcs <- min(co1, co2)

# cell type selection----

t <- table(DatLab$tissues, Idents(DatLab))

t <- table(DatLab$samples, Idents(DatLab))[c(1,2),]

# ctslct <- colnames(t)[apply(t, 2, function(x) all(x > 500))] # cell type selected

# DatSub <- subset(DatLab, cellType %in% ctslct)

til <- subset(DatLab, tissues == "TIL", 
              idents = c("B cells", "Tregs", "Mono&macs", "CD8 T cells", "CD4 T cells", "NK cells")) 
pbmc <- subset(DatLab, tissues == "PBL", 
               idents = c("B cells", "Tregs", "Mono&macs", "CD8 T cells", "CD4 T cells", "NK cells"))

til <- subset(DatLab, samples == "IDC", 
              idents = c("B cells", "Mono&macs", "CD8 T cells", "CD4 T cells", "NK cells")) 
pbmc <- subset(DatLab, samples == "ILC", 
               idents = c("B cells", "Mono&macs", "CD8 T cells", "CD4 T cells", "NK cells"))

til <- subset(DatLab, samples == "PBL.IDC", 
              idents = c("B cells", "Mono&macs", "CD8 T cells", "CD4 T cells", "Naive T cells")) 
pbmc <- subset(DatLab, samples == "PBL.ILC", 
               idents = c("B cells", "Mono&macs", "CD8 T cells", "CD4 T cells", "Naive T cells"))

pbmc.split <- SplitObject(pbmc,split.by="cellType")
pbmc.embed <- lapply(pbmc.split,function(x) Embeddings(x,reduction="pca"))
til.split <- SplitObject(til,split.by="cellType")
til.embed <- lapply(til.split,function(x) Embeddings(x,reduction="pca"))

bdDf <- data.frame()
Db.samp <- c()
Db.rand <- c()

for (a in unique(pbmc$cellType)) {
  
  for (i in 1:100) {
    
    Db.samp[i] <- Db(pca.mat.x=pbmc.embed[[a]],pca.mat.y=til.embed[[a]],pcs.use=1:pcs,num.sample=500,random.sample=FALSE)
    
    Db.rand[i] <- Db(pca.mat.x=pbmc.embed[[a]],pca.mat.y=til.embed[[a]],pcs.use=1:pcs,num.sample=500,random.sample=TRUE)
    
  }
  
  bdDf <- rbind(bdDf,
                data.frame(celltype = a,
                           value = Db.samp,
                           type = "sampled"),
                data.frame(celltype = a,
                           value = Db.rand,
                           type = "random"))
  
}

bdDf$type <- factor(bdDf$type, level = c("sampled", "random"))


# data.frame(cbind(bdDf[bdDf$type == "sampled",c(1,2)],
#       y = bdDf[bdDf$type == "random",]$value)) %>%
# ggplot(aes(x)) +
#   geom_density(aes(y = value, x = ..density..), col = "green")+
#   geom_density(aes(y = y, x = -..density..), col = "red")+
#   ylim(0, 10) + 
#   facet_wrap(~celltype, nrow = 1)+
#   theme_void()

summaryDf <- cbind(bdDf  %>% 
              group_by(celltype) %>% 
              summarize(results = wilcox.test(value ~ type)$p.value) %>% 
                rename(pval = results),
            
            bdDf  %>% 
              group_by(celltype, type) %>% 
              summarize(results =mean(value)) %>%
              pivot_wider(names_from = type, values_from = results) %>%
              mutate(times.mean = sampled/random) %>%
              select(times.mean),
            
            bdDf  %>%
              group_by(celltype) %>%
              summarize(results = max(value)) %>% 
              rename(max = results),
            
            bdDf  %>%
              group_by(celltype) %>%
              summarize(results = min(value)) %>% 
              rename(min = results), 
            
            bdDf  %>% 
              group_by(celltype, type) %>% 
              summarize(results =median(value)) %>%
              pivot_wider(names_from = type, values_from = results) %>%
              mutate(times.median = sampled/random) %>%
              select(times.median)) %>% 
  select(-c(3, 5, 7, 9)) %>% 
  mutate(FDR = p.adjust(pval, method = "fdr"))
summaryDf

bdDf$celltype <- factor(bdDf$celltype,
                        # levels = summaryDf[order(summaryDf$times.mean, decreasing = T),1])
                        levels = summaryDf[order(summaryDf$max, summaryDf$times.median, decreasing = T),1])

summaryDf$celltype <- factor(summaryDf$celltype,
                             # levels = summaryDf[order(summaryDf$times.mean, decreasing = T),1])
                             levels = summaryDf[order(summaryDf$max, summaryDf$times.median, decreasing = T),1])

# png("../3.Result/Bh_distance_250_TIL(IDCvsILC)1.png", width = 12, height = 5, units = "in", res = 300)
ggplot(bdDf %>% group_by(celltype, type) %>%
         mutate(meanV = median(value)), 
       aes(x = type, y = value)) +
  geom_point(aes(fill = value), shape = 21, 
             position=position_jitterdodge(jitter.width = 1), 
              alpha = 0.3, size = 1.8)+
  geom_boxplot(fill = NA, lwd = 0.8, col = "black") +
  facet_grid(~celltype)+
  coord_fixed(ratio = 3.6/max(bdDf$value)) +
  geom_text(data = summaryDf, aes(x = 1.5, y = max(bdDf$value) + 0.1,
                                  label = paste0(round(times.median, 2), "-fold")),
            fontface = "italic",
            size = 3,
            col = "gray50")+
  # geom_text(data = summaryDf, aes(x = 1.5, y = 2, 
  #                                 label = paste0("p value: ", scientific(results))), size = 3) +
  # stat_signif(comparisons=list(c('sampled','random')),map_signif_level = T, test = wilcox.test)+
  # scale_color_gradient(low=hcl(195,100,75), high=hcl(15,100,75))+ 
  scale_fill_viridis(option = "plasma", # use virdis' plasma scale
                      # begin = 0.2, # limit range of used hues
                      name = "Distance")+
  ggtitle("TILvsPBL")+
  ylim(0, max(bdDf$value) + 0.1)+
  xlab(NULL)+
  ylab("Distance") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 12, face = 2),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 15),
        # strip.text.x = element_blank(),
        strip.background.x = element_blank()
        )
# dev.off()
ggsave(paste0("../3.Result/", Sys.Date(), "_bhatta_dist_TILvsPBL.png"), 
       width = 12, height = 4, units = "in", dpi = 300)
