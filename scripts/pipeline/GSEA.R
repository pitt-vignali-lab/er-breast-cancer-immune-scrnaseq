library(tidyverse)
library(msigdbr)
library(fgsea)
library(viridis)

de <- read.csv("../3.Result/DE_TIL_ILC2vsILC4.csv") %>% 
  dplyr::select(X, avg_log2FC) %>%
  arrange(desc(avg_log2FC)) %>% 
  na.omit() %>% 
  deframe()

for(cat in paste0("C", c(2, 5, 7))){
  m_df <- msigdbr(species = "Homo sapiens", category = cat)
  
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
  
  fgseaRes <- fgsea(pathways = fgsea_sets,
                    stats = de,
                    minSize=15,
                    maxSize=500)
  
  data.table::fwrite(fgseaRes, paste0("../3.Result/", "gseaResult_ILC2vsILC4_", cat, ".csv"))
  # dir.create(paste0("../3.Result//gseaFigure_IDC2vsIDC4", i))
  
  selPath <- fgseaRes[padj < 0.05,]$pathway
  gs2id <- fgsea_sets[selPath]
  
  gs2val <- lapply(gs2id, function(id) {
    res <- de[id]
    res <- res[!is.na(res)]
  })
  
  nn <- names(gs2val)
  i <- match(nn, fgseaRes$pathway)
  
  j <- order(fgseaRes[["padj"]][i], decreasing = F)
  len <- sapply(gs2val, length)
  gs2val.df <- data.frame(category = rep(nn, times=len),
                          color = rep(as.data.frame(fgseaRes)[i,"padj"], times=len),
                          value = unlist(gs2val))
  
  colnames(gs2val.df)[2] <- "padj"
  
  gs2val.df <- gs2val.df[gs2val.df$category %in% selPath,]
  
  pdf(paste0("../3.Result/gseaFigure_ILC2vsILC4_", cat, ".pdf"), 
      width = 15, 
      height = length(unique(gs2val.df$category))/5)
  ggplot(gs2val.df, aes(x=value, y=reorder(category, padj), fill=padj)) +
    ggridges::geom_density_ridges(alpha = 0.5) +
    scale_fill_viridis(name = "Adjusted P value",
                       guide=guide_colorbar(reverse=T)) +
    # scale_fill_continuous(low="yellow", high="orange", name = "padj", 
    #                       guide=guide_colorbar(reverse=T)) +
    # facet_grid(Cat~., switch = "y", scales = "free_y", space = "free_y") +
    xlab("Gene enrichment") + 
    ylab(NULL)  +
    theme_bw() +
    theme(text = element_text(size=12, face = "bold"),
          legend.position = "bottom",
          legend.key.width = unit(2, 'cm'),
          strip.background = element_rect(colour = "black", fill = "white"))
  dev.off()
}


