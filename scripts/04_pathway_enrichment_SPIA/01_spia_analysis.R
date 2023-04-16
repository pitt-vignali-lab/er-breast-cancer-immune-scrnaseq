library(clusterProfiler)
library(SPIA)
# library(EnrichmentBrowser)
# library(pathview)
library(org.Hs.eg.db)

# kegg.gs <- getGenesets(org="hsa", db="kegg")
# dir.create("./pathway.dir")
# download.kegg(pathway.id = substr(names(kegg.gs),4,8), 
#               species = "hsa", kegg.dir = "./pathway.dir", file.type=c("xml"))
# makeSPIAdata(kgml.path = "pathway.dir/", organism = "hsa", out.path = "./")

allGene <- unique(as.data.frame(org.Hs.egACCNUM2EG)$gene_id)

deGs <- read.csv("../3.Result/DE_TIL_ILC2vsIDC2.csv")

geneDf <- bitr(deGs$X, fromType = "SYMBOL",
               toType = c("ENTREZID"),
               OrgDb = org.Hs.eg.db)

geneRank <- merge(geneDf, deGs, by.x = "SYMBOL", by.y = "X")

ranks <- geneRank$avg_log2FC
names(ranks) <- geneRank$ENTREZID  
ranks <- sort(ranks, decreasing = T)

spiaRes <- spia(de = ranks, all = allGene, organism = "hsa", data.dir = "./", nB = 2000)

write.csv(spiaRes, "../3.Result/spia_TIL_ILC2vsIDC2.csv")
