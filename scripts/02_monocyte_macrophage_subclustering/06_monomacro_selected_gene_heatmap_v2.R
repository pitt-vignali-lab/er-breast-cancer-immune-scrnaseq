library(tidyverse)
library(ComplexHeatmap)

# cluster 2 ----
df1 <- read.csv("../3.Result/new_DeG_IDC_Clst2.csv") %>% 
  filter(gene %in% tarG[1:85]) %>% 
  select(gene, avg_log2FC) %>% 
  rename("IDC" = 2)

df2 <- read.csv("../3.Result/new_DeG_ILC_Clst2.csv") %>% 
  filter(gene %in% tarG[1:85]) %>% 
  select(gene, avg_log2FC) %>% 
  rename("ILC" = 2 )

df <- merge(df1, df2, by = "gene", all = T) %>% 
  column_to_rownames(var = "gene") %>% 
  as.matrix()

df[is.na(df)] <- 0

dfTrans <- apply(df, 2, function(x) (x - min(x))/(max(x) - min(x)))

htmp <- unique(dfTrans)

htGenes <- strsplit("IRF7
MYOF
CEBPD
SIRPG
TNFSF10
MS4A6A
SLC11A1
CLEC10A
JUN
IL10RA
RAB24
WDR54
BCL2L13
TIMP1
PLSCR1
DOCK7
TGFBI
LST1
FCER1G
MILR1
CSTA", split = "\n")[[1]]

ha = rowAnnotation(foo = anno_mark(at = which(rownames(htmp) %in%  htGenes), 
                                   labels_gp = gpar(fontsize = 10, fontface = 2),
                                   labels = rownames(htmp)[rownames(htmp) %in%  htGenes]))

pdf("../3.Result/mono_macor_re-cluster_selected_genes_in_cluster2.pdf", 
    width = 4,height = 8.5)
draw(Heatmap(htmp,
             na_col = "black",
             rect_gp = gpar(col = "gray50", lwd = 0.8),
             row_dend_width = unit(4, "cm"),
             col = circlize::colorRamp2(c(0, 0.5, 1), c("yellow", "gray", "red")),
             right_annotation = ha, 
             row_names_gp = gpar(fontsize = 8),
             cluster_rows = T,
             show_row_dend = T,
             cluster_columns = F,
             column_title = NULL,
             row_names_side = "left",
             column_names_rot = 0,
             column_names_centered = T,
             column_names_gp = gpar(fontsize = 10, fontface = 2),
             name = "Normalized Gene Expression",
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter",
                                         legend_width = unit(5, "cm"))), 
     heatmap_legend_side = "bottom")

dev.off()

# Heatmap(unique(dfTrans),
#         na_col = "black",
#         rect_gp = gpar(col = "gray50", lwd = 0.8),
#         row_dend_width = unit(2, "cm"),
#         col = circlize::colorRamp2(c(0, 0.5, 1), c("yellow", "gray", "red")),
#         name = "Normalized logFC",
#         # col = col_fun,
#         # right_annotation = ha,
#         row_names_gp = gpar(fontsize = 7.5),
#         cluster_rows = T,
#         show_row_dend = T,
#         cluster_columns = F)

# Cluster 3 ----
df1 <- read.csv("../3.Result/new_DeG_IDC_Clst3.csv") %>% 
  filter(gene %in% tarG[86:145]) %>% 
  select(gene, avg_log2FC) %>% 
  rename("IDC" = 2)

df2 <- read.csv("../3.Result/new_DeG_ILC_Clst3.csv") %>% 
  filter(gene %in% tarG[86:145]) %>% 
  select(gene, avg_log2FC) %>% 
  rename("ILC" = 2 )

df <- merge(df1, df2, by = "gene", all = T) %>% 
  column_to_rownames(var = "gene") %>% 
  as.matrix()

df[is.na(df)] <- 0

dfTrans <- apply(df, 2, function(x) (x - min(x))/(max(x) - min(x)))

htmp <- unique(dfTrans)

htGenes <- strsplit("VCAN
KLF10
STXBP2
PLBD1
SYNJ1
LGALS3
ARFIP1
PIK3IP1
CD74
CXCR3
CXCR4
FTH1
NEAT1
GZMK
RETN
CTSB", split = "\n")[[1]]

ha = rowAnnotation(foo = anno_mark(at = which(rownames(htmp) %in%  htGenes), 
                                   labels_gp = gpar(fontsize = 10, fontface = 2),
                                   labels = rownames(htmp)[rownames(htmp) %in%  htGenes]))

pdf("../3.Result/mono_macor_re-cluster_selected_genes_in_cluster3.pdf", 
    width = 4,height = 6)
draw(Heatmap(htmp,
             na_col = "black",
             rect_gp = gpar(col = "gray50", lwd = 0.8),
             row_dend_width = unit(4, "cm"),
             col = circlize::colorRamp2(c(0, 0.5, 1), c("yellow", "gray", "red")),
             right_annotation = ha, 
             row_names_gp = gpar(fontsize = 8),
             cluster_rows = T,
             show_row_dend = T,
             cluster_columns = F,
             column_title = NULL,
             row_names_side = "left",
             column_names_rot = 0,
             column_names_centered = T,
             column_names_gp = gpar(fontsize = 10, fontface = 2),
             name = "Normalized Gene Expression",
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter",
                                         legend_width = unit(5, "cm"))), 
     heatmap_legend_side = "bottom")

dev.off()
