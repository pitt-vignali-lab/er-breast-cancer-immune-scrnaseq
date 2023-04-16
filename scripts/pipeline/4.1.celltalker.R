library(Seurat)
library(tidyverse)
library(celltalker)

# scRNA dataset ----
Dat <- readRDS("1.initialDat(res0.5,0.6).rds")
ClusterId <- c("CD4 T cells", # 0
               "Naive T cells", # 1
               "CD8 T cells", # 2
               "CD8 T cells", # 3
               "NK cells", # 4
               "monocytes and macs", # 5
               "B cells", # 6
               "monocytes and macs", # 7
               "Effector T cells", # 8
               "Tregs", # 9
               "monocytes and macs", # 10
               "DCs", # 11
               "pDCs", # 12
               "monocytes and macs" # 13
)

names(ClusterId) <- levels(Dat)

DatLab <- RenameIdents(Dat, ClusterId)

DimPlot(DatLab) +
  coord_fixed()

datSub <- subset(DatLab, idents = c("CD4 T cells", "CD8 T cells", "monocytes and macs"))

# add sub mono macro cluster ----
# newClst <- read.csv("../1.Data/monomcro_new_cluster_id.csv") %>% 
#   # column_to_rownames(var = "X") %>% 
#   mutate(x = paste("mm", x, sep = " ")) %>% 
#   rename(newClst = x)
# 
# orgClst <- data.frame(Idents(datSub)) %>% 
#   rename(newClst = Idents.datSub.) %>% 
#   rownames_to_column(var = "X")
# 
# newClst <- rbind(orgClst[orgClst$newClst != "monocytes and macs",],
#                  newClst[newClst$newClst %in% c("mm 2", "mm 3"),]) %>% 
#   remove_rownames() %>% 
#   column_to_rownames(var = "X")
# 
# datSub <- AddMetaData(subset(datSub, cells = rownames(newClst)), newClst)
# Idents(datSub) <- "newClst"
# 
# DimPlot(datSub) +
#   coord_fixed()

#---- 

datSub[['cell_types']] <- Idents(datSub)

# LR pair ----
# cell talk DB
refCt <- data.table::fread("../../RefData/LigandReceptorDBs/celltalkDB_human_lr_pair.txt")[,c(2, 3)]
colnames(refCt) <- c("ligand", "receptor")


# Bader DB
refBa <- data.table::fread("../../RefData/LigandReceptorDBs/receptor_ligand_interactions_mitab.txt")[,c(1,2)]
refBa <- unique(refBa)
colnames(refBa) <- c("ligand", "receptor")

# ramilowski DB
refRp <- ramilowski_pairs[c(1, 2)]

## combination 
ref <- rbind(refBa, refCt, refRp)
ref <- unique(data.frame(t(apply(ref,1,sort))))
names(ref) <- c("ligand", "receptor")
ref$pair <- paste(ref$ligand, ref$receptor, sep = "_")

# cell talker ----
DimPlot(datSub)

hca_bm_interactions <- celltalk(input_object=subset(datSub, samples == "ILC"),
                                metadata_grouping="cell_types",
                                ligand_receptor_pairs=ref,
                                number_cells_required=100,
                                min_expression=100,
                                max_expression=20000,
                                scramble_times=10)

## Identify top statistically significant interactions
top_stats <- hca_bm_interactions %>%
  # mutate(fdr=p.adjust(p_val,method= "fdr")) %>%
  filter(p_val < 0.05) %>%
  group_by(cell_type1) %>%
  top_n(10,interact_ratio) %>%
  ungroup()

write.csv(top_stats[c(1, 2)], 
          paste0("../3.Result/", Sys.Date(), "_List_of _receptor_ligand_ILC.csv"))
## Generate a circos plot
colors_use <- RColorBrewer::brewer.pal(n=length(unique(datSub$cell_types)),"Set2")

png("../3.Result/cmbBD_mm_IDC.png", 
    width = 5, height = 5, units = "in", res = 300)
circos_plot(ligand_receptor_frame=top_stats,
            cell_group_colors=colors_use,
            ligand_color="blue",
            receptor_color="red",
            cex_outer=0.5,
            cex_inner=0.4)
dev.off()

ilc <- top_stats$interaction
idc <- top_stats$interaction

intersect(ilc, idc)
ilc[! ilc %in% intersect(ilc, idc)]
idc[! idc %in% intersect(ilc, idc)]
