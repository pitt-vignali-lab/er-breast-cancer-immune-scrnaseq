library(tidyverse)
pdf.options(encoding='ISOLatin2.enc')

idc <- read.csv("../3.Result/2022-06-23_List_of _receptor_ligand_IDC.csv") %>% 
  separate(interaction_pairs, c("typeA", "typeB"), sep = "_") %>% 
  separate(interaction, c("receptor", 'ligand'), sep = "_") %>% 
  arrange(p_val)

idc <- idc[idc$typeA != idc$typeB,]
idc <- idc[idc$receptor != idc$ligand,]
idc$type <- ifelse(idc$typeA %in% c("CD4 T cells", "CD8 T cells"),
                   "rl", "lr")

idc$cd <- paste(idc$typeA, idc$receptor, sep = "::")
idc$mono <- paste(idc$typeB, idc$ligand, sep = "::")

idc[c(9, 10)] <- t(apply(idc[c(9, 10)], 1,sort))

idc$mono <- gsub("monocytes and macs::", "", idc$mono)

idc$cd <- paste0(idc$cd, ifelse(idc$type == "lr", ">", "<"))
idc$mono <- paste0(idc$mono, ifelse(idc$type == "lr", "<", ">"))
idc$cd <- factor(idc$cd, levels = sort(unique(idc$cd), decreasing = T))

idc <- idc %>% separate(cd, into = c("cell", "lr"), sep = "::",
                        remove = F)

ggplot(idc, aes(mono, lr)) +
  geom_point(shape = 21, aes(size = interact_ratio,
                             alpha = -log10(p_val+1E-5)),
             fill = "#005493") +
  geom_point(shape = 21, aes(size = interact_ratio), 
             col = "black", stroke = 1) +
  facet_grid(cell~., scales = "free_y", space = "free") +
  scale_color_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
  scale_size(breaks = seq(0, 20, 4), limits =  c(1, 16))+
  # coord_fixed() +
  labs(x = "",
       y = "", 
       alpha = "Significance \n(-logPvalue)")+
  theme_bw() +
  theme(legend.position = 'none',
          axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5,
                                   face = 2,
                                   size = 15),
        axis.text.y = element_text(face = 2,
                                   size = 15),
        axis.title = element_text(face = 2,
                                  size = 15),
        legend.text = element_text(face = 2,
                                   size = 15),
        legend.title = element_text(face = 2,
                                    size = 15),
        strip.text = element_text(face = 2,
                                  size = 15),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(0, "lines")) 

ggsave("../3.Result/ligandReceptorCDtoMonoIDC.pdf", width = 7, height = 7, units = "in")

ilc <- read.csv("../3.Result/2022-06-23_List_of _receptor_ligand_ILC.csv")%>% separate(interaction_pairs, c("typeA", "typeB"), sep = "_") %>% 
  separate(interaction, c("receptor", 'ligand'), sep = "_") %>% 
  arrange(p_val)

ilc <- ilc[ilc$typeA != ilc$typeB,]
ilc <- ilc[ilc$receptor != ilc$ligand,]
ilc$type <- ifelse(ilc$typeA %in% c("CD4 T cells", "CD8 T cells"),
                   "rl", "lr")

ilc$cd <- paste(ilc$typeA, ilc$receptor, sep = "::")
ilc$mono <- paste(ilc$typeB, ilc$ligand, sep = "::")

ilc[c(9, 10)] <- t(apply(ilc[c(9, 10)], 1,sort))

ilc$mono <- gsub("monocytes and macs::", "", ilc$mono)

ilc$cd <- paste0(ilc$cd, ifelse(ilc$type == "lr", ">", "<"))
ilc$mono <- paste0(ilc$mono, ifelse(ilc$type == "lr", "<", ">"))
ilc$cd <- factor(ilc$cd, levels = sort(unique(ilc$cd), decreasing = T))

ilc <- ilc %>% separate(cd, into = c("cell", "lr"), sep = "::",
                        remove = F)

ggplot(ilc, aes(mono, lr)) +
  geom_point(shape = 21, aes(size = interact_ratio,
                             alpha = -log10(p_val+1E-5)),
             fill = "#942193") +
  geom_point(shape = 21, aes(size = interact_ratio), 
             col = "black", stroke = 1) +
  facet_grid(cell~., scales = "free_y", space = "free") +
  scale_color_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
  scale_size(breaks = seq(0, 20, 4), limits =  c(1, 16))+
  # coord_fixed() +
  labs(x = "",
       y = "", 
       alpha = "Significance \n(-logPvalue)")+
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1, 
                                   vjust = 0.5,
                                   face = 2,
                                   size = 15),
        axis.text.y = element_text(face = 2,
                                   size = 15),
        axis.title = element_text(face = 2,
                                  size = 15),
        legend.text = element_text(face = 2,
                                   size = 15),
        legend.title = element_text(face = 2,
                                    size = 15),
        strip.text = element_text(face = 2,
                                  size = 15),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        panel.spacing = unit(0, "lines")) 

ggsave("../3.Result/ligandReceptorCDtoMonoILC.pdf", 
       width = 7, height = 7, units = "in")

