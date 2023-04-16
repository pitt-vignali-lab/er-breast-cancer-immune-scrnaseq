library(Seurat)
library(tidyverse)
library(ggrepel)

gList <- read.csv("../3.Result/DE_TIL_IDC2vsILC2.csv")



dat <- readRDS("MonoMacro.RDS")
DefaultAssay(dat) <- "RNA"
dat$des <- paste0(dat$samples, dat$integrated_snn_res.0.2)
Idents(dat) <- "des"

# Clsuter 2 ----
DEgene2 <- FindMarkers(dat, ident.1 = "IDC2", ident.2 = "ILC2")

gnSelect2 <- c(head(rownames(DEgene2[DEgene2$avg_log2FC > 0,]), 1),
              head(rownames(DEgene2[DEgene2$avg_log2FC < 0,]), 1),
              rownames(DEgene2[DEgene2$avg_log2FC == max(DEgene2$avg_log2FC),]),
              rownames(DEgene2[DEgene2$avg_log2FC == min(DEgene2$avg_log2FC),]))

gnSelect2 <- c(gnSelect2,
               strsplit("CCL3L3
                      CCL5
                      IL7R
                      IL8
                      CD36
                      VCAN
                      IL32
                      LTB
                      FCN1
                      MALAT1
                      TNFSF13
                      VSIR
                      CXCL8
                      HLADQA2
                      APOC1
                      CCL3L1
                      CTSB
                      IL18
                      LAIR1
                      CD81", "
                      ")[[1]])

# Cluster 4 ----
DEgene4 <- FindMarkers(dat, ident.1 = "IDC4", ident.2 = "ILC4")

gnSelect4 <- unique(c(head(rownames(DEgene4[DEgene4$avg_log2FC > 0,]), 4),
                     head(rownames(DEgene4[DEgene4$avg_log2FC < 0,]), 5),
                     rownames(DEgene4[DEgene4$avg_log2FC == max(DEgene4$avg_log2FC),]),
                     rownames(DEgene4[DEgene4$avg_log2FC == min(DEgene4$avg_log2FC),])))


p1 <- ggplot(DEgene2, aes(avg_log2FC, -log10(p_val_adj+1e-50), 
                          label = ifelse(rownames(DEgene2) %in% gnSelect2, rownames(DEgene2), ""))) +
  # scale_y_continuous(trans = "log10") +
  geom_point(aes(col = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "UP",
                              ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "DN",
                                     "black"))), 
             alpha = 0.8,
             size = 5, shape = 21, stroke = 2)+
  geom_label_repel(size = 10, max.overlaps = 20) +
  scale_color_manual(name = "Gene expressin status", 
                     values=c("red", "blue"), 
                     labels = c("Up-regulated genes in IDC", 
                                "Up-regulated genes in ILC"),
                     limits = c("UP", "DN")) +
  ggtitle("Monomacro cluster 2 IDC vs ILC") +
  theme_bw()+
  xlab("Log2 folder change") +
  xlim(c(-6, 6))+
  scale_y_continuous(name = "Significnace \n(adjusted p value)",
                     labels = c(0, 1E-20, 1E-40, 1E-60),
                     breaks = c(0, 20, 40, 60))+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20, face="bold"),
        panel.grid = element_blank(),
        axis.text = element_text(size=20, face="bold", colour = "black"),
        axis.ticks = element_line(color="black", size = 2),
        panel.border = element_rect(colour = "black", fill=NA, size=2))

p2 <- ggplot(DEgene4, aes(avg_log2FC, -log10(p_val_adj+1e-50), 
                         label = ifelse(rownames(DEgene4) %in% gnSelect4, rownames(DEgene4), ""))) +
  # scale_y_continuous(trans = "log10") +
  geom_point(aes(col = ifelse(p_val_adj < 0.05 & avg_log2FC > 0, "UP",
                              ifelse(p_val_adj < 0.05 & avg_log2FC < 0, "DN",
                                     "black"))), 
             alpha = 0.8,
             size = 5, shape = 21, stroke = 2)+
  geom_label_repel(size = 10, max.overlaps = 20) +
  scale_color_manual(name = "Gene expressin status", 
                     values=c("red", "blue"), 
                     labels = c("Up-regulated genes in IDC", 
                                "Up-regulated genes in ILC"),
                     limits = c("UP", "DN")) +
  ggtitle("Monomacro cluster 4 IDC vs ILC") +
  theme_bw()+
  xlab("Log2 folder change")+
  xlim(c(-6, 6))+
  scale_y_continuous(name = "Significnace \n(adjusted p value)",
                     labels = c(0, 1E-20, 1E-40, 1E-60),
                     breaks = c(0, 20, 40, 60))+
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20, face="bold"),
        panel.grid = element_blank(),
        axis.text = element_text(size=20, face="bold", colour = "black"),
        axis.ticks = element_line(color="black", size = 2),
        panel.border = element_rect(colour = "black", fill=NA, size=2))

png("../3.Result/volcanoPlotRegular.png", width = 15, height = 25, units = "in", res = 300)
ggpubr::ggarrange(p1, p2, ncol=1, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

#-------------------------------------------------------------------------------

DEgene4$Clst <- 4
DEgene2$Clst <- 2

deG <- merge(DEgene2[DEgene2$p_val_adj < 0.05,][2],
             DEgene4[DEgene4$p_val_adj < 0.05,][2],
             by = 0, all = T)

ggplot() +
  geom_point(data = DEgene2[DEgene2$p_val_adj<0.05,], aes(avg_log2FC, p_val_adj))+
  geom_point(data = DEgene4[DEgene4$p_val_adj<0.05,], aes(p_val_adj, avg_log2FC))
