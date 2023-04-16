library(tidyverse)
library(stringr)

dfAll2 <- read.csv("../3.Result/new_SPIA_KEGG_Clst2_IDC_vs_ILC.csv")
dfAll3 <- read.csv("../3.Result/new_SPIA_KEGG_Clst3_IDC_vs_ILC.csv")
# cluster 2 IDC up ----
df <- dfAll2[dfAll2$Name %in% strsplit("NOD-like receptor signaling pathway
PD-L1 expression and PD-1 checkpoint pathway in cancer
Cell cycle
Apoptosis
Lipid and atherosclerosis
Th17 cell differentiation
Cytosolic DNA-sensing pathway
NF-kappa B signaling pathway
Regulation of actin cytoskeleton
Toll-like receptor signaling pathway
Focal adhesion
Pathways in cancer
PI3K-Akt signaling pathway
T cell receptor signaling pathway
Chemokine signaling pathway
Antigen processing and presentation
TGF-beta signaling pathway", split = "\n")[[1]],]

df %>%
  arrange(tA) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=tA) ) +
  geom_bar(stat="identity", fill = "#005493") +
  geom_text(hjust = 1, size = 5.5, fontface = "bold", 
            aes(y = 0, label = paste(Name," "))) +
  # theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  coord_polar(theta = "y", clip = "off") +
  ylim(0, ceiling(max(abs(range(df$tA)))) /0.75) +
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "gray90", size = 15),
        panel.grid = element_line(colour = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"))+
     expand_limits(x = c(0, 0))

ggsave(paste0("../3.Result/", Sys.Date(), "_new_spia_circle_bar_IDC_2_up.pdf"), 
       width = 13, height = 8, units = "in",
       dpi = 300, bg = "white")

# cluster 2 IDC dn ----
df <- dfAll2[dfAll2$Name %in% strsplit("Th1 and Th2 cell differentiation
Human immunodeficiency virus 1 infection
Wnt signaling pathway
Endocrine resistance
Proteoglycans in cancer
Fc epsilon RI signaling pathway
Apoptosis - multiple species
Non-alcoholic fatty liver disease
Estrogen signaling pathway
Central carbon metabolism in cancer
HIF-1 signaling pathway
Complement and coagulation cascades
Acute myeloid leukemia
Regulation of lipolysis in adipocytes
Fc gamma R-mediated phagocytosis
Endocrine and other factor-regulated calcium reabsorption
Protein processing in endoplasmic reticulum", split = "\n")[[1]],]

df %>%
  arrange(desc(tA)) %>%
  # filter(tA < 0) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=-tA) ) +
  geom_bar(stat="identity", fill = "#942193") +
  geom_text(hjust = 0, size = 5.5, fontface = "bold", 
            aes( y = 0, label = paste(" ", Name))) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  coord_polar(theta = "y", clip = "off") +
  ylim(ceiling(max(abs(range(df$tA)))) /0.74, 0) +
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "gray90", size = 15),
        panel.grid = element_line(colour = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"))+
  expand_limits(x = c(0, 0))
ggsave(paste0("../3.Result/", Sys.Date(), "_new_spia_circle_bar_IDC_2_dn.pdf"), 
       width = 13, height = 8, units = "in",
       dpi = 300, bg = "white")

# cluster 3 IDC up ----
df <- dfAll3[dfAll3$Name %in% strsplit("NOD-like receptor signaling pathway
HIF-1 signaling pathway
PI3K-Akt signaling pathway
Leukocyte transendothelial migration
Central carbon metabolism in cancer
Fc gamma R-mediated phagocytosis
Natural killer cell mediated cytotoxicity
NF-kappa B signaling pathway
RIG-I-like receptor signaling pathway
Endocrine resistance
Fc epsilon RI signaling pathway
Neutrophil extracellular trap formation
Toll-like receptor signaling pathway
Cytosolic DNA-sensing pathway
TGF-beta signaling pathway
Vascular smooth muscle contraction
Regulation of lipolysis in adipocytes", split = "\n")[[1]],]

df %>%
  arrange(tA) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=tA) ) +
  geom_bar(stat="identity", fill = "#005493") +
  geom_text(hjust = 1, size = 5.5, fontface = "bold", 
            aes(y = 0, label = paste(Name," "))) +
  # theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  coord_polar(theta = "y", clip = "off") +
  ylim(0, ceiling(max(abs(range(df$tA)))) /0.77) +
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "gray90", size = 15),
        panel.grid = element_line(colour = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"))+
  expand_limits(x = c(0, 0))

ggsave(paste0("../3.Result/", Sys.Date(), "_new_spia_circle_bar_IDC_3_up.pdf"), 
       width = 13, height = 8, units = "in",
       dpi = 300, bg = "white")

# cluster 3 IDC dn ----
df <- dfAll3[dfAll3$Name %in% strsplit("Th1 and Th2 cell differentiation
Apoptosis
Staphylococcus aureus infection
Human immunodeficiency virus 1 infection
Chemokine signaling pathway
Growth hormone synthesis, secretion and action
PD-L1 expression and PD-1 checkpoint pathway in cancer
Non-alcoholic fatty liver disease
JAK-STAT signaling pathway
Cytokine-cytokine receptor interaction
Choline metabolism in cancer
Complement and coagulation cascades
Focal adhesion
Cell cycle
Wnt signaling pathway
mTOR signaling pathway
Transcriptional misregulation in cancer", split = "\n")[[1]],]

df %>%
  arrange(desc(tA)) %>%
  # filter(tA < 0) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=-tA) ) +
  geom_bar(stat="identity", fill = "#942193") +
  geom_text(hjust = 0, size = 5.5, fontface = "bold", 
            aes( y = 0, label = paste(" ", Name))) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.text = element_blank()
  ) +
  xlab("") +
  ylab("") +
  coord_polar(theta = "y", clip = "off") +
  ylim(ceiling(max(abs(range(df$tA)))) /0.75, 0) +
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "gray90", size = 15),
        panel.grid = element_line(colour = "gray90"),
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(0,0,0,0),"cm"))+
  expand_limits(x = c(0, 0))
ggsave(paste0("../3.Result/", Sys.Date(), "_new_spia_circle_bar_IDC_3_dn.pdf"), 
       width = 13, height = 8, units = "in",
       dpi = 300, bg = "white")
