library(tidyverse)
library(ggpubr)
# ----
ptw <- strsplit("Pathways upregulated in ER+ IDC cluster 4
Antigen processing and presentation
Chemokine signaling pathway
Transcriptional misregulation in cancer
Toll-like receptor signaling pathway
B cell receptor signaling pathway
Fc gamma R-mediated phagocytosis
Natural killer cell mediated cytotoxicity
PD-L1 expression and PD-1 checkpoint pathway in cancer
Leukocyte transendothelial migration
Estrogen signaling pathway
Human immunodeficiency virus 1 infection
mTOR signaling pathway
Phospholipase D signaling pathway
Staphylococcus aureus infection
Cell cycle
Viral protein interaction with cytokine and cytokine receptor
Endocrine resistance
C-type lectin receptor signaling pathway
Tuberculosis
Neutrophil extracellular trap formation
Pathways upregulated in ER+ ILC cluster 4
Th17 cell differentiation
Th1 and Th2 cell differentiation
Apoptosis
NOD-like receptor signaling pathway
T cell receptor signaling pathway
TNF signaling pathway
Cholesterol metabolism
Adipocytokine signaling pathway
Cytosolic DNA-sensing pathway
JAK-STAT signaling pathway
Pathways in cancer
Thermogenesis
Cholinergic synapse
p53 signaling pathway
Central carbon metabolism in cancer
ErbB signaling pathway
Hepatitis B
Hepatitis C
Fc epsilon RI signaling pathway
Cellular senescence
Pathways upregulated in ER+ IDC cluster 2
Antigen processing and presentation
Lipid and atherosclerosis
Toll-like receptor signaling pathway
Th17 cell differentiation
Cytokine-cytokine receptor interaction
Chemokine signaling pathway
Leukocyte transendothelial migration
NF-kappa B signaling pathway
MAPK signaling pathway
Osteoclast differentiation
Human immunodeficiency virus 1 infection
Kaposi sarcoma-associated herpesvirus infection
Non-alcoholic fatty liver disease
Focal adhesion
Viral protein interaction with cytokine and cytokine receptor
Renal cell carcinoma
Pathogenic Escherichia coli infection
Renal cell carcinoma
Pertussis
Focal adhesion
Pathways upregulated in ER+ ILC cluster 2
NOD-like receptor signaling pathway
HIF-1 signaling pathway
Apoptosis
PD-L1 expression and PD-1 checkpoint pathway in cancer
Fc gamma R-mediated phagocytosis
IL-17 signaling pathway
Cholesterol metabolism
T cell receptor signaling pathway
Th1 and Th2 cell differentiation
TNF signaling pathway
Human immunodeficiency virus 1 infection
Kaposi sarcoma-associated herpesvirus infection
Non-alcoholic fatty liver disease
Focal adhesion
Viral protein interaction with cytokine and cytokine receptor
Renal cell carcinoma
Pathogenic Escherichia coli infection
Renal cell carcinoma
Pertussis
Focal adhesion", "\n")[[1]]

# ----

ptwDf <- data.frame()

titleOd <- grep("Pathways upregulated", ptw)

for(i in 1:length(titleOd)){
    if (i == length(titleOd)) {
      t <- data.frame(Cat = ptw[titleOd[i]],
                      pthw = ptw[(titleOd[i]+1) : length(ptw)])
    } else {
      t <- data.frame(Cat = ptw[titleOd[i]],
                      pthw = ptw[(titleOd[i]+1) : (titleOd[i+1]-1)])
    }
  
  ptwDf <- rbind(t, ptwDf)
}

spia4 <- read.csv("../3.Result/spia_TIL_IDC4vsILC4.csv")

spia4 <- spia4[spia4$Name %in% ptwDf[grep(4, ptwDf$Cat), ]$pthw,] %>% 
  merge(ptwDf[grep(4, ptwDf$Cat), ], by.x = "Name", by.y = "pthw") %>% 
  select(c("Name", "NDE", "tA", "pGFdr", "Cat", "Status"))

spia4$Cat <- 4

spia2 <- read.csv("../3.Result/spia_TIL_IDC2vsILC2.csv")

spia2 <- spia2[spia2$Name %in% ptwDf[grep(2, ptwDf$Cat), ]$pthw,] %>% 
  merge(ptwDf[grep(2, ptwDf$Cat), ], by.x = "Name", by.y = "pthw") %>% 
  select(c("Name", "NDE", "tA", "pGFdr", "Cat", "Status"))

spia2$Cat <- 2
spia24 <- rbind(spia2, spia4)

spia24[spia24$tA == 0, ]$tA <- -0.02
spia24[spia24$tA > 0 &spia24$tA < 0.02, ]$tA <- 0.08
# ptwDf[!ptwDf$pthw %in% spia24$Name, ]$pthw

# old ----
p1 <- ggplot(spia24[spia24$Cat == 2,], 
       aes(x = tA,  
           y = reorder(Name, ifelse(tA > 0, tA, -tA)), 
           fill=ifelse(tA < 0, "UP in ILC", "UP in IDC"),
           # alpha = -log10(pGFdr)
           )) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Name, 
                x = ifelse(tA > 0, -0.05, 0.05), 
                hjust = ifelse(tA > 0, 1, 0)), 
            size = 4, 
            alpha = 1) +
  xlab("Relative enrichment (tA) of pathways in mono/macro 2") +
  facet_grid(rows = vars(Status), scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(position="bottom")+
  scale_fill_discrete(name = "Pathway status") +
  # theme_minimal() +
  scale_x_continuous(limits=c(-30, 30), breaks=c(seq(-30,30,5), 0))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_blank())

p2 <- ggplot(spia24[spia24$Cat == 4,], 
             aes(x = tA,  
                 y = reorder(Name, ifelse(tA > 0, tA, -tA)), 
                 fill=ifelse(tA < 0, "UP in ILC", "UP in IDC"),
                 # alpha = -log10(pGFdr)
             )) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Name, 
                x = ifelse(tA > 0, -0.1, 0.1), 
                hjust = ifelse(tA > 0, 1, 0)), 
            size = 4, 
            alpha = 1) +
  xlab("Relative enrichment (tA) of pathways in mono/macro 4") +
  facet_grid(rows = vars(Status), scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(position="bottom")+
  scale_fill_discrete(name = "Pathway status") +
  # theme_minimal() +
  scale_x_continuous(limits=c(-45, 45), breaks=c(seq(-45,45,10), 0))+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_blank())

# new ----
options(scipen=10000)
df1 <- unique(spia24[spia24$Cat == 2,])
p1 <- ggplot(df1, 
             aes(x = tA,  
                 y = reorder(Name, ifelse(tA > 0, tA, -tA)), 
                 col=ifelse(tA < 0, "UP in ILC", "UP in IDC"),
                 # size = -log2(pGFdr)
                 # alpha = -log10(pGFdr)
             )) +
  geom_point(aes(size = -log2(pGFdr))) +
  geom_text(aes(label = Name), 
            col = "black", 
            size = 4.5,
            fontface = "bold",
            hjust = ifelse(df1$tA < 0, 0, 1),
            nudge_x = ifelse(df1$tA < 0, 0.5, -0.5)) +
  facet_grid(rows = vars(Status), scales = "free_y", 
             space = "free_y", switch = "y")+
  guides(size = guide_legend(title="FDR"),
         color = guide_legend(title="Pathway status")) +
  scale_x_continuous(limits=c(-30, 30), breaks=c(seq(-30,30,5), 0),
                     position="bottom") +
  scale_size_continuous(breaks = as.numeric(-log2(summary(df1$pGFdr)[c(2,4,6)])), 
                        labels = format(summary(df1$pGFdr)[c(2,4,6)], 
                                        digits = 3,
                                        scientific = T),
                        range = c(3, 4)) +
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_blank())

df2 <- unique(spia24[spia24$Cat == 4,])
p2 <- ggplot(df2, 
             aes(x = tA,  
                 y = reorder(Name, ifelse(tA > 0, tA, -tA)), 
                 col=ifelse(tA < 0, "UP in ILC", "UP in IDC"),
                 # alpha = -log10(pGFdr)
             )) +
  geom_point(aes(size = -log2(pGFdr))) +
  geom_text(aes(label = Name), 
            col = "black", 
            size = 4.5,
            fontface = "bold",
            hjust = ifelse(df2$tA < 0, 0, 1),
            nudge_x = ifelse(df2$tA < 0, 0.8, -0.8)) +
  facet_grid(rows = vars(Status), scales = "free_y", 
             space = "free_y", switch = "y")+
  guides(size = guide_legend(title="FDR"),
         color = guide_legend(title="Pathway status")) + 
  scale_x_continuous(limits=c(-45, 45), breaks=c(seq(-45,45,10), 0),
                     position="bottom") +
  scale_size_continuous(breaks = as.numeric(-log2(summary(df2$pGFdr)[c(2,4,6)])), 
                        labels = format(summary(df2$pGFdr)[c(2,4,6)], 
                                        digits = 3,
                                        scientific = T),
                        range = c(1, 6)) +
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        legend.key = element_rect(colour = NA, fill = NA),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        panel.background = element_blank())

png("../3.Result/spiaVisulBar.png", width = 15, height = 15, units = "in", res = 300)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom", widths = c(9, 5))
dev.off()

# newly new ----
df1 <- unique(spia24[spia24$Cat == 2,])
df2 <- unique(spia24[spia24$Cat == 4,])
df <- unique(df1)

df <- rbind(df %>% filter(tA > 0) %>% arrange(desc(tA)),
            df %>% filter(tA < 0) %>% arrange(tA))

df$pos1 <- 1:nrow(df)

df$tA <- ifelse(abs(df$tA) > 0.1, df$tA, sign(df$tA)*0.1)

df %>%
  arrange(tA) %>%
  filter(tA > 0) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=tA) ) +
  geom_bar(stat="identity", fill = "#cb4154") +
  geom_text(hjust = 1, size = 10, fontface = "bold", 
            aes( y = 0, label = paste(Name," "))) +
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
  ylim(0,ceiling(max(abs(range(df$tA)))) + 3) +
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 30),
        panel.grid = element_line(colour = "black"))

wd <- df %>%
  arrange(tA) %>%
  filter(tA > 0) %>% 
  nrow()

ggsave("circleBarClust2IDC.png", width = wd + 5, height = wd, units = "in", bg = "white")

df %>%
  arrange(desc(tA)) %>%
  filter(tA < 0) %>%
  mutate(Name=factor(Name, Name)) %>%
  ggplot( aes(x=Name, y=-tA) ) +
  geom_bar(stat="identity", fill = "skyblue") +
  geom_text(hjust = 0, size = 10, fontface = "bold", 
            aes( y = 0, label = paste(" ", Name))) +
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
  ylim(ceiling(max(abs(range(df$tA)))) + 3, 0) +
  theme_minimal()+
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 30),
        panel.grid = element_line(colour = "black"))

sd <- df %>%
  arrange(tA) %>%
  filter(tA < 0) %>% 
  nrow()

ggsave("circleBarClust2ILC.png", width = wd + 5, height = wd, units = "in", bg = "white")