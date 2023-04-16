require(tidyverse)
require(igraph)
require(ggraph)
require(pracma)

# rm(list = ls()[ls() != "top_stats"])
# 
# gc()
ctRes <- read.csv("../3.Result/2022-06-23_List_of _receptor_ligand_IDC.csv")

ctRes <- ctRes %>% 
  select(2:4) %>% 
  separate(interaction, into = c("fromG", "toG"), sep = "_", remove = F) %>%
  separate(interaction_pairs, into = c("fromT", "toT"), sep = "_", remove = F) %>% 
  rename('type' = 1 ) 

ctRes$from <- paste(ctRes$fromT, ctRes$fromG, sep = "_")
ctRes$to <- paste(ctRes$toT, ctRes$toG, sep = "_")
ctRes$allInteraction <- paste(ctRes$from, ctRes$to, sep = "_")

ctRes$type <- sapply(ctRes$type, function(x) paste(
  sort(strsplit(x, split = "_")[[1]]), collapse = "_"))


ctRes <- ctRes[ctRes$from != ctRes$to, ]

d1 <- data.frame(from="origin", to = unique(c(ctRes$fromT, ctRes$toT)))
d2 <- unique(rbind(ctRes[c("fromT", "from")] %>% rename(from = 1, to = 2),
                 ctRes[c("toT", "to")] %>% rename(from = 1, to = 2)))

edges <- rbind(d1, d2)

connect <- ctRes[ctRes$from != ctRes$to, ]


vertices  <-  data.frame(
  name = unique(c(as.character(edges$from), as.character(edges$to)))
) 

vertices$id <- NA
myleaves <- which(is.na( match(vertices$name, edges$from) ))

nleaves <- length(myleaves)
vertices$id[ myleaves ] <- seq(1:nleaves)


vertices <- vertices %>% separate(name, into = c("tissue", "gene"), sep = "_", 
                                  remove = F)


mygraph <- graph_from_data_frame(edges, vertices = vertices )

from = match(connect$from, vertices$name)
to = match(connect$to, vertices$name)

connDf <- data.frame(from = from, to = to,
                     fromG = vertices$name[from],
                     toG = vertices$name[to])

connDf$allInteraction = paste(connDf$fromG, connDf$toG, sep = "_")

connDf <- merge(ctRes[c("type", "allInteraction")], connDf, by = "allInteraction")

# figure ----
colors_use <- RColorBrewer::brewer.pal(n=5,"Set2")[c(1, 4, 5)] # ilc

colors_use <- RColorBrewer::brewer.pal(n=5,"Set2") # idc
p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) 

for (i in 1:length(unique(connDf$type))) {
  type = unique(connDf$type)[i]
  p = p + geom_conn_bundle(data = get_con(from = connDf[connDf$type == type,]$from, 
                                          to = connDf[connDf$type == type,]$to),
                           aes(x = x * .6, y = y * .6),
                           colour = colors_use[i], tension = 0.5,
                           width = 0.5)
}

p + geom_node_point(aes(x = .65 * x, y = .65 * y,
                        filter = leaf, col = tissue), 
                  size= 5, alpha = 0.6) +
  # geom_node_text(aes(filter = leaf, label=name)) +
  geom_node_text(aes(filter = leaf, label=gene, 
                     x = .7 * x, y = .7 * y,
                     angle = ifelse(atan(-(x/y))*(180/pi) < 0,
                                    90 + atan(-(x/y))*(180/pi),
                                    270 + atan(-x/y)*(180/pi)),
                     hjust = ifelse(x > 0, 0 ,1)), 
                 size=4, fontface = 2) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + coord_fixed()

ggsave("../3.Result/celltalkerIDC_1.pdf", width = 5, height = 5, units = "in",
       dpi = 300, bg = "white")




pdf("../3.Result/celltalkerILCLgd.pdf", width = 10, height = 8)

plot.new()
legend(x = 0.5, y = 0.5,
       legend = unique(connDf$type),
       col = colors_use, lwd = 5, box.col = "white")
legend(x = 0.54, y = 0.3,
       legend = c("CD4 T cells", "CD8 T cells", "Mono/Macs"),
       col =  scales::hue_pal()(3),
       pt.cex = 2,
       pch = 16, 
       box.col = 'white')
dev.off()
