edges <- ctRes[c("from", "to", "type")]
vertices <- data.frame(name = unique(c(ctRes$from, ctRes$to)))

vertices <- vertices |> 
  separate(name, into = c("tissue", "gene"), sep = "_", remove = F) |>
  mutate(tissue =tissue[order(tissue)])

mygraph <- graph_from_data_frame(edges, vertices = vertices )


ggraph(mygraph, layout = "linear", circular=TRUE) +
  geom_edge_arc(aes(col = type)) +
  geom_node_point(aes(x = x * 1.05, y = y* 1.05, col = tissue)) +
  geom_node_text(aes(label = gene, x = x * 1.1, y = y* 1.1,
                     angle = ifelse(atan(-(x/y))*(180/pi) < 0,
                                    90 + atan(-(x/y))*(180/pi),
                                    270 + atan(-x/y)*(180/pi)),
                     hjust = ifelse(x > 0, 0 ,1)), size = 3.5)+
  theme(
    legend.position="bottom",
    plot.margin=unit(c(0,0,0,0),"cm"),
  ) +
  theme_graph() +
  expand_limits(x = c(-1.2, 1.2), y = c(-1.2, 1.2)) + coord_fixed()

  

