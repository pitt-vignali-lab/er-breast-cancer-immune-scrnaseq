library(Seurat)
library(scales)
library(tidyverse)

dat <- readRDS("../TIL_mono_macro.RDS")

# dat <- subset(dat, samples %in% c("IDC", "ILC")) #IDC, ILC, PBL.IDC, PBL.ILC
# 
# dat <- FindNeighbors(dat)
# 
# dat <- FindClusters(dat, resolution = 0.4)
# 
# dat <- subset(dat, samples == "ILC")

plt <- DimPlot(dat, label = T, label.size = 8) + 
  coord_fixed(ratio = 1) +
  NoLegend()

plt

# x-range
round(diff(layer_scales(plt)$x$range$range)/2, 1)
# y-range
round(diff(layer_scales(plt)$y$range$range)/2, 1)

dat$cluster <- Idents(dat)


# data infor ----
cells <- data.frame(id = Cells(dat),
                    type = dat$samples,
                    types = dat$samples,
                    Embeddings(dat, reduction = "umap"),
                    dat@meta.data["cluster"])


cells$id <- gsub("m107", "BC1", 
            gsub("m108", "BC2",
            gsub("m194", "BC3",
            gsub("m226", "BC4",
            gsub("m234", "BC5",
            gsub("m268", "BC6",
            gsub("m272", "BC7",
            gsub("m302", "BC8",
            gsub("m349", "BC9",
            gsub("m365", "BC10",
            gsub("m533", "BC11",
            gsub("m571", "BC12",
            gsub("m615", "BC13",
            gsub("m888", "BC14",
            gsub("m899", "BC15",
                      cells$id
                 )))))))))))))))

cells$type <- ifelse(cells$type %in% c("IDC", "ILC"),
                     cells$type <- "TIL",
                     cells$type <- "PBMC")

cells$newId <- paste(ifelse(sapply(strsplit(cells$id, '_\\s*'), `[`, 1) %in% 
                               paste0("BC", 1:10),
                             paste(sapply(strsplit(cells$id, '_\\s*'), `[`, 1), "PBMC_TIL", sep = "_"), 
                             paste(sapply(strsplit(cells$id, '_\\s*'), `[`, 1), cells$type, sep = "_")),
                      sapply(strsplit(cells$id, '_\\s*'), `[`, 2),
                     sep = "_")


cells <- cells %>% remove_rownames() %>% 
  column_to_rownames(var = "newId") %>% 
  select(-c("id", "type"))

write.csv(rownames(cells), file = "cellID.csv", row.names = FALSE)
write.csv(cells[2:3], file = "umap.csv")
write.csv(cells[4], file = "clusters.csv")
write.csv(cells[1], file = "types.csv")

cols <- RColorBrewer::brewer.pal(n=length(unique(table(Idents(dat)))),"Set2")
cols

cols <- hue_pal()(length(table(Idents(dat))))
show_col(hue_pal()(length(table(Idents(dat)))))
names(cols) <- levels(dat$cluster)
paste(cols, sep = ", "  )

cols <- RColorBrewer::brewer.pal(n=length(unique(table(Idents(dat)))),"Set2")

# x-range
round(diff(layer_scales(plt)$x$range$range), 1)

# y-range
round(diff(layer_scales(plt)$y$range$range), 1)