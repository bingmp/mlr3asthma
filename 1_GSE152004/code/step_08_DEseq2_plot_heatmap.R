if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("ggplot2")) BiocManager::install("ggplot2")
if (!require("pheatmap")) BiocManager::install("pheatmap")

this.path::this.dir() %>%
  dirname() %>%
  setwd()

rm(list = ls())

exp <- readRDS('RDS/vsdmat.RDS')
group <- readRDS('RDS/group.RDS')

co_DEGs <- c(
  "SCGB1A1",
  "CYP24A1",
  "CCL26",
  "TFF3",
  "KCNJ16",
  "CSTA",
  "C3",
  "CD44"
)

heatmap_func <- function(data, group){
  
  choose_matrix1 <- data[co_DEGs, ]
  
  annotation_col = data.frame(
    group = group$condition
  )
  rownames(annotation_col) <- colnames(choose_matrix1)
  
  d <- choose_matrix1
  annotation_col$group <- factor(annotation_col$group)
  p <- pheatmap::pheatmap(d,  cluster_cols = F, scale = "row",
                          show_colnames = F,
                          cutree_rows = 2,
                          silent = F,
                          treeheight_col = 20,
                          treeheight_row = 20,
                          fontsize = 6,
                          annotation_col = annotation_col ) 
  return(p)
}

heatmap_func(data = exp, group = group)

ggsave(filename = 'result/GSE152004 asthma co-DEG heatmap.png',
       height = 8, width = 20,units = 'cm',dpi = 300)
