
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(pheatmap)) install.packages("pheatmap")

this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())

exp <- readRDS('RDS/norm_exp.RDS')
group <- readRDS('RDS/group.RDS')

co_DEGs <-c("SCGB1A1","CSTA","CYP24A1","CCL26", "TFF3","KCNJ16","C3","CD44") # 凌耀政

heatmap_func <- function(data, group, main=''){
  
  choose_matrix1 <- data[co_DEGs, group$geo_accession]
  
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
                          main= main,
                          treeheight_col = 20,
                          treeheight_row = 20,
                          fontsize = 6,
                          annotation_col = annotation_col ) 
  return(p)
}


p <- heatmap_func(exp, group, main = 'co-DEGs heatmap of GSE67472')

ggsave(p, filename = 'result//GSE67472_co-DEG_heatmap.PDF', width = 11, height = 8, units = 'cm',dpi = 300)
