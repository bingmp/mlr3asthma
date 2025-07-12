
if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("ComplexUpset")) BiocManager::install("ComplexUpset")
if (!require("ggplot2")) BiocManager::install("ggplot2")

this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())

GSE152004_DEGs <- readRDS("../1_GSE152004/RDS/GSE152004_DEGs.RDS")
GSE67472_DEGs <- readRDS("../2_GSE67472/RDS/GSE67472_DEGs.RDS")
DEGs_Copper <- read.csv('rawdata/2067-GeneCards-SearchResults.csv')

coDEGs <- intersect(GSE152004_DEGs$ID, GSE67472_DEGs$ID) %>% intersect(DEGs_Copper$Gene.Symbol)
saveRDS(coDEGs,file = "RDS/coDEGs.RDS")

UpSet <- UpSetR::fromList(list(GSE67472=GSE67472_DEGs$ID, 
                               GSE152004=GSE152004_DEGs$ID,
                               Copper=DEGs_Copper$Gene.Symbol))

head(UpSet)

p <- ComplexUpset::upset(UpSet, colnames(UpSet), 
                         name='group', min_degree = 2,
                         width_ratio=0.5, wrap = T,
                         queries=list(
                           upset_query(
                             intersect=c('GSE67472', 'GSE152004',"Copper"),
                             color='orange', fill='orange'
                           )
                         ) ) +
  ggtitle("Upset plot of co-DEGs")+
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size=15 ),
    plot.background = element_rect(fill = NA)
    
  )

ggsave(p, filename = 'result/upset.PDF',width = 16, height = 8, units = 'cm',dpi = 300)

