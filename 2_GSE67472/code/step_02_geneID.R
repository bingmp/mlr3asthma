
if(!require(org.Hs.eg.db)) BiocManager::install("org.Hs.eg.db")
if(!require(tidyverse)) install.packages("tidyverse")

this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())


# GPL: probe to symbol gene names
raw_exp <- readRDS('RDS/raw_exp.RDS')

gpl <- getGEO("GPL16311",destdir = "rawdata")
colnames(Table(gpl))

Table(gpl)[1:20,c(1,2)] 
## get Symbol names
probe2gene <-   Table(gpl)[,c(1,2)] 

colnames(probe2gene) <- c("ID","Entrez")

g_symbol <- mapIds(x = org.Hs.eg.db,#注释包
                   keys = probe2gene$Entrez, #需要转换的基因ID
                   keytype = "ENTREZID", #需要转换的类型
                   column = "SYMBOL")  #需要转换为的类型

probe2gene$Symbol <- as.character(g_symbol )
probe2gene$name <- names(g_symbol)
which(probe2gene$name!=probe2gene$Entrez)

probe2gene <- na.omit( probe2gene[,c(1,3)] )

exp <- as.data.frame( raw_exp )
exp <- exp %>% 
  mutate(ID=rownames(exp)) %>% 
  inner_join(probe2gene, by="ID") %>%
  dplyr::select(ID, Symbol, everything())

exp <- exp[!duplicated(exp$Symbol),]  # duplicate gene symbol
rownames(exp) <- exp$Symbol

exp <- exp[,-(1:2)] # remove column of probe & symbol 

exp <- as.matrix(exp)

saveRDS( probe2gene,file = 'RDS/probe2gene.RDS')
saveRDS( exp, file = 'RDS/exp_id.RDS')

