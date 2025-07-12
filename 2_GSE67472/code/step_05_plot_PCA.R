if(!require(FactoMineR)) install.packages("FactoMineR")
if(!require(factoextra)) install.packages("factoextra")
if(!require(tidyverse)) install.packages("tidyverse")

this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())

exp <- readRDS('RDS/norm_exp.RDS')
group <- readRDS('RDS/group.RDS')

pca <- as.data.frame(t(exp ))
pca$group <-  group$condition

df.pca <- FactoMineR::PCA( pca[,-ncol(pca)], graph = F,scale=F) 

p <- fviz_pca_ind(df.pca,
                  geom = "point",
                  mean.point  = T,
                  repel = F,
                  # palette   = c("#9ed048", "#1685a9"),
                  col.ind = pca[,"group"]
)+
  scale_color_discrete("")+
  scale_shape_discrete("")

ggsave(p, filename = 'result/PCA.png',height = 10,width = 12,units = 'cm',dpi = 300)

