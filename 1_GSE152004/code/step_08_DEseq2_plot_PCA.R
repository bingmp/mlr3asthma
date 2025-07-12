if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("ggplot2")) BiocManager::install("ggplot2")
if (!require("FactoMineR")) install.packages("FactoMineR")
if (!require("factoextra")) install.packages("factoextra")
if (!require("DESeq2")) install.packages("DESeq2")

this.path::this.dir() %>%
  dirname() %>%
  setwd()

rm(list = ls())

exp <- readRDS('RDS/vsdmat.RDS')
group <- readRDS('RDS/group.RDS')

# 7、PCA
pca <- as.data.frame(t(exp))
pca$group <- group$cluster

df.pca <- FactoMineR::PCA(pca[, -ncol(pca)], graph = F, scale = F)
fviz_pca_ind(df.pca,
                  geom = "point",
                  mean.point = T,
                  repel = F,
                  # palette   = c("#9ed048", "#1685a9"),
                  col.ind = pca[, "group"]
) +
  scale_color_discrete("") +
  scale_shape_discrete("")

ggsave(filename = 'result/PCA.png',height = 10,width = 12,units = 'cm',dpi = 300)
