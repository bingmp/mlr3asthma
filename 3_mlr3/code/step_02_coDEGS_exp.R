
if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("limma")) BiocManager::install("limma")

this.path::this.dir() %>%
  dirname() %>%
  setwd()

rm(list = ls())

GSE152004_group <- readRDS("../1_GSE152004/RDS/group_cluster.RDS")
GSE152004_group <- subset(GSE152004_group, cluster != "Th2_low")
GSE152004_vsdmat <- readRDS("../1_GSE152004/RDS/vsdmat.RDS")[,GSE152004_group$title]

GSE67472_group <- readRDS("../2_GSE67472/RDS/group.RDS")
GSE67472_exp <- readRDS("../2_GSE67472/RDS/norm_exp.RDS")[,GSE67472_group$geo_accession]

ids <- intersect(rownames(GSE67472_exp),rownames(GSE152004_vsdmat ))

all_exp <- cbind(GSE67472_exp[ids,],GSE152004_vsdmat[ids,])

norm_exp <- normalizeBetweenArrays(all_exp )  # limma: get normlization data

# boxlplot of norm_exp
if(F){
  if (!require("reshape2")) BiocManager::install("reshape2")
  
  exp_L <- melt( norm_exp )     # 
  colnames(exp_L)=c('probe','sample','value') # set colnames
  sample <- c(rep("GSE67472",nrow(GSE67472_group)),rep("GSE152004",nrow(GSE152004_group)))
  exp_L$group = rep(sample, each = nrow(norm_exp))  # set group
  
  if (!require("ggplot2")) BiocManager::install("ggplot2")
  p <- ggplot(exp_L,aes(x=sample,y=value,fill=group,color=group))+  
    geom_boxplot(alpha=0.5)+
    ggtitle("Boxplot of normlization data")+ 
    # ylim(0,10)+
    theme(
      plot.title = element_text(hjust = 0.5,size = 20),
      axis.title = element_text(size = 15),
      axis.text.x=element_blank(),

      # axis.text.x=element_text(angle=90,hjust = 1,size=7)
      legend.title = element_text(size = 15)
    ) 
  p
  ggsave(p, height = 6, width = 12,file='result/two_dataset_normlization_boxplot.PDF')
}

co_DEGs <- readRDS(file = 'RDS/coDEGs.RDS')
co_DEGs <- co_DEGs[order(co_DEGs)]

# get machlearn datasets
two_datasets <- as.data.frame( t(norm_exp[co_DEGs,]) )
two_datasets$group <-c(GSE67472_group$condition, GSE152004_group$condition)
two_datasets$sample <-  c(rep("GSE67472",nrow(GSE67472_group)),rep("GSE152004",nrow(GSE152004_group))) %>% as.factor()

two_datasets <- dplyr::select(two_datasets,sample,group,everything())
two_datasets

saveRDS(two_datasets,'RDS/two_dataset_co_DGEs_exp.RDS')
