
if(!require(limma)) install.packages("limma")
if(!require(reshape2)) install.packages("reshape2")
if(!require(magrittr)) install.packages("magrittr")
if(!require(ggplot2)) install.packages("ggplot2")

this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())

norm_exp <- readRDS('RDS/exp_id.RDS') %>% normalizeBetweenArrays( ) # limma: get normlization data

saveRDS(norm_exp, file = 'RDS/norm_exp.RDS')
saveRDS(norm_exp, file = 'RDS/GSE67472_matrix.RDS')

# boxlplot of norm_exp
exp_L <- melt( norm_exp )     # 
colnames(exp_L)=c('probe','sample','value') # set colnames

group <- readRDS('RDS/group.RDS')
exp_L$group = rep(group$condition, each = nrow(norm_exp))  # set group

p <- ggplot(exp_L,aes(x=sample,y=value,fill=group))+  
  geom_boxplot()+
  ggtitle("Boxplot of normlization data")+ 
  ylim(0,10)+
  
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x=element_text(angle=90,hjust = 1,size=7)
  ) 

p

ggsave(p, height = 10, width = 20,file='result/normlization_boxplot.PDF')

