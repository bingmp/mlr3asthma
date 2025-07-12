if(!require(ggplot2)) install.packages("ggplot2")
if(!require(pROC)) install.packages("pROC")
if(!require(tidyverse)) install.packages("tidyverse")

this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())

load("RDS/th2_h_vs_con_diff_all.Rdata")

co_DEGs <-c("SCGB1A1","CSTA","CYP24A1","CCL26", "TFF3","KCNJ16","C3","CD44") # 凌耀政

exp <- readRDS('RDS/norm_exp.RDS')
group <- readRDS('RDS/group.RDS')

roc_func <- function(data, group, name_gene){
  roc.m <- data.frame(data= data[name_gene, ],
                      group=group$condition
  )
  colnames(roc.m) <- c(name_gene,'group')
  roc.fit <- pROC::roc(group~.,data=roc.m,aur=TRUE,
                       ci=TRUE, # show 95%CI
                       # percent=TRUE, 
                       smooth=F,
                       levels=c('Control','Asthma')
  )
  
  if(!require(ggplot2)) install.packages("ggplot2")
  p <- pROC::ggroc(roc.fit, 
                   legacy.axes = T, # set x.axis: 1- spe
                   color='pink') + 
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                 color="darkgrey", linetype=4) +
    annotate("text", x=0.75, y=0.36, size=4,
             label=paste(roc.fit[["predictor.name"]],"- AUC = ", round(roc.fit$auc,3) ) )+
    theme_bw()+
    theme(
      axis.title = element_text(size=15),
      axis.text = element_text(size=10),
      legend.text = element_text(size=8),
      legend.title = element_blank()
    )+
    scale_color_discrete()
  return(p)
}


for (i in co_DEGs) {
  ggsave( roc_func(exp, group, i), height = 6, width = 6 ,filename = paste0('result/ROC/',i,'.PDF'))
}


