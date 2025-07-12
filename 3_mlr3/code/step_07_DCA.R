
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("rmda")) install.packages("rmda")
if (!require("ggDCA")) remotes::install_github("yikeshu0611/ggDCA")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("rms")) install.packages("rms")
if (!require("ranger")) install.packages("ranger")


rm(list = ls())
# 读取数据
df1 <- readRDS("RDS/two_dataset_co_DGEs_exp.RDSS")[,-1]
df1$group <- as.numeric(df1$group) - 1


if (F) {
  ddist <- datadist(df1)
  options(datadist = "ddist")
  
  `DE-CMRG` <- lrm(group ~ ., df1)
  CSTA <- lrm(group ~ CSTA, df1)
  C3 <- lrm(group ~ C3, df1)
  CD44 <- lrm(group ~ CD44, df1)
  
  dca(`DE-CMRG`,CSTA,C3,CD44, model.names = c( "DE-CMRG", "CSTA", "C3","CD44") ) %>%
    ggplot()+ggtitle("Decision curve analysis of Logistic regression")+
    theme(
      plot.title =  element_text(size = 15,hjust = 0.5,vjust = 0.5)
    )+
    ggsci::scale_color_lancet()
  
  # ggsave(filename = "result/step_07_dca.pdf", width = 16, height = 8, dpi = 300, units = "cm")
  

  CD44 <- lrm(group ~ CD44, df1)
  CCL26 <- lrm(group ~ CCL26, df1)
  CSTA <- lrm(group ~ CSTA, df1)
  KCNJ16 <- lrm(group ~ KCNJ16, df1)
  TFF3 <- lrm(group ~ TFF3, df1)
  
  dca(CD44,CCL26,CSTA,KCNJ16,TFF3, model.names = c( "CD44","CCL26", "CSTA", "KCNJ16","TFF3") ) %>%
    ggplot()+ggtitle("Decision curve analysis of Logistic regression")+
    theme(
      plot.title =  element_text(size = 15,hjust = 0.5,vjust = 0.5)
    )+
    ggsci::scale_color_lancet()
  
  ggsave(filename = "result/step_07_dca_5gene.pdf", width = 16, height = 8, dpi = 300, units = "cm")
  
}

