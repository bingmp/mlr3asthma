if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("ggplot2")) BiocManager::install("ggplot2")
if (!require("pROC")) install.packages("pROC")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("DESeq2")) install.packages("DESeq2")

this.path::this.dir() %>%
  dirname() %>%
  setwd()

rm(list = ls())

exp <- readRDS('RDS/vsdmat.RDS')
group <- readRDS('RDS/group_cluster.RDS')
group <- subset(group, cluster!="Th2_low")

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

# 8、ROC
roc_func <- function(data, group, name_gene) {
  roc.m <- data.frame(
    data = data[name_gene, ],
    group = group$cluster
  )
  colnames(roc.m) <- c(name_gene, "group")
  roc.fit <- pROC::roc(group ~ .,
                       data = roc.m, aur = TRUE,
                       ci = TRUE, # 显示95%CI
                       # percent=TRUE, ##是否需要以百分比显示
                       smooth = F,
                       levels = c("Control", "Th2_high")
  )
  
  if (!require(ggplot2)) install.packages("ggplot2")
  p <- pROC::ggroc(roc.fit, legacy.axes = T, color = "red") + # 将X轴改为0-1，（默认是1-0）
    geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
                 color = "darkgrey", linetype = 4
    ) +
    annotate("text",
             x = 0.75, y = 0.36, size = 4,
             label = paste(roc.fit[["predictor.name"]], "- AUC = ", round(roc.fit$auc, 3))
    ) +
    theme_bw() + # 设置背景
    theme(
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),
      legend.text = element_text(size = 6),
      legend.title = element_blank()
    ) +
    scale_color_discrete()
  return(p)
}

roc_func(data = exp,group = group, "CD44")

for (i in co_DEGs) {
  ggsave(roc_func(data = exp,group = group, i),
         height = 10, width = 10, units = "cm", dpi = 300,
         filename = paste0("result/ROC/511_control_Th2_high_roc_", i, ".png")
  )
}
