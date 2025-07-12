
if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("ggplot2")) BiocManager::install("ggplot2")

this.path::this.dir() %>%
  dirname() %>%
  setwd()

rm(list = ls())

load("RDS/511_control_Th2_high_diff.Rdata")
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
# 6、火山图 -------------------------------------------------------------------
logFC_cutoff <- log2(1.5)

df <- na.omit(temoutput)
# 设置火山图的标题
this_tile <- paste(
  "Cutoff for logFC is ", round(logFC_cutoff, 3),
  "\nThe number of up gene is ", nrow(tem[tem$change == "up", ]),
  "\nThe number of down gene is ", nrow(tem[tem$change == "down", ])
)

temoutput$label <- NA
for (i in co_DEGs) {
  temoutput[which(temoutput$ID == i), "label"] <- i
}

# 画火山图
ggplot(
  data = temoutput,
  aes(
    x = log2FoldChange,
    y = -log10(pvalue), # 这里将 pvalue 取负对数
    color = change
  )
) +
  geom_point(alpha = 0.4, size = 1) + #  绘制点图
  ggrepel::geom_text_repel(
    data = temoutput,
    aes(
      x = log2FoldChange,
      y = -log10(pvalue),
      label = label
    ),
    max.overlaps = 30,
    size = 3,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"),
    segment.color = "black",
    show.legend = FALSE
  ) +
  xlab("log2 fold change") +
  ylab("-log10 pvalue") + # 轴标签
  # ggtitle(this_tile) +
  theme_bw(base_size = 10) +
  theme(
    axis.text  = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 15, hjust = 0.5),
    plot.title = element_text(size = 8, hjust = 0.5)
  ) +
  scale_color_manual(values = c("blue", "black", "red")) # 设定颜色

ggsave(height = 6, width = 8, filename = "result/511_control_Th2_high_volcano.PDF") # 输出

  