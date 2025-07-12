
if (!require("magrittr")) BiocManager::install("magrittr")

this.path::this.dir() %>%
  dirname() %>%
  setwd()

# 4、提取DEseq2分析的结果 ---------------------------------------------------------
rm(list = ls())
load("RDS/511_control_Th2_high_dds.Rdata")

table(group$cluster)
# 5.1 set contrast. contrast = c('condition','test','control')
res <- results(dds, contrast = c("condition", "Th2_high", "Control"))

temoutput <- as.data.frame(res)
temoutput <- temoutput[order(temoutput$log2FoldChange), ]

temoutput$change <- as.factor(
  ifelse(
    temoutput$pvalue < 0.05 & abs(temoutput$log2FoldChange) > log2(1.5),
    ifelse(
      temoutput$log2FoldChange > log2(1.5), "up", "down"
    ),
    "stable"
  )
)

table(temoutput$change) # 查看基因上、下调情况

temoutput$ID <- rownames(temoutput)
temoutput <- temoutput[order(temoutput$log2FoldChange), ]
tem <- subset(temoutput, change != "stable")

temoutput$ID <- rownames(temoutput)
temoutput <- temoutput[order(temoutput$log2FoldChange), ]
tem <- subset(temoutput, change != "stable")

# write.csv(temoutput, row.names = F, file = "result/511_control_Th2_high_all_diff.csv")
save(temoutput, tem, file = "RDS/511_control_Th2_high_diff.Rdata")

saveRDS(tem, file = "RDS/GSE152004_DEGs.RDS")
 