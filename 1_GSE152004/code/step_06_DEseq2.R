
if (!require("ggplot2")) BiocManager::install("ggplot2")
if (!require("DESeq2")) install.packages("DESeq2")
if (!require("BiocParallel")) install.packages("BiocParallel")
if (!require("magrittr")) BiocManager::install("magrittr")

this.path::this.dir() %>%
  dirname() %>%
  setwd()
rm(list = ls())

# 1、DEseq2

group <- readRDS("RDS/group_cluster.RDS")
raw_counts <- read.table("rawdata/GSE152004_695_raw_counts.txt.gz", header = T, sep = "\t", row.names = 1)
raw_counts <- as.matrix(raw_counts)

exp <- raw_counts
good_genes <- rowSums(raw_counts >= 10) >= (ncol(raw_counts) * .15)
group <- subset(group, cluster != "Th2_low")

exp <- raw_counts[good_genes, group$title]
colnames(exp) == group$title
condition <- group$cluster
head(exp[, 1:5])

# 2.2 构建分组数据框 colData
colData <- as.data.frame(group$cluster)
rownames(colData) <- colnames(exp) # 设置行名为样本名信息
colnames(colData) <- "condition" # 列名
colData$condition <- factor(colData$condition)
head(colData)

# 3.1 构建DEseq2对象
dds <- DESeqDataSetFromMatrix(
  countData = exp,
  colData = colData,
  design = ~condition
)
dds

# 3.2 DEseq2
dds <- DESeq(dds, parallel = T)

save(group, dds, file = "RDS/511_control_Th2_high_dds.Rdata")

# 标准化
# 方差稳定变换，The variance stabilizing transformation(vst)
vsd <- vst(dds)
class(vsd)
vsdmat <- assay(vsd) # 提取转化后的矩阵
plotPCA(vsd, intgroup = c("condition"))
saveRDS(vsdmat, file = "RDS/vsdmat.RDS")

