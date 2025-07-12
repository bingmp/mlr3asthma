
if(!require(WGCNA)) BiocManager::install("WGCNA")
if (!require("magrittr")) BiocManager::install("magrittr")

rm(list = ls())

vstMat <- readRDS("RDS/vstMat.RDS")
group <- readRDS("RDS/group.RDS")

load('RDS/wgcna_gene2module.RData')


# rownames(vstMat)
GSE152004 <- vstMat[subset(gene2module, module== "brown4")[,1],
                    which(group$condition=="Asthma")]

GSE152004 <- t(GSE152004)
dim(GSE152004); head(GSE152004[,1:6])

# 利用碎石图确定聚类个数
wss <- (nrow(GSE152004)-1)*sum(apply(GSE152004,2,var)) # 计算离均差平方和
for (i in 2:15) wss[i] <- sum(kmeans(GSE152004, 
                                     centers=i)$withinss) #计算不同聚类个数的组内平方和
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares") # 绘图

# Ward层次聚类
# d <- dist(GSE152004, method = "euclidean") # 计算各个样本点之间的欧氏距离
# fit2 <- hclust(d, method="ward.D") #进行Ward层次聚类

fit2 <- hclust(dist(GSE152004, method = "euclidean"), method="ward.D")
groups <- cutree(fit2, k=2) # 设定聚类个数为2
table(groups)

pdf("IMAGES/441.Cluster_Dendrogram.hclust.pdf")
plot(fit2) # 绘制树状图展示聚类结果
# 给聚成的2个类别加上红色边框
rect.hclust(fit2, k=2, border="red")
dev.off()

library(ggtree)
group_h <- data.frame(label=names(groups ), cluster = factor(groups ))  # 创建hclust的分组信息用于ggtree绘图
group_h$cluster <- factor(group_h$cluster,
                          levels = c(1,2),
                          labels = c("Th2-high: 257","Th2-low: 184"))
p <- ggtree(ape::as.phylo( fit2 ), 
            linetype='dashed', 
            color = "#487AA1", 
            layout = "fan") %<+% group_h +
  geom_tiplab(aes(color = cluster), size=4) + 
  theme(legend.position = "right",
        legend.title = element_text(size = 18),
        legend.text =  element_text(size = 16))

ggsave(p,file='IMAGES/441.hclust_cic.pdf',width = 20*1.414,height = 20)

group$cluster <- "Control"
group[which(group$condition=="Asthma"),"cluster"] <- paste0("cluster", groups)
group$cluster <- factor(group$cluster,
                        levels = c("Control","cluster1", "cluster2" ),
                        labels = c("Control","Th2_high","Th2_low"))

table(group$cluster)
saveRDS(group,file = "RDS/group_cluster.RDS")
