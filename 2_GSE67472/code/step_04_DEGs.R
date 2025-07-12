
if(!require(limma)) install.packages("limma")
if(!require(tidyverse)) install.packages("tidyverse")

this.path::this.dir() %>% dirname() %>% setwd()

# 2、limma
rm(list = ls())

exp <- readRDS('RDS/norm_exp.RDS')
group <- readRDS('RDS/group.RDS')
condition <- group$cluster

# make contrast matrix & data.frame
design = model.matrix(~0+factor(condition)) 
colnames(design) = levels(factor(condition))
rownames(design) = colnames( condition )

dematrix <- as.data.frame(design)
head(dematrix) ;  table(dematrix )

# makeContrasts(Test-Control,levels = design) ---------------------------------------------------------------
contrast.matrix <- makeContrasts(Th2_high-Control,levels = design)
contrast.matrix

# 第一步lmFit，# lmFit为每个基因给定一系列的阵列来拟合线性模型
fit <- lmFit(exp,design)
# 第二步eBayes，# eBayes给出了一个微阵列线性模型拟合，通过经验贝叶斯调整标准误差到一个共同的值来计算修正后的t统计量、修正后的f统计量和微分表达式的对数概率。
fit1 <- contrasts.fit(fit, contrast.matrix)
fit1 <- eBayes(fit1)

# 3、diff.gene
temoutput <- topTable(fit1, coef=1, adjust="BH", n=Inf) #
temoutput <- na.omit(temoutput) # remove NA

# get diff.genes ---------------------------------------------------------------
temoutput$change = ifelse(temoutput$P.Value>0.05,'stable', # stable gene
                          ifelse( temoutput$logFC> log2(1.5),'up', # up gene
                                  ifelse( temoutput$logFC < -log2(1.5),'down','stable') ) # down or stable gene
)
table(temoutput$change)
temoutput <- temoutput[order(temoutput$logFC ), ]
temoutput$ID <- rownames(temoutput)

tem <- temoutput[!temoutput$change=='stable',]
table(tem$change)

diff.m <- exp[rownames(tem),]
diff.m <- as.data.frame(t(diff.m  ))
diff.m$group <- group$condition

save(tem,temoutput, diff.m, group,file = 'RDS/th2_h_vs_con_diff_all.Rdata')
saveRDS(tem, file = 'RDS/GSE67472_DEGs.RDS')

