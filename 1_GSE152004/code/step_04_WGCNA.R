
if(!require(WGCNA)) BiocManager::install("WGCNA")
if (!require("magrittr")) BiocManager::install("magrittr")

this.path::this.dir() %>% dirname() %>% setwd()

# WGCNA analysis
rm(list = ls())

asthma_expMat <- readRDS("RDS/asthma_expMat.RDS")
group <- readRDS("RDS/group.RDS")

# pick soft threshold = 9
softPower <- 9 # sft$powerEstimate + 1

####step by step WGCNA

#1. Create Similarity Matrix
pearson <- WGCNA::cor(as.matrix(asthma_expMat),method="pearson")

#2. Convert Similarity Matrix to Adjacency Matrix using Adjacency Function
adjacency.p <- adjacency.fromSimilarity(pearson,type = "signed",power=softPower)

#3. Convert Adjacency to TOM dissimilarity
# This process will spend a lot of time.
TOM <- TOMsimilarity(adjacency.p,TOMType = "signed",TOMDenom = "min")

saveRDS(TOM,file = "RDS/TOM.RDS")

TOMdissim.p <- 1 - TOM

#4. Perform hierarchical clustering on the dissimilarity matrix
geneTree = hclust(as.dist(TOMdissim.p), method = "average")

deepsplit <- 0.82

#5. cut tree based on hierarchical clustering
modules = cutreeDynamic(dendro = geneTree,method='hybrid', 
                        distM = TOMdissim.p,
                        pamStage = F, 
                        pamRespectsDendro = F,
                        maxCoreScatter = min(geneTree$height)+deepsplit*(max(geneTree$height)-min(geneTree$height)),
                        minGap = (1-deepsplit)*(max(geneTree$height)-min(geneTree$height)),
                        cutHeight = quantile(geneTree$height,.99), minClusterSize=30)

modcolors = labels2colors(modules)

table(modcolors)

# dendogram
pdf("result/441.res.dendrogram_signed_wgcna.pdf",width=14,height=7)
#plotDendroAndColors(dendro=geneTree,colors=cbind(modcolors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE,main=' ')
plotDendroAndColors(dendro=geneTree,colors=modcolors, dendroLabels = FALSE,main=' ')
dev.off()

#######
MEs <- moduleEigengenes(asthma_expMat, modcolors)$eigengenes
rownames(MEs) <- rownames(asthma_expMat)
MEs_no_grey <- MEs[,-which(colnames(MEs) %in% c("MEgrey"))]

write.table(MEs_no_grey,  file="TABLES/441.res.WGCNA.ME.txt", sep='\t', quote=F)

# create table for genes and their corresponding module
gene2module = data.frame(gene=colnames(asthma_expMat), module=modcolors)
gene2module<-na.omit(gene2module)
write.table(gene2module, file="TABLES/441.res.WGCNA.gene2module.txt",sep="\t",quote=F,row.names=F)

datKME <- signedKME(asthma_expMat, MEs)
write.table(datKME, file="TABLES/441.res.WGCNA.KME.txt",sep="\t",quote=F)

gene2module_with_cor <- gene2module
gene2module_with_cor$cor <- NA

for(i in unique(gene2module_with_cor$module)) {
  kME_name <- paste0("kME",i)
  idx <- which(gene2module_with_cor$module==i)
  gene.idx <- as.character(gene2module_with_cor[idx,"gene"])
  gene.idx <- gsub("-",".",gene.idx)
  gene2module_with_cor$cor[idx] <- datKME[gene.idx,kME_name]
  #print(kME_name)
}
write.table(gene2module_with_cor, file="TABLES/441.res.WGCNA.gene2module.with.cor.txt",sep="\t",quote=F,row.names=F)


### hub genes
ADJ1 <- abs(cor(asthma_expMat, use="p"))^softPower

if(F){
  ADJ1 <- abs(cor(mydata,use = 'p'))^softPower
  
  pdf("result/441.hist_k.pdf",width=14,height=7)
  k <- apply(ADJ1,2,sum) -1
  hist(k)
  dev.off()
  
  pdf("result/441.hist_k_lm.pdf",width=14,height=7)
  cut1=cut(k,10)
  binned.k=tapply(k,cut1,mean)
  freq1=tapply(k,cut1,length)/length(k)
  plot(log10(binned.k),log10(freq1+.000000001),xlab="log10(k)",ylab="log10(p(k))")
  xx= as.vector(log10(binned.k))
  lm1=lm(as.numeric(log10(freq1+.000000001))~ xx )
  lines(xx,predict(lm1),col=1)
  title(paste( "scale free R^2=",as.character(round(summary(lm1)$adj.r.squared,2)),", slope=", round(lm1$coefficients[[2]],2)))
  dev.off()
  
}

Alldegrees1 <- intramodularConnectivity(ADJ1, gene2module$module)
hub_genes <- data.frame()
for(i in 1:length(unique(gene2module$module))){
  Alldegrees1$Module = gene2module$module
  tmp = Alldegrees1[Alldegrees1$Module == unique(gene2module$module)[i], ]
  hub_genes<-rbind(hub_genes, head(tmp[order(tmp$kWithin, decreasing=T),], n=nrow(tmp)))
}

write.table(hub_genes, file="TABLES/441.res.WGCNA.hub_genes.txt",sep="\t",quote=F)

##### module dendogram
moduleTree <- hclust(dist(t(MEs_no_grey)), method = "average")

pdf("result/441.res.ModuleTree.wgcna.pdf")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(moduleTree, main = "Module clustering", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

save(gene2module,ADJ1,  gene2module_with_cor, hub_genes, 
     file = "RDS/wgcna_result.RData")

save(gene2module,ADJ1, file = "RDS/wgcna_gene2module.RData")

