if(!require(DESeq2)) BiocManager::install("DESeq2")
if (!require("magrittr")) BiocManager::install("magrittr")

this.path::this.dir() %>% dirname() %>% setwd()

#  choose a set of candidate powers
rm(list = ls())
# library(WGCNA)

asthma_expMat <- readRDS("RDS/asthma_expMat.RDS")
group <- readRDS("RDS/group.RDS")

# Now cluster donors (in contrast to clustering genes later on...)
sampleTree = hclust(dist(asthma_expMat), method = "average")

# Plot the sample tree; tips will be ordered to some extent by library size,
# Look for outliers along this continuum
pdf("result/441.res.SampleTree_outliers.wgcna.pdf", height=48, width=64)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

enableWGCNAThreads()

#First, choose a set of candidate powers to look at
powers = c(1:15)

# Call the network topology analysis function
sft = pickSoftThreshold(asthma_expMat, powerVector = powers, verbose = 5, networkType="signed")
sft$powerEstimate
# Plot the results:
pdf("result//441.res.pickSoftThresholdPlots.wgcna.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

