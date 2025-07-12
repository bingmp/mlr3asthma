
if(!require(DESeq2)) BiocManager::install("DESeq2")
if (!require("magrittr")) BiocManager::install("magrittr")

this.path::this.dir() %>% dirname() %>% setwd()

### Load raw GALA II count matrix
rm(list = ls())

### Load raw GALA II count matrix
raw_counts <- read.table("rawdata/GSE152004_695_raw_counts.txt.gz", header=T, sep="\t", row.names=1)
raw_counts <- as.matrix(raw_counts)

group <- readRDS('RDS/group.RDS')
raw_counts <- raw_counts[,group$title]

### Vst normalize data
design <- data.frame(row.names = colnames(raw_counts), subject = colnames(raw_counts))
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = design,
  design = ~1)
vst <- varianceStabilizingTransformation(dds)
vstMat <- assay(vst)

saveRDS(vstMat, file = "RDS/vstMat.RDS")

### Size factor normalize data
expr_norm <- estimateSizeFactors(dds)
expr_norm <- counts(expr_norm, normalized=T)

### Remove lowly expressed genes
### 17473 genes
good_genes <- rowSums(raw_counts >= 10) >= (ncol(raw_counts) * .15)
vstMat_rm <- vstMat[good_genes,]
dim(vstMat_rm)

### Load phenotype data
phen <- read.table("rawdata/fixed_phenotype_age_gender_bmi_asthma.csv", header=T, sep=",")
rownames(phen) <- phen$SubjectID

### 695 samples
int.samples <- intersect(rownames(phen), colnames(vstMat_rm))
traits <- phen[int.samples,]

### compute residual matrix (regressing out age and gender)
vstMat_res <- vstMat_rm
for(gene in rownames(vstMat_rm)) {
  fit <- lm(vstMat_rm[gene,]~ traits$age + as.factor(traits$Male))
  res <- residuals(fit)
  vstMat_res[gene,] <- res
}

####check genes
expMat <- t(vstMat_res)
gsg <- goodSamplesGenes(expMat, verbose=3)

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(colnames(expMat)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(expMat)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  expMat = expMat[gsg$goodSamples, gsg$goodGenes]
}

asthma_expMat <- expMat[group$title,]  
asthma_expMat <- asthma_expMat[which(group$condition=="Asthma"),]
saveRDS(asthma_expMat, file = "RDS/asthma_expMat.RDS")

save(expMat,expr_norm,vstMat_res,vstMat_rm ,vstMat,asthma_expMat, file = "RDS/expMat_all.Rdata")
