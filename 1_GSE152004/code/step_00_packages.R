if(!require(BiocManager)) install.packages("BiocManager")
if(!require(GEOquery)) BiocManager::install("GEOquery")
if(!require(DESeq2)) BiocManager::install("DESeq2")
if(!require(WGCNA)) BiocManager::install("WGCNA")
if (!require("magrittr")) BiocManager::install("magrittr")

rm(list = ls())

this.path::this.dir() %>% dirname() %>% setwd()