
if(!require(GEOquery)) BiocManager::install("GEOquery")
if (!require("magrittr")) BiocManager::install("magrittr")

this.path::this.dir() %>% dirname() %>% setwd()

### Load GALA II group data
rm(list = ls())
gset <- getGEO('GSE152004', getGPL = F, AnnotGPL = F,destdir = 'rawdata')
group <- pData( gset[[1]] )
group$condition <- group$`asthma status:ch1`
group$condition <- factor(group$condition,
                          levels =c("asthmatic","healthy control"), 
                          labels = c("Asthma","Control"))  
table(group$condition)

group <- group[order(group$condition),]
saveRDS(group,file = 'RDS/group.RDS')
