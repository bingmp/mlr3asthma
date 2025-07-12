
if(!require(GEOquery)) BiocManager::install("GEOquery")
this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())

gset <- getGEO(filename = 'rawdata/GSE67472_series_matrix.txt.gz',
               getGPL = F, AnnotGPL = F,destdir = '.')

raw_exp <- exprs(gset)

group <- pData( phenoData(gset) ) ; colnames(group )

group$condition <- factor(group$`disease state:ch1`,
                          levels = c("healthy","asthma"),
                          labels = c("Control","Asthma"))

group$`th2 group:ch1`<- ifelse(is.na(group$`th2 group:ch1`),"Control",group$`th2 group:ch1`)
table(group$`th2 group:ch1`)

group$cluster <- group$`th2 group:ch1`
group[which(group$cluster=='high'),'cluster'] <- "Th2_high"
group[which(group$cluster=='low'),'cluster'] <- "Th2_low"
table(group$condition)

table(group$cluster)
group <- group[order(group$condition),]

raw_exp <- raw_exp[, group$geo_accession]

group$age <- as.numeric(group$`age:ch1` )
group$gender <- as.factor(group$`gender:ch1`)
summary(group[,c("age","gender")])

saveRDS(group, file = 'RDS/group.RDS')
saveRDS(raw_exp, file = 'RDS/raw_exp.RDS')
