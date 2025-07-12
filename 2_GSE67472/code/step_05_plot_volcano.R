if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyverse)) install.packages("tidyverse")

this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())

load("RDS/th2_h_vs_con_diff_all.Rdata")
co_DEGs <-c("SCGB1A1","CSTA","CYP24A1","CCL26", "TFF3","KCNJ16","C3","CD44") # 凌耀政

logFC_cutoff <- log2(1.5)
df <- na.omit(temoutput)

# set plot title
this_tile <- paste('Cutoff for logFC is ',round(logFC_cutoff,3),
                   '\nThe number of up gene is ',nrow(tem[tem$change=='up',]),
                   '\nThe number of down gene is ',nrow(tem[tem$change=='down',]))

temoutput$label <- NA
for (i in co_DEGs) {
  temoutput[which(temoutput$ID==  i ),'label'] <- i
}


# plot volcano
p <- ggplot(data= temoutput,aes(x = logFC, 
                                y = -log10(P.Value), 
                                color = change ) ) +
  geom_point( alpha=0.4, size=1) + 
  ggrepel::geom_text_repel(data = temoutput,
                           aes(x = logFC, 
                               y = -log10(P.Value), 
                               label = label),
                           max.overlaps = 30,
                           size = 3, 
                           box.padding = unit(0.5, "lines"),
                           point.padding = unit(0.8, "lines"), 
                           segment.color = "black", 
                           show.legend = FALSE) +
  xlab("log2 fold change") + 
  ylab("-log10 pvalue") + 
  # ggtitle(this_tile) +
  ggtitle("DEGs volcano of GSE67472") +
  theme_bw(base_size=10) +
  # xlim(-20,20)+
  theme(
    axis.text  = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 15, hjust = 0.5),
    plot.title = element_text(size = 15, hjust = 0.5)
  ) +
  scale_color_manual(values=c('blue','black','red'))   #设定颜色

ggsave(p, width = 4, height = 4, filename='result/GSE67472_volcano.PDF') # 输出
