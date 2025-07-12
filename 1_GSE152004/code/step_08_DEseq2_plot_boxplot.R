
if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("ggplot2")) BiocManager::install("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("ggsignif")) install.packages("ggsignif")

this.path::this.dir() %>%
  dirname() %>%
  setwd()

rm(list = ls())

exp <- readRDS('RDS/vsdmat.RDS')
group <- readRDS('RDS/group_cluster.RDS')
group <- subset(group, cluster!="Th2_low")

# 10、boxplot
boxplot_func <- function(data, group, name_y){
  df <- as.data.frame( t( data )[,c(name_y )] )
  df$group <- group$cluster
  colnames(df)[1] <- name_y
  p <- ggplot()+
    geom_boxplot(data = df,aes(y =  eval(parse(text = name_y)),
                               x = group, 
                               color=group ),
                 alpha=0.2, show.legend = F )+
    geom_jitter(data = df,aes(y =  eval(parse(text = name_y)),
                              x = group, 
                              color=group ),
                show.legend = F )+
    xlab("") + ylab(name_y)+
    ggsignif::geom_signif(comparisons = list(c("Control","Th2_high")),
                          step_increase = 0.3,
                          map_signif_level = T,
                          test = wilcox.test )+
    theme_classic()+
    theme(
      axis.title = element_text(size=15),
      axis.text = element_text(size=10)
    )
  return(p)
  
}

boxplot_func( data = exp, group = group, c("CD44"))

ggsave(filename = 'result/CD44_boxplot.png', height = 6, width = 6, units = 'cm',dpi = 300)


