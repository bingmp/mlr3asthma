
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(tidyverse)) install.packages("tidyverse")

this.path::this.dir() %>% dirname() %>% setwd()

rm(list = ls())

exp <- readRDS('RDS/norm_exp.RDS')
group <- readRDS('RDS/group.RDS')


# 7、correlatin
plot_func <- function(data, group, name_x, name_y){
  if(!require(ggpubr)) install.packages("ggpubr")
  library(ggplot2)
  df <- as.data.frame( t( data )[,c(name_x,name_y)] )
  df$group <- group$condition
  p <- ggplot(data = df,aes(x =  eval(parse(text = name_x)),
                            y = eval(parse(text = name_y)), 
                            color=group ))+
    geom_point(size=0.8,show.legend = T  )+
    geom_smooth(method = 'lm', formula = 'y ~ x',show.legend = F)+
    ggpubr::stat_cor(method = "pearson",color='black',size=4,p.accuracy = 0.001)+
    xlab(name_x) + ylab(name_y)+ 
    theme_classic()+
    theme(
      axis.title = element_text(size=12),
      axis.text = element_text(size=10),
      legend.text = element_text(size=10)
    )+
    scale_color_discrete("")
  return(p)
}

p <- plot_func(exp, group, name_x="CD44",name_y="CSTA")

ggsave(p, filename = 'result/CD44_correlation.png', height = 6, width = 6, units = 'cm',dpi = 300)

