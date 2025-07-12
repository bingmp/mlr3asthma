if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("ggplot2")) BiocManager::install("ggplot2")
if (!require("ggpubr")) install.packages("ggpubr")


this.path::this.dir() %>%
  dirname() %>%
  setwd()

rm(list = ls())

exp <- readRDS('RDS/vsdmat.RDS')
group <- readRDS('RDS/group.RDS')

# 9、correlation
plot_func <- function(data, group, name_x, name_y) {
  library(ggplot2)
  df <- as.data.frame(t(data)[, c(name_x, name_y)])
  df$group <- group$cluster
  p <- ggplot(data = df, aes(
    x = eval(parse(text = name_x)),
    y = eval(parse(text = name_y)),
    color = group
  )) +
    geom_point(size = 0.8, show.legend = T) +
    geom_smooth(method = "lm", formula = "y ~ x", show.legend = F) +
    ggpubr::stat_cor(method = "pearson", color = "black", size = 4, p.accuracy = 0.001) +
    xlab(name_x) +
    ylab(name_y) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.text = element_text(size = 10)
    ) +
    scale_color_discrete("")
  return(p)
}

plot_func(data = exp, group = group, name_x = "CD44", name_y = "CSTA")
ggsave( filename = "result/511_control_Th2_high_correlation_CD44.png", height = 10, width = 10, units = "cm", dpi = 300)

