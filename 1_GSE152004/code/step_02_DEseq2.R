if (!require("magrittr")) BiocManager::install("magrittr")
if (!require("ggplot2")) BiocManager::install("ggplot2")
if (!require("pROC")) install.packages("pROC")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("FactoMineR")) install.packages("FactoMineR")
if (!require("factoextra")) install.packages("factoextra")
if (!require("DESeq2")) install.packages("DESeq2")

this.path::this.dir() %>% dirname() %>% setwd()

# 1、DEseq2
if(F){
  
  rm(list = ls())
  group <- readRDS("RDS/group_cluster.RDS")
  raw_counts <- read.table("rawdata/GSE152004_695_raw_counts.txt.gz", header=T, sep="\t", row.names=1)
  raw_counts <- as.matrix(raw_counts)

  good_genes <- rowSums(raw_counts >= 10) >= (ncol(raw_counts) * .15)
  
  group <- subset(group, cluster!="Th2_low")
  
  exp <- raw_counts[good_genes, group$title]

  colnames(exp)==group$title
  condition <- group$cluster
  head(exp[,1:5])

  # 2.2 构建分组数据框 colData
  colData <- as.data.frame( group$cluster )
  rownames(colData) <- colnames(exp)  # 设置行名为样本名信息
  colnames(colData) <- 'condition' # 列名
  colData$condition <- factor(colData$condition)
  head(colData)
  
  # 三、DEseq2处理 --------------------------------------------------------------
  # 3.1 构建DEseq2对象
  dds <- DESeqDataSetFromMatrix(countData = exp,
                                colData = colData,
                                design = ~ condition)
  dds
  
  # 3.2 DEseq2
  library(BiocParallel)
  dds <- DESeq(dds,parallel = T)

  save(group,dds, file = 'RDS/511_control_Th2_high_dds.Rdata')
  }

# 3、标准化
if(F){
    # 方差稳定变换，The variance stabilizing transformation(vst)
    vsd <- vst(dds)
    class(vsd)
    vsdmat <- assay(vsd) # 提取转化后的矩阵
    plotPCA(vsd, intgroup=c('condition'))
    saveRDS(vsdmat,file = "RDS/vsdmat.RDS")
  }
  

# 4、提取DEseq2分析的结果 ---------------------------------------------------------
if(F){
  rm(list = ls())
  load('RDS/511_control_Th2_high_dds.Rdata')

  table(group$cluster)
  # 5.1 set contrast. contrast = c('condition','test','control')
  res <- results(dds, contrast = c("condition", "Th2_high", "Control") )
  
  temoutput <- as.data.frame( res )
  temoutput <- temoutput[order(temoutput$log2FoldChange),]
  
  temoutput$change <- as.factor(
    ifelse(
      temoutput$pvalue < 0.05 & abs(temoutput$log2FoldChange)> log2(1.5),
      ifelse(
        temoutput$log2FoldChange> log2(1.5),'up','down'),
      'stable'))
  
  table(temoutput$change) # 查看基因上、下调情况
  
  temoutput$ID <- rownames(temoutput)
  temoutput <- temoutput[order(temoutput$log2FoldChange),]
  tem <- subset(temoutput, change!="stable")
  
  write.csv(temoutput,row.names = F,file = 'control_Th2_high_all_diff.csv')
  save(temoutput, tem,file = 'RDS/Th2_high_vs_Control_diff.Rdata')

# 5、热图 --------------------------------------------------------------------

  load('./control_Th2_high/vsdmat.Rdata')
  
  choose_gene <- c(head(rownames(tem ), 25), tail(rownames(tem ), 25) )
  # choose_gene <-  c("ALOX15", "SLC39A8","SLC40A1","CA2","C3")
  
  choose_matrix1 <- vsdmat[choose_gene,]
  
  annotation_col = data.frame(
    group = as.vector(group$cluster)
  )
  rownames(annotation_col) <- colnames(choose_matrix1)
  
  d <- apply(choose_matrix1, 2, function(x){log2(x+1)})
  
  p <- pheatmap::pheatmap(d, scale = 'row', cluster_cols = T,
                          cutree_rows = 2,cutree_cols = 2,
                          show_colnames = F,
                          silent = F,
                          annotation_col = annotation_col ) 

  ggsave(p,height = 8,width = 8,filename='./control_Th2_high/heatmap_top25_up_down.PDF',) # 输出
}

# 6、火山图 -------------------------------------------------------------------
if(F){
  logFC_cutoff <- log2(1.5)
  
  df <- na.omit(temoutput)
  # 设置火山图的标题
  this_tile <- paste('Cutoff for logFC is ',round(logFC_cutoff,3),
                     '\nThe number of up gene is ',nrow(tem[tem$change=='up',]),
                     '\nThe number of down gene is ',nrow(tem[tem$change=='down',]))
  
  temoutput$label <- NA
  for (i in co_DEGs) {
    temoutput[which(temoutput$ID==  i ),'label'] <- i
  }
  
  # 画火山图

  p <- ggplot(data= temoutput,
              aes(x = log2FoldChange,
                  y = -log10(pvalue),   # 这里将 pvalue 取负对数
                  color = change ) ) +
    geom_point(alpha=0.4,size=1) +     #  绘制点图
    ggrepel::geom_text_repel(data = temoutput,
                             aes(x = log2FoldChange, 
                                 y = -log10(pvalue), 
                                 label = label),
                             max.overlaps = 30,
                             size = 3, 
                             box.padding = unit(0.5, "lines"),
                             point.padding = unit(0.8, "lines"), 
                             segment.color = "black", 
                             show.legend = FALSE) +
    xlab("log2 fold change") + 
    ylab("-log10 pvalue") +    # 轴标签
    # ggtitle(this_tile) +
    theme_bw(base_size=10) +
    theme(
      axis.text  = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 15, hjust = 0.5),
      plot.title = element_text(size = 8, hjust = 0.5)
    ) +
    scale_color_manual(values=c('blue','black','red'))   #设定颜色
  p
  ggsave( p, height = 6,width = 8,filename='./control_Th2_high/511_control_Th_high_volcano.PDF') # 输出
  
}

# 7、PCA
if(F){
  pca <- as.data.frame(t(vsdmat ))
  pca$group <-  group$cluster
  
  library(FactoMineR)
  library(factoextra)
  df.pca <- FactoMineR::PCA( pca[,-ncol(pca)], graph = F,scale=F) 
  p <- fviz_pca_ind(df.pca,
                    geom = "point",
                    mean.point  = T,
                    repel = F,
                    # palette   = c("#9ed048", "#1685a9"),
                    col.ind = pca[,"group"]
                    )+
    scale_color_discrete("")+
    scale_shape_discrete("")
  p
  
  # ggsave(p,filename = 'pca.png',height = 10,width = 12,units = 'cm',dpi = 300)
  
}

# 8、ROC
if(F){

  roc_func <- function(data,name_gene){
    roc.m <- data.frame(data= data[name_gene, ],
                        group=group$cluster
    )
    colnames(roc.m) <- c(name_gene,'group')
    roc.fit <- pROC::roc(group~.,data=roc.m,aur=TRUE,
                         ci=TRUE, # 显示95%CI
                         # percent=TRUE, ##是否需要以百分比显示
                         smooth=F,
                         levels=c('Control','Th2_high')
    )
    
    if(!require(ggplot2)) install.packages("ggplot2")
    p <- pROC::ggroc(roc.fit, legacy.axes = T, color='red') + # 将X轴改为0-1，（默认是1-0）
      geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
                   color="darkgrey", linetype=4) +
      annotate("text",x=0.75,y=0.36,size=4,
               label=paste(roc.fit[["predictor.name"]],"- AUC = ", round(roc.fit$auc,3) ) )+
      theme_bw()+ # 设置背景
      theme(
        axis.title = element_text(size=10),
        axis.text = element_text(size=8),
        legend.text = element_text(size=6),
        legend.title = element_blank()
      )+
      scale_color_discrete()
    return(p)
  }
  
  roc_func(vsdmat,"MUC5B")
  roc_func(vsdmat,"CCL26")
  roc_func(vsdmat,"SERPINB2")
  roc_func(vsdmat,"POSTN")
  roc_func(vsdmat,"PHLDB2")
  roc_func(vsdmat,"CSTA")
  roc_func(vsdmat,"CDC42EP5")

  for (i in co_DEGs) {
    ggsave( roc_func(vsdmat,i), height = 10, width = 10 ,units = 'cm',dpi = 300,
            filename = paste0('./control_Th2_high/roc/',i,'_roc.png'))
  }
  
  p
  # ggsave(p,filename = 'roc.png',height = 10,width = 10,units = 'cm',dpi = 300)
  
}

# 9、correlatin
if(F){
  
  if(!require(ggpubr)) install.packages("ggpubr")
  plot_func <- function(data,name_x,name_y){
    
    library(ggplot2)
    df <- as.data.frame( t( data )[,c(name_x,name_y)] )
    df$group <- group$cluster
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
  

  
  plot_func(vsdmat,name_x="CCL26",name_y="IL1RL1")
  plot_func(vsdmat,name_x="S100A9",name_y="IL6")
  
  p <-   plot_func(vsdmat,name_x="CCL26",name_y="IL1RL1")
  p
  
  ggsave(p,filename = './control_Th2_high/il6correlation.png',height = 10,width = 10,units = 'cm',dpi = 300)
  
}

# 10、boxplot
if(F){
  library(ggplot2)
  boxplot_func <- function(data,name_y){
    df <- as.data.frame( t( data )[,c(name_y )] )
    df$group <- group$cluster
    colnames(df)[1] <- name_y
    p <- ggplot(data = df,aes(y =  eval(parse(text = name_y)),
                              x = group, 
                              color=group ))+
      geom_boxplot(alpha=0.2, show.legend = F )+
      geom_jitter(show.legend = F )+
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

  boxplot_func( vsdmat ,c("CCL26"))
  boxplot_func( vsdmat ,"VWF")
  boxplot_func( vsdmat ,"PHLDB2")
  boxplot_func( vsdmat ,"PCSK6")
  
}
