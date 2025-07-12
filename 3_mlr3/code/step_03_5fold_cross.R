
if(!require("mlr3verse")) install.packages("mlr3verse")
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2")

rm(list = ls())

## 准备 -------------------
if(T){
  # 读取数据
  df_all <- readRDS("RDS/two_dataset_co_DGEs_exp.RDS")[,-1]
  
  # colnames(df_all) <- make.names(names(df_all) )
  # 建立任务
  task <- as_task_classif(df_all, target="group")
  
  # # 数据划分
  set.seed(123)
  split_task <- partition(task, ratio=0.75)
  task_train <-  task$clone()$filter(split_task$train)
  task_test <- task$clone()$filter(split_task$test) 
  
  mymsr = msrs(c("classif.acc", "classif.auc", "classif.npv","classif.ppv","classif.sensitivity", "classif.specificity", "classif.bbrier"))
  
  # 数据预处理
  graph_preb <- po("filter",
                   filter = mlr3filters::flt("find_correlation"),
                   filter.cutoff = 0.3
  ) %>>% # 去除高度相关的列
    po("scale", center =T, scale = F) %>>% # 中心化
    po("removeconstants") %>>% # 去掉零方差变量
    po("encode") %>>%
    po("imputemedian")
  
  # 4象限图函数
  four_quadrant <- function(x, col_quad="gray60", col_text="white") {    
    nx <- length(x)
    sqx <- sqrt(x) 
    df <- data.frame(x=c(sqx[1],-sqx[2],-sqx[3],sqx[4])/2, 
                     y=c(sqx[1],sqx[2],-sqx[3],-sqx[4])/2, 
                     size=sqx, label=x)
    mm <- max(df$size)*1.1
    
    p <- ggplot(data=df, aes(x=x, y=y, width=size, height=size, 
                             group=factor(size))) +
      geom_tile(fill=col_quad) +
      geom_text(aes(label=label), col=col_text, size=5) +
      geom_hline(aes(yintercept=0), size=0.8,color='gray60') +
      geom_vline(aes(xintercept=0), size=0.8,color='gray60') +
      coord_fixed() +
      xlim(c(-mm,mm)) + ylim(c(-mm,mm)) +
      theme(legend.position = "none")+
      labs(y="Truth",x="Prediction")+
      ggthemes::theme_pander() +
      theme(
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20)
      )
    
    return(p)
  }
  
  # 数据准备
  data_5fold <- data.frame(matrix(nrow = 4,ncol = 7))
  colnames(data_5fold) <- c("acc", "auc", "npv","ppv","sensitivity", "specificity", "bbrier")
  rownames(data_5fold) <- c("rf","logis","dtree","knn")
  
}

## 5.1 随机森林模型 -----------------
if(T){
  # 选择随机森林模型
  rf_glr <- as_learner(graph_preb %>>% lrn("classif.ranger", predict_type="prob"))
  rf_glr$id <- "randomForest"
  set.seed(123)
  
  # 5折交叉验证
  rr <- resample(task = task,
                 learner = rf_glr,
                 resampling = rsmp("cv",folds = 5),
                 store_models = T
  )

  # 混淆矩阵：
  rr$prediction()$confusion
  df_pred <- rr$prediction()
  
  autoplot(rr,type = "roc") +
    ggtitle("5-fold cross validation of randomForest") +
    theme(
      plot.title =  element_text(size = 30,hjust = 0.5),
      text = element_text(size = 20),
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 15)
    ) +
    ggsci::scale_color_aaas()+
    ggsci::scale_fill_aaas()
  
  ggsave(filename = 'result/step_03_rf_5fold_roc.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  # 混淆矩阵可视化：
  rr$prediction()$confusion %>% 
    as.vector() %>% 
    .[c(2,1,3,4)] %>% 
    four_quadrant(col_quad="salmon2") + ggtitle("5-fold cross validation of randomForest")
  
  ggsave(filename = 'result/step_03_rf_5fold_four_quadrant.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  # 查看其他结果：
  data_5fold[1,] <- mymsr %>% rr$aggregate()

}

## 5.2 选择 logistics 模型 ---------------
if(T){
  
  # 选择 logistics 模型
  rf_glr <- as_learner(graph_preb %>>% lrn("classif.log_reg", predict_type="prob"))
  rf_glr$id <- "logistic"
  set.seed(123)
  
  # 5折交叉验证
  rr <- resample(task = task,
                 learner = rf_glr,
                 resampling = rsmp("cv",folds = 5),
                 store_models = T
  )

  # 混淆矩阵：
  rr$prediction()$confusion
  df_pred <- rr$prediction()

  autoplot(rr,type = "roc") +
    ggtitle("5-fold cross validation of Logistic") +
    theme(
      plot.title =  element_text(size = 30,hjust = 0.5),
      text = element_text(size = 20),
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 15)
    ) +
    ggsci::scale_color_aaas()+
    ggsci::scale_fill_aaas()
  
  ggsave(filename = 'result/step_03_logis_5fold_roc.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  
  # 混淆矩阵4象限图
  rr$prediction()$confusion %>% 
    as.vector() %>% 
    .[c(2,1,3,4)] %>% 
    four_quadrant(col_quad="salmon2") + ggtitle("5-fold cross validation of Logistics")
  
  ggsave(filename = 'result/step_03_logis_5fold_four_quadrant.pdf',width = 28,height = 16,dpi = 300,units = 'cm')

  # 查看其他结果：
  data_5fold[2,] <- mymsr %>% rr$aggregate()
}

## 5.3 决策树 ------------
if(T){
  # 决策树
  rf_glr <- as_learner(graph_preb %>>% lrn("classif.rpart", predict_type="prob"))
  rf_glr$id <- "decisionTree"
  set.seed(123)
  
  # 5折交叉验证
  rr <- resample(task = task,
                 learner = rf_glr,
                 resampling = rsmp("cv",folds = 5),
                 store_models = T
  )
  
  
  # 混淆矩阵：
  rr$prediction()$confusion
  df_pred <- rr$prediction()
  
  autoplot(rr,type = "roc") +
    ggtitle("5-fold cross validation of decisionTree") +
    theme(
      plot.title =  element_text(size = 30,hjust = 0.5),
      text = element_text(size = 20),
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 15)
    ) +
    ggsci::scale_color_aaas()+
    ggsci::scale_fill_aaas()
  
  ggsave(filename = 'result/step_03_dtree_5fold_roc.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  # 混淆矩阵4象限图
  rr$prediction()$confusion %>% 
    as.vector() %>% 
    .[c(2,1,3,4)] %>% 
    four_quadrant(col_quad="salmon2") + ggtitle("5-fold cross validation of decisionTree")
  
  ggsave(filename = 'result/step_03_dtree_5fold_four_quadrant.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  # 查看其他结果：
  data_5fold[3,] <- mymsr %>% rr$aggregate()
  
}

## 5.4 kknn模型 -------------
if(T){
  # kknn模型
  rf_glr <- as_learner(graph_preb %>>% lrn("classif.kknn", predict_type="prob"))
  rf_glr$id <- "kknn"
  set.seed(123)
  
  # 5折交叉验证
  rr <- resample(task = task,
                 learner = rf_glr,
                 resampling = rsmp("cv",folds = 5),
                 store_models = T
  )
  
  # 混淆矩阵：
  rr$prediction()$confusion
  df_pred <- rr$prediction()
  
  autoplot(rr,type = "roc") +
    ggtitle("5-fold cross validation of K-Nearest Neighbors") +
    theme(
      plot.title =  element_text(size = 30,hjust = 0.5),
      text = element_text(size = 20),
      axis.title = element_text(size = 25),
      axis.text = element_text(size = 15),
      legend.title = element_text(size = 25),
      legend.text = element_text(size = 15)
    ) +
    ggsci::scale_color_aaas()+
    ggsci::scale_fill_aaas()
  
  ggsave(filename = 'result/step_03_knn_5fold_roc.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  # 混淆矩阵4象限图
  rr$prediction()$confusion %>% 
    as.vector() %>% 
    .[c(2,1,3,4)] %>% 
    four_quadrant(col_quad="salmon2") + ggtitle("5-fold cross validation of K-Nearest Neighbors")
  
  ggsave(filename = 'result/step_03_knn_5fold_four_quadrant.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  # 查看其他结果：
  data_5fold[4,] <- mymsr %>% rr$aggregate()
  
  }

write.csv(data_5fold,file = 'result/step_03_data_5fold.csv')
