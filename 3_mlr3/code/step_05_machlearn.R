if (!require("tidyverse")) install.packages("tidyverse")
if (!require("mlr3verse")) install.packages("mlr3verse")
if (!require("ranger")) install.packages("ranger")
if (!require("kknn")) install.packages("kknn")
if (!require("ggsci")) install.packages("ggsci")
showtext::showtext_auto()

rm(list = ls())

## 1 读取数据 --------
df_all <- readRDS("RDS/two_dataset_co_DGEs_exp.RDS")[,-1]
task <- as_task_classif(df_all, target = "group") # 确定任务

# 学习算法准备
if (T) {
  
  # 图学习器
  graph_preb <- po("filter",
                   filter = mlr3filters::flt("find_correlation"), filter.cutoff = 0.3
  ) %>>% # 去除高度相关的列
    po("scale", center =T, scale = F) %>>% # 中心化
    po("removeconstants") %>>% # 去掉零方差变量
    po("encode") %>>%
    po("imputemedian")
  
  mymsr = msrs(c("classif.acc", "classif.auc", "classif.npv","classif.ppv","classif.sensitivity", "classif.specificity", "classif.bbrier"))
  
  # 随机森林
  rf_glr <- as_learner(graph_preb %>>% lrn("classif.ranger", predict_type = "prob"))
  rf_glr$id <- "randomForest"
  
  # 逻辑回归
  log_glr <- as_learner(graph_preb %>>% lrn("classif.log_reg", predict_type = "prob"))
  log_glr$id <- "logistic"
  
  # 决策树
  tree_glr <- as_learner(graph_preb %>>% lrn("classif.rpart", predict_type = "prob"))
  tree_glr$id <- "decisionTree"

  # k近邻
  kknn_glr <- as_learner(graph_preb %>>% lrn("classif.kknn", predict_type = "prob"))
  kknn_glr$id <- "kknn"
}

# 3 训练模型

# 3.1 4种机器学习模型
if(T){
  # 3.1 建立多个模型
  design <- benchmark_grid(tasks = task, # 任务
                           learners = list(rf_glr, log_glr, tree_glr, kknn_glr), # 算法
                           resampling = rsmp("cv", folds = 5) # 抽样方法
                           ) # design
  
  set.seed(123)
  bmr <- benchmark(design) # 执行基准测试

  # 可视化: 对比性能
  autoplot(bmr, type = "roc") + # ROC 曲线
    ggtitle("")+
    labs(fill = "Learner", color = "Learner") +
    ggsci::scale_color_aaas()+
    ggsci::scale_fill_aaas()+
    theme_classic()+
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 10)
    ) 
  
  ggsave(filename = 'result/step_05_learner_4models_roc.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  autoplot(bmr, measure = msr("classif.auc")) + # AUC 箱线图
    theme(
      text = element_text(size = 0),
      plot.title = element_text(size = 30, hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size =15)
    ) +
    ggtitle("AUC of Learner") +
    ylab("AUC")
  
  ggsave(filename = 'result/step_05_learner_4models_boxplot.pdf',width = 20,height = 10,dpi = 300,units = 'cm')
  
  # bmr$aggregate() # 默认输出
  bmr_msrs <- mymsr %>% bmr$aggregate() %>% 
    subset(select=c('learner_id', # 'resampling_id', 'iters', 
                    "classif.acc", "classif.auc", "classif.npv","classif.ppv","classif.sensitivity", "classif.specificity", "classif.bbrier"))
  bmr_msrs
  readr::write_excel_csv(bmr_msrs,file = 'result/step_05_bmr_4models_msrs.csv')
  
}

# 3.2 随机森林与logistics回归
if(T){
  # 建立多个模型
  design <- benchmark_grid(tasks = task, # 任务
                           learners = list(rf_glr, log_glr), # 算法
                           resampling = rsmp("cv", folds = 5) # 抽样方法
  ) # design
  
  set.seed(123)
  bmr <- benchmark(design) # 执行基准测试

  ## 3.2 可视化: 对比性能 ----------------
  autoplot(bmr, type = "roc") + # ROC 曲线
    ggtitle("")+
    labs(fill = "Learner", color = "Learner") +
    ggsci::scale_color_aaas()+
    ggsci::scale_fill_aaas()+
    theme_classic()+
    theme(
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 10)
    ) 
  
  ggsave(filename = 'result/step_05_logis_rf_roc.pdf',width = 28,height = 16,dpi = 300,units = 'cm')
  
  autoplot(bmr, measure = msr("classif.auc")) + # AUC 箱线图
    theme(
      text = element_text(size = 0),
      plot.title = element_text(size = 30, hjust = 0.5),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 20),
      legend.text = element_text(size =15)
    ) +
    ggtitle("AUC of Learner") +
    ylab("AUC")
  
  ggsave(filename = 'result/step_05_logis_rf_boxplot.pdf',width = 20,height = 10,dpi = 300,units = 'cm')
  
  # bmr$aggregate() # 默认输出
  bmr_msrs <- mymsr %>% bmr$aggregate() %>% 
    subset(select=c('learner_id', # 'resampling_id', 'iters', 
                    "classif.acc", "classif.auc", "classif.npv","classif.ppv","classif.sensitivity", "classif.specificity", "classif.bbrier"))
  bmr_msrs
  readr::write_excel_csv(bmr_msrs,file = 'result/step_05_bmr_logis_rf_msrs.csv')
  
}
