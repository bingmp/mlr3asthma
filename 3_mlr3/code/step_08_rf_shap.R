
# 加载包
library(tidyverse)
library(foreach)
library(mlr3) 
library(mlr3learners)
library(kernelshap)    # 计算shap值
library(shapviz)       # shap值可视化
library(doFuture)      # 并行计算
library(foreach)


rm(list = ls())
# 读取数据
df_all <- readRDS("RDS/two_dataset_co_DGEs_exp_lyz.RDS")[,-1]
# colnames(df_all) <- make.names(names(df_all) )

# 建立任务
task <- as_task_classif(df_all, target="group")

# # 数据划分
set.seed(123)
split_task <- partition(task, ratio=0.75)

# 数据预处理
graph_preb <- po("filter",
                 filter = mlr3filters::flt("find_correlation"),
                 filter.cutoff = 0.3
) %>>% # 去除高度相关的列
  po("scale", center =T, scale = F) %>>% # 中心化
  po("removeconstants") %>>% # 去掉零方差变量
  po("encode") %>>%
  po("imputemedian")

# 设置学习器
lrn_rf = as_learner(graph_preb %>>% lrn("classif.ranger", predict_type="prob") )
lrn_rf$id <- "Random Forest"
set.seed(123)
# lrn("classif.ranger", id = "Random Forest", seed = 123, predict_type = "prob",  predict_sets = c("train", "test"))
# lrn_rf

lrn_rf$train(task, row_ids = split_task$train)

# 以下是单个模型逐一计算评价指标
mymsr = msrs(c("classif.auc", "classif.acc","classif.sensitivity", "classif.specificity"))

lrn_rf$predict(task, row_ids = split_task$train)$score(measures = mymsr)
lrn_rf$predict(task, row_ids = split_task$test)$score(measures = mymsr)

# 绘制ROC曲线
pred_rf_train = lrn_rf$predict(task, row_ids = split_task$train)
# 画训练集上预测结果的ROC曲线 ?autoplot.PredictionClassif
autoplot(pred_rf_train, type = "roc")


# 开启并行运算
registerDoFuture()
plan(multisession, workers = 50)  # Windows

# 需要解释的观测数据(仅包括特征变量，不包括结局变量)
X_explain = task$data(cols = task$feature_names,rows = split_task$test)

# background dataset，可以是样本子集也可以是全部观测
bg_X = task$data(rows = split_task$train)

pshap_xgb <- kernelshap(
  lrn_rf, X = X_explain, bg_X = bg_X, predict_type = "prob",
  feature_names = task$feature_names, 
  parallel = TRUE, 
  # 可导入ML算法所使用的包
  parallel_args = list(.packages = c("mlr3verse", "mlr3extralearners"))
)
# plan(sequential)  # Windows

# saveRDS(pshap_xgb,"RDS/pshap_xgb_test.RDS")

# 构建可视化对象
sv_xgb <- shapviz(pshap_xgb, interactions = T)
# 条形图（取多个观测的shap值的绝对值再求平均值）
sv_importance(sv_xgb, kind = "bar")
# ggsave(filename = "result//step_08_shap_bar_test.pdf", width = 16, height = 8, dpi = 300, units = "cm")

# 蜂群图（小提琴散点图）
sv_importance(sv_xgb, kind = "beeswarm")
# ggsave(filename = "result/step_08_shap_beeswarm_test.pdf", width = 16, height = 8, dpi = 300, units = "cm")

# 降序输出shap value
# mean(abs(pshap_xgb[["S"]][["Control"]]))  # 计算mean(|SHAP|)
pshap_xgb[["S"]][["Control"]]
pshap_xgb[["S"]][["Asthma"]]

# write.csv(as.data.frame(pshap_xgb[["S"]]),file = 'result/shap_value_test.csv')

sv_importance(sv_xgb, kind = "no")
# write.csv(sv_importance(sv_xgb, kind = "no"),file = 'result/shap_importance_test.csv')

# 依赖图
sv_dependence(sv_xgb, v = "CD44")
# ggsave(filename = "result/step_08_shap_CD44_test.pdf", width = 16, height = 8, dpi = 300, units = "cm")

# 瀑布图（单个样本）
sv_waterfall(sv_xgb, row_id = 1)
# ggsave(filename = "result/step_08_shap_waterfall_test.pdf", width = 16, height = 8, dpi = 300, units = "cm")

# 力图（单个样本）
sv_force(sv_xgb, row_id = 1)
# ggsave(filename = "result/step_08_shap_force_test.pdf", width = 16, height = 10, dpi = 300, units = "cm")
