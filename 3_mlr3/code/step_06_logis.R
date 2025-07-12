if (!require("mlr3verse")) install.packages("mlr3verse")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("rms")) install.packages("rms")

rm(list = ls())
# 读取数据
df_all <- readRDS("RDS/two_dataset_co_DGEs_exp.RDS")[,-1]

df_all$group <- as.numeric(df_all$group) - 1

ddist <- datadist(df_all)
options(datadist = "ddist")

logis_lrm <- lrm(group ~ ., data = df_all,x=TRUE, y=TRUE)
summary(logis_lrm)

# 矫正曲线
cal1<-calibrate(logis_lrm,method="boot", B=1000)

pdf(file = "result/step_06_logis_calibrate.pdf", width = 10, height = 5)
plot(cal1)
dev.off()

# 列线图
nom1 <- nomogram(logis_lrm,
  fun = plogis,
  fun.at = c(0.001, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99),
  lp = T, # 是否显示线性概率
  funlabel = "Risk of PH"
)

pdf(file = "result/step_06_logis_nom.pdf", width = 10, height = 5)
plot(nom1)
dev.off()

# 全部差异因素回归
logis_glm <- glm(group ~ ., data = df_all, family = binomial)
summary(logis_glm)$coef


anova(object = logis_glm, test = "Chisq")
coef(logis_glm)

exp(coef(logis_glm)) # OR

write.csv(cbind(
  summary(logis_glm)$coef,
  OR=exp(coef(logis_glm))
),'result/step_06_logis_glm.csv')


# library("MASS") # 逐步回归，与全部变量回归的模型比较
# step.model <- step(object = logis_glm, trace = 0)
# summary(step.model)
# 
# anova(object = step.model, test = "Chisq")
# 
# anova(logis_glm, step.model, test = "Chisq")


# S型 模型图
# Fill predicted values using regression model
df_all$pred = predict( logis_glm, df_all , type="response")

ggplot(df_all, aes(x=pred, y=group)) + 
  geom_point(alpha=0.2,size=0.5) +
  stat_smooth(method="glm", color="grey60", se=FALSE,
              method.args = list(family=binomial))+
  ylab("PH (Mild or Severe)")+
  xlab("Predict")+
  scale_y_continuous( breaks = c(0,1))+
  scale_x_continuous( breaks = seq(0.1,1,0.2),limits = c(0.1,0.9))+
  theme_classic()

ggsave(filename = 'result/step_06_logistic_s.pdf',width = 16,height = 8,dpi = 300,units = 'cm')


