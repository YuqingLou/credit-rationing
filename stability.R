####### 高级微观经济学小论文代码-2.1 #######
####### 稳健性分析、工具变量法、中介效应模型 #######
rm(list = ls())  #清空环境
# 导入分析过程中需要的宏包
library(readr) # 用于读入dta/xlsx数据表
library(haven) # 用于读入dta数据表
library(psych)  # 用于描述性统计
library(ggplot2) # 绘图
library(margins) # 用于数据分析
library(gtsummary) # 生成规范的差异性分析表格
library(dplyr) #用于表格构造
library(modelsummary) #用于模型结果归纳
library(marginaleffects) #用于计算边际效益
library(rstudioapi) #用于智能反馈
# 设置随机种子用于复现
set.seed(123)
# setwd() 
# 导入原始数据表和基础回归模型数据表
df_initial.dta <- read_dta("df_initial.dta", encoding = "UTF-8")
df_probit.dta <- read_dta("df_probit.dta", encoding = "UTF-8")
initial <- as.data.frame(df_initial.dta)
probit <- as.data.frame(df_probit.dta)
# View(initial)
# View(probit)
edit <- cbind(probit,initial[,c(15,24:30)])
# View(edit)

# 标准化处理（min-max标准化到0-1区间）
df_digital <- probit %>%
  mutate(third = (third - min(third)) / (max(third) - min(third)),
         #loan = (loan - min(loan)) / (max(loan) - min(loan)),
         care = (care - min(care)) / (max(care) - min(care)),
         across(c(total, age, education, health, income, size), scale)
  )
# 熵值法计算三部分权重
calculate_entropy_weights <- function(data_matrix) {
  # 数据准备
  X <- as.matrix(data_matrix)
  n <- nrow(X)
  m <- ncol(X)
  # 1. 数据非负化（由于我们已经标准化到0-1，可以跳过）
  # 2. 计算第j项指标下第i个样本的比重
  P <- X / rowSums(X)
  # 处理可能的NA/Inf
  P[is.na(P)] <- 0
  # 3. 计算第j项指标的熵值
  k <- 1 / log(n)
  e_j <- -k * colSums(P * log(P + 1e-10), na.rm = TRUE)  # 加一个小数避免log(0)
  # 4. 计算信息熵冗余度
  d_j <- 1 - e_j
  # 5. 计算各项指标的权重
  w_j <- d_j / sum(d_j)
  return(w_j)
}
micro_vars <- df_digital[, c("third","care")]
entropy_weights <- calculate_entropy_weights(micro_vars)
names(entropy_weights) <- c("third", "care")
# print(entropy_weights)
# 使用熵值法权重计算家庭微观指数并标准化
df_digital <- df_digital %>%
  mutate(
    micro_digital = 
      third * entropy_weights["third"] +
      #loan * entropy_weights["loan"] + 
      care * entropy_weights["care"],
    micro_digital_std = 
      (micro_digital - min(micro_digital)) / 
      (max(micro_digital) - min(micro_digital)) * 100
  )
# summary(df_digital$micro_digital)
df_digital <- df_digital %>%
  mutate(
    # 方法一：采用加权平均将微观指数与宏观指数结合
    macro_index_std = (total - min(total)) / (max(total) - min(total)) * 100,
    digital_participation = 0.5 * micro_digital_std + 0.5 * macro_index_std,
    # 方法二：使用微观指数作为调节项：综合指标 = 宏观指数 × (1 + 微观参与度)
    digital_participation_interaction = total * (1 + micro_digital_std/100)
  )
# 方法三：熵权法构建指标
vars <- df_digital[, c("micro_digital_std","macro_index_std")]
entropy_weights <- calculate_entropy_weights(vars)
names(entropy_weights) <- c("micro_digital_std","macro_index_std")
df_digital <- df_digital %>%
  mutate(
    dp_entropy = 
      micro_digital_std * entropy_weights["micro_digital_std"] +
      macro_index_std * entropy_weights["macro_index_std"],
    # 三个指数全部标准化
    digital_participation = (digital_participation - min(digital_participation)) / (max(digital_participation) - min(digital_participation)) * 100,
    dp_entropy = (dp_entropy- min(dp_entropy)) / (max(dp_entropy) - min(dp_entropy)) * 100,
    digital_participation_interaction = (digital_participation_interaction- min(digital_participation_interaction)) / (max(digital_participation_interaction) - min(digital_participation_interaction)) * 100,
    )
edit$micro_std <- df_digital$micro_digital_std
edit$macro_std <- df_digital$macro_index_std
edit$dp_interact <- df_digital$digital_participation_interaction
edit$dp_entropy <- df_digital$dp_entropy
edit$digital_participation <- df_digital$digital_participation
# 重新整理表格为可分析形式
newdata <- edit[,c("operation","digital_participation","dp_interact","dp_entropy","micro_std","macro_std","loan",
                   "third","care","coverage1","depth1",
                   "digital1","age","sex","education","marriage","health",
                   "prefer","income","size","rural","provcode")]
# View(newdata)
summary(newdata)


# 稳健性分析部分
# 首先定义共同的模型和控制变量
controls <- c("age", "sex", "education", "marriage", "health", 
              "prefer", "income", "size")
formula_simple <- as.formula(paste("operation ~", paste(controls, collapse = " + "),"+ as.factor(rural)"))
# 第一部分：分析digital_participant -> micro_std + macro_std
model_micro <- glm(update(formula_simple, . ~ . + micro_std), 
                      data = newdata, family = binomial(link = "probit"))
model_macro <- glm(update(formula_simple, . ~ . + macro_std), 
                   data = newdata, family = binomial(link = "probit"))
modelsummary(
  list("家庭微观指数" = model_micro, "地区宏观指数"= model_macro),
  stars = TRUE,
  gof_map = c("nobs", "r.squared"),
  title = "替换核心解释变量的稳健性检验-1",
  # 计算并添加每个变量的平均边际效应
  # Probit系数本身不易解释，边际效应的解释是“...概率提高X个百分点”
  estimates = "dydx",
  # 只显示核心变量的系数，让表格更简洁
  coef_map = c("micro_std", "macro_std"),
  fmt = 6  # 有效数字
)
# 第二部分：分析micro_std -> micro -> third\loan\care
model_third <- glm(update(formula_simple, . ~ . + third), 
                      data = newdata, family = binomial(link = "probit"))
model_care <- glm(update(formula_simple, . ~ . + care), 
                            data = newdata, family = binomial(link = "probit"))
modelsummary(
  list("第三方账户确认" = model_third, "经济关注度" = model_care),
  stars = TRUE,
  gof_map = c("nobs", "r.squared"),
  title = "替换核心解释变量的稳健性检验-2",
  # 计算并添加每个变量的平均边际效应
  # Probit系数本身不易解释，边际效应的解释是“...概率提高X个百分点”
  estimates = "dydx",
  # 只显示核心变量的系数，让表格更简洁
  coef_map = c("third","care"),
  fmt = 6  # 有效数字
)
# 第三部分：分析macro_std -> total -> coverage1\depth1\digital1
model_coverage <- glm(update(formula_simple, . ~ . + coverage1), 
                      data = newdata, family = binomial(link = "probit"))
model_depth <- glm(update(formula_simple, . ~ . + depth1), 
                   data = newdata, family = binomial(link = "probit"))
model_digitalization <- glm(update(formula_simple, . ~ . + digital1), 
                            data = newdata, family = binomial(link = "probit"))
modelsummary(
  list("覆盖广度" = model_coverage, "使用深度" = model_depth, "数字化程度" = model_digitalization),
  stars = TRUE,
  gof_map = c("nobs", "r.squared"),
  title = "表4：替换核心解释变量的稳健性检验-3",
  # 计算并添加每个变量的平均边际效应
  # Probit系数本身不易解释，边际效应的解释是“...概率提高X个百分点”
  estimates = "dydx",
  # 只显示核心变量的系数，让表格更简洁
  coef_map = c("coverage1", "depth1", "digital1"),
  fmt = 6  # 有效数字
)
# 第四部分：选取不同的方法组装这个综合指标
model_equality <- glm(update(formula_simple, . ~ . + digital_participation), 
                      data = newdata, family = binomial(link = "probit"))
model_interact <- glm(update(formula_simple, . ~ . + dp_interact), 
                   data = newdata, family = binomial(link = "probit"))
model_entropy <- glm(update(formula_simple, . ~ . + dp_entropy), 
                            data = newdata, family = binomial(link = "probit"))
modelsummary(
  list("等权重" = model_equality, "交互" = model_interact, "熵权" = model_entropy),
  stars = TRUE,
  gof_map = c("nobs", "r.squared"),
  title = "替换核心解释变量的稳健性检验-4",
  # 计算并添加每个变量的平均边际效应
  # Probit系数本身不易解释，边际效应的解释是“...概率提高X个百分点”
  estimates = "dydx",
  # 只显示核心变量的系数，让表格更简洁
  coef_map = c("digital_participation", "dp_interact", "dp_entropy"),
  fmt = 6  # 有效数字
)


