####### 高级微观经济学小论文代码-3 #######
####### 异质性分析与深入 #######
rm(list = ls())  #清空环境
# 导入分析过程中需要的宏包
library(readr) # 用于读入dta/xlsx数据表
library(haven) # 用于读入dta数据表
library(psych)  # 用于描述性统计
library(ggplot2) # 绘图
library(margins) # 用于数据分析
library(gtsummary) # 生成规范的差异性分析表格
library(dplyr) # 用于表格构造、数据分组
library(modelsummary) # 用于模型结果归纳
library(marginaleffects) # 用于计算边际效益
library(rstudioapi) # 用于智能反馈
library(AER) # 用于工具变量检验
library(GJRM) # 用于IV-probit模型搭建
library(purrr) # 用于回归循环设计

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
# 修正income系数全为0的问题
newdata$income <- newdata$income/10000
#View(newdata)
#summary(newdata)


# 预处理：创建分组ID-地区平均数字金融参与度（除样本户）
newdata <- newdata %>%
  mutate(
    # 假设 prov_code 是省份代码, rural 是城乡虚拟变量
    # 如果没有 prov_code，请替换为你数据中代表省份的变量名
    # paste0 用于将两个变量拼接成一个字符串ID
    regiongroup = paste0(provcode, "_", rural) 
  )
newdata <- newdata %>%
  group_by(regiongroup) %>%
  mutate(
    # 1. 计算该组内的总参与度总和
    group_sum = sum(digital_participation, na.rm = TRUE),
    # 2. 计算该组内的样本数量
    group_n = n(),
    # 3. 计算 Leave-one-out Mean (核心步骤)
    # 公式: (组总和 - 自己) / (组人数 - 1)
    iv_region = (group_sum - digital_participation) / (group_n - 1)
  ) %>%
  ungroup() %>%
  # 去除那些极其罕见的孤立样本（如果某组只有1个人，分母为0，会报错）
  filter(group_n > 1)
describe(newdata)

# Step 1:数据分组 (四类指标: rural, education, income, prefer)
# 按照rural分类 (rural==0 城镇 ==1 乡村)
nd_urban <- newdata %>% 
  filter(rural == 0)
nd_rural <- newdata %>% 
  filter(rural == 1)
# 按照education分类 (education>=4 超过人均教育水平 <=3 低于人均水平)
nd_highedu <- newdata %>% 
  filter(education >=4)
nd_lowedu <- newdata %>% 
  filter(education <=3)
# 按照income分类 (income>=5.182 超过收入中位数 income<5.182 低于收入中位数)
nd_highinc <- newdata %>% 
  filter(income >=5.182)
nd_lowinc <- newdata %>% 
  filter(income < 5.182)
# 按照prefer分类 (prefer==1 风险偏好 prefer==0 风险厌恶)
nd_like <- newdata %>% 
  filter(prefer==1)
nd_dislike <- newdata %>% 
  filter(prefer==0)
# 将所有子样本数据框放入一个命名的列表中
data_list <- list(
  "乡村样本" = nd_rural,
  "城镇样本" = nd_urban,
  "低收入组" = nd_lowinc,
  "高收入组" = nd_highinc,
  "低学历组" = nd_lowedu,
  "高学历组" = nd_highedu,
  "风险偏好组" = nd_like,
  "风险规避组" = nd_dislike
)

# Step 2:分别进行IV_2SLS回归并汇总结果





