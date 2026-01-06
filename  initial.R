####### 高级微观经济学小论文代码-1 #######
####### 初步解释性分析 #######
rm(list = ls())  #清空环境
# 导入分析过程中需要的宏包
library(readr) # 用于读入dta/xlsx数据表
library(haven) # 用于读入dta数据表
library(dplyr)  # 用于数据清洗
library(VIM)    # 用于KNN插补
library(mice)   # 用于多重插补
library(psych)  # 用于描述性统计
library(naniar)  # 缺失数据可视化
library(ggplot2) # 绘图
library(readstata13) # 用于读入dta13格式的数据
library(margins) # 用于数据分析
library(gtsummary) # 生成规范的差异性分析表格
# 设置随机种子用于复现
set.seed(123)
# setwd() 
# 导入数据表并初步了解、重命名各列
combined.dta <- read_dta("combined.dta", encoding = "UTF-8")
df <- as.data.frame(combined.dta)
# summary(df)
colnames(df) <- c("hhid","pid","marriage","prefer","year",
                  "province","proveng","provcode","total",
                  "coverage1","depth1","payment2","insurance2","investment2",
                  "credit2","digital1","operation","loan","third",
                  "care","rural","region","income","hhead",
                  "sex","birthyear","education","health","size",
                  "age")
# 区分来源不同的数据
hhpid <- as.data.frame(df[,c(1:2)])
year <- df[,5]
mergebase_prov <- as.data.frame(df[,c(6:8,21,22)])
digital <- as.data.frame(df[,c(9:16)])
invest <- as.data.frame(df[,c(4,17:20)])
baseline <- as.data.frame(df[,c(25,30,26,27,28,3,24,29,23)])
# View(invest)
#重新整理成易于分析和理解的格式
edit_df <- as.data.frame(cbind(year,hhpid,baseline,
                               mergebase_prov,invest,digital))
# 处理情况特殊的loan列缺失
for (j in 1:nrow(edit_df)) {
  if (!is.na(edit_df[j,c("operation")]) && edit_df[j, c("operation")] == 0) {
    # 将下一列的值复制到当前位置
    edit_df[j, c("loan")] <- edit_df[j, c("operation")]
    edit_df[j, c("third")] <- edit_df[j, c("third")]
  }
}
# 调整家庭收入水平的表现形式
#edit_df$income <- log(edit_df$income)
initial_df <- edit_df 
# View(initial_df)
# summary(initial_df)
#存档初步处理表格（含有缺失值）
write_dta(initial_df, "df_initial.dta")


# 数据预处理
# 调整变量类型便于插补
df_factor <- edit_df %>%
  mutate(
    education = as.factor(education),
    health = as.factor(health),
    marriage = as.factor(marriage),
    # income = as.factor(income),
    prefer = as.factor(prefer),
    operation = as.factor(operation),
    loan = as.factor(loan),
    third = as.factor(third),
    care = as.factor(care)
  )
# 查看缺失值概况
# summary(df_factor)
# View(df_factor)
missing_proportion <- sapply(df_factor, function(x) {
  round(sum(is.na(x)) / length(x) * 100, 2)
})
missing_summary <- data.frame(
  variable = names(missing_proportion),
  missing_percent = missing_proportion,
  missing_count = sapply(df_factor, function(x) sum(is.na(x)))
)
# missing_proportion
# missing_summary
# 可视化缺失模式
# vis_miss(df_factor) +
#   theme(axis.text.x = element_text(angle = 45))
# md.pattern(df_factor, plot = TRUE, rotate.names = TRUE)
# 采用KNN插补（k=5）处理缺失值并移除辅助变量
df_missing <- df_factor[,c(7:9,12,18:22)]
df_imputed <- kNN(df_missing, k = 5)
df_imputed <- df_imputed[, 1:ncol(df_missing)]
df_imputed <- df_imputed %>%
  mutate(
    education = as.numeric(education),
    health = as.numeric(health),
    care = as.numeric(care),
    education = ifelse(education < 1, 1, 
                       ifelse(education > 9, 9, education)),
    health = ifelse(health < 1, 1, 
                    ifelse(health > 5, 5, health)),
    care = ifelse(care < 1, 1, 
                  ifelse(care > 5, 5, care))
)
df_final <- df_imputed
# View(df_final)
edit_df$education <- df_final$education
edit_df$income <- df_final$income
edit_df$health <- df_final$health
edit_df$marriage <- df_final$marriage
edit_df$prefer <- df_final$prefer
edit_df$operation <- df_final$operation
edit_df$loan <- df_final$loan
edit_df$third <- df_final$third
edit_df$care <- df_final$care
# View(initial_df)
# summary(edit_df)
# summary(initial_df)
df_NoNA <- edit_df # 存档处理缺失值后的表格
# View(df_NoNA)

# 正式建模分析: 提取分析所需变量
# 描述性统计
df_probit <- df_NoNA[,c("operation","total","age","sex","education",
                        "marriage","health","prefer","loan","third",
                        "care","income","size","rural")]
df_probit$operation <- as.numeric(df_probit$operation)-1
df_probit$marriage <- as.numeric(df_probit$marriage)-1
df_probit$prefer <- as.numeric(df_probit$prefer)-1
df_probit$loan <- as.numeric(df_probit$loan)-1
df_probit$third <- as.numeric(df_probit$third)-1
df_probit <- df_probit %>%
  mutate(
    # 重新编码loan变量
    loan = case_when(
      loan == 1 ~ 2,  # 有经营且有贷款 → 最高分
      loan == 2 ~ 1,  # 有经营但无贷款 → 中等
      loan == 0 ~ 0   # 无经营 → 最低分
    ),    # 反向编码care变量（1-5 → 5-1）
    care = 6 - care  # 现在5=非常关注，1=从不关注
  )
# View(df_probit)
detailed_stats <- describe(df_probit)
# print(detailed_stats)

# 标准化处理（min-max标准化到0-1区间）
df_digital <- df_probit %>%
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
# 采用加权平均将微观指数与宏观指数结合
df_digital <- df_digital %>%
  mutate(
    # 首先标准化宏观指数到0-100（如果尚未标准化）
    macro_index_std = (total - min(total)) / (max(total) - min(total)) * 100,
    digital_participation = 0.5 * micro_digital_std + 0.5 * macro_index_std,
    # 或者使用微观指数作为调节项：综合指标 = 宏观指数 × (1 + 微观参与度)
    # digital_participation_interaction = total * (1 + micro_digital_std/100)
  )
# 检查最终的综合指标
# summary(df_digital$digital_participation)
# describe(df_digital$digital_participation)

# 根据是否从事工商业经营做差异性分析
df_probit$digital_participation <- df_digital$digital_participation
# View(df_probit)
# 首先，定义我们需要进行检验的变量列表
vars_to_test <- c(
  "digital_participation", "age", "sex", "education",
  "marriage", "health", "prefer", "income", "size", "rural"
)

# 使用 dplyr 计算两组的均值
group_means <- df_probit%>%
  group_by(operation) %>%
  summarise(across(all_of(vars_to_test), mean, .names = "mean_{.col}")) %>%
  t() %>% # 转置矩阵，方便后续处理
  as.data.frame() %>%
  rename(mean_group0 = V1, mean_group1 = V2) %>%
  slice(-1) # 删除第一行(operation)

# 使用循环/map函数对每个变量执行T检验，并提取P值
p_values <- purrr::map_df(vars_to_test, ~{
  test_result <- t.test(df_probit[[.x]] ~ df_probit$operation)
  data.frame(variable = .x, p_value = test_result$p.value)
})

# 将均值和P值合并到一个结果表中
results_manual <- cbind(
  variable = vars_to_test,
  group_means,
  p_value = p_values$p_value
)

# 计算差值并增加显著性星号
results_final <- results_manual %>%
  mutate(
    difference = mean_group1 - mean_group0,
    significance = case_when(
      p_value < 0.01 ~ "***",
      p_value < 0.05 ~ "**",
      p_value < 0.10 ~ "*",
      TRUE ~ ""
    ),
    # 格式化数字，方便阅读
    across(c(mean_group0, mean_group1, difference), ~ round(.x, 4)),
    p_value = round(p_value, 10)
  ) %>%
  # 重新排列列的顺序
select(variable, mean_group0, mean_group1, difference, p_value, significance)
# print(results_final)


# 运行基准Probit模型
# 公式：operation ~ digital_participation + 所有控制变量
probit_1 <- glm(operation ~ digital_participation + age + sex + education + marriage + 
                    health + prefer + income + size + rural,
                  family = binomial(link = "probit"), data = df_digital)
summary(probit_1)

# 存储已经处理好的df_probit
write_dta(df_probit, "df_probit.dta")




