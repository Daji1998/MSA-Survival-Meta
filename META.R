getwd()
#进展轨迹图----

library(ggplot2)
library(dplyr)
library(ggprism)
library(survival)
library(survminer)
library(meta)
library(cowplot)
"#0099B3"
"#55A0FB"
"#3861ac"

data <- read.csv("progression.csv")
head(data,20)
str(data)

p1 <- ggplot(data, aes(x = duration, y = UMSARS_I)) +
  geom_line(aes(group = study),color = "#3861ac", size = 1, alpha = 0.8) +   # 所有队列用黑色线
  geom_point(color = "#3861ac", alpha = 0.8) +  # 黑色点
  geom_smooth(aes(weight = sample), 
              method = "lm", se = TRUE, color = "#EC0000", size = 1.2) + 
  labs(x = "Disease duration (years)",
       y = "UMSARS I score",
       size = "Sample size") +
  theme_prism(base_size = 14)

p2 <- ggplot(data, aes(x = duration, y = UMSARS_II)) +
  geom_line(aes(group = study),color = "#3861ac", size = 1, alpha = 0.8) +   # 所有队列用黑色线
  geom_point(color = "#3861ac", alpha = 0.8) +  # 黑色点
  geom_smooth(aes(weight = sample), 
              method = "lm", se = TRUE, color = "#EC0000", size = 1.2) + 
  labs(x = "Disease duration (years)",
       y = "UMSARS II score",
       size = "Sample size") +
  theme_prism(base_size = 14)

p3 <- ggplot(data, aes(x = duration, y = UMSARS_total)) +
  geom_line(aes(group = study),color = "#3861ac", size = 1, alpha = 0.8) +   # 所有队列用黑色线
  geom_point(color = "#3861ac", alpha = 0.8) +  # 黑色点
  geom_smooth(aes(weight = sample), 
              method = "lm", se = TRUE, color = "#EC0000", size = 1.2) + 
  labs(x = "Disease duration (years)",
       y = "UMSARS total score",
       size = "Sample size") +
  theme_prism(base_size = 14)
plot_grid(p1,p2,p3,ncol = 3)


#年评分进展meta----
#风险偏倚评估
library(robvis)
data <- read.csv("rob2.csv")
data
p1 <- rob_traffic_light(data, tool = "ROB2", psize = 15)
p1 <- p1 +
  coord_flip() +
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15)
  )
data$Weight <- 1
p2 <- rob_summary(data, tool = "ROB2")
p2 <- p2 +
  theme(
    text = element_text(size = 16),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15)
  )
plot_grid(p2,p1,ncol = 1)

ggsave(
  "RoB2_traffic_light.tiff",
  plot = p1,
  width = 7.5,
  height = 4,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
#导入数据，至少有ID，mean，sd三列（或LCI，UCI四列）
data <- read.csv("meta.csv")
head(data)
str(data)
data$rigion <- as.factor(data$rigion)
data$type <- as.factor(data$type)
Idata <- subset(data, !is.na(UMSARS_II_diff) & study != "Meissner WG, 2012")
Idata <- subset(data, !is.na(UMSARS_I_diff))
Idata <- data
str(Idata)
metaI <- metagen(
  TE = UMSARS_I_diff,
  lower = UMSARS_I_lci,
  upper = UMSARS_I_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE,
  subgroup = rigion  #分组变脸
)
metaII <- metagen(
  TE = UMSARS_II_diff,
  lower = UMSARS_II_lci,
  upper = UMSARS_II_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE,
  subgroup = rigion  #分组变脸
)
metatotal <- metagen(
  TE = UMSARS_total_diff,
  lower = UMSARS_total_lci,
  upper = UMSARS_total_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE,
  subgroup = rigion #分组变脸
)
# 绘制随机效应森林图----
forest(
  metaI,
  sortvar = TE,           # 根据效应量排序
  xlab = "UMSARS I score change",
  leftcols = c("studlab","sample"),  # 左侧显示研究名和分组
  leftlabs = c("Study", "NO"),
  rightcols = c("effect", "ci"),  # 只显示效应量和置信区间
  rightlabs = c("Effect", "95% CI"),
  col.study = "black",    # 研究名颜色
  col.square = "#3861ac",   # 方块颜色
  col.diamond = "red",    # 合并效应颜色
  col.diamond.lines = "red",
  print.subgroup.name = F, # 打印分组名
  digits = 2,
  colgap.left = unit(45, "mm")
)

forest(
  metaII,
  sortvar = TE,           # 根据效应量排序
  xlab = "UMSARS II score change",
  leftcols = c("studlab","sample"),  # 左侧显示研究名和分组
  leftlabs = c("Study", "NO"),
  rightcols = c("effect", "ci"),  # 只显示效应量和置信区间
  rightlabs = c("Effect", "95% CI"),
  col.study = "black",    # 研究名颜色
  col.square = "#3861ac",   # 方块颜色
  col.diamond = "red",    # 合并效应颜色
  col.diamond.lines = "red",
  print.subgroup.name = F, # 打印分组名
  digits = 2,
  colgap.left = unit(40, "mm")
)

forest(
  metatotal,
  sortvar = TE,           # 根据效应量排序
  xlab = "UMSARS total score change",
  leftcols = c("studlab","sample"),  # 左侧显示研究名和分组
  leftlabs = c("Study", "NO"),
  rightcols = c("effect", "ci"),  # 只显示效应量和置信区间
  rightlabs = c("Effect", "95% CI"),
  col.study = "black",    # 研究名颜色
  col.square = "#3861ac",   # 方块颜色
  col.diamond = "red",    # 合并效应颜色
  col.diamond.lines = "red",
  print.subgroup.name = F, # 打印分组名
  digits = 2,
  colgap.left = unit(40, "mm")
)

library(meta)
library(gridExtra)
library(gridGraphics)
library(grid)
library(cowplot)
# 函数：把 forest 图保存成 grob
grab_forest <- function(metaobj, xlab, leftcols = TRUE, col.square = "black") {
  grid.grabExpr(
    forest(metaobj,
           sortvar = metaI$TE,
           xlab = xlab,
           leftcols = if (leftcols) c("studlab","sample") else FALSE,
           leftlabs = if (leftcols) c("Study","N") else NULL,
           rightcols = c("effect", "ci"),  # 只显示效应量和置信区间
           rightlabs = c("Effect", "95% CI"),
           col.square = col.square,
           col.diamond = "red",
           print.subgroup.name = F, # 打印分组名
           digits = 2,
           colgap.left = unit(45, "mm"))
  )
}

# 分别生成三个 grob
g1 <- grab_forest(metaI, "UMSARS I", leftcols = TRUE, col.square = "#3861ac")
g2 <- grab_forest(metaII, "UMSARS II", leftcols = F, col.square = "#3861ac")
g3 <- grab_forest(metatotal, "UMSARS Total", leftcols = F, col.square = "#3861ac")

# 横向拼接
grid.arrange(g1, g2, g3, ncol = 3)

plot_grid(g1,g2,g3,ncol = 3)
#敏感性分析及森林图----
# forest(metainf(hrmeta,pooled = "random"),layout = "RevMan5")
forest(metainf(metatotal,pooled = "random"),layout = "meta",
       col.bg = rgb(248,118,109,maxColorValue = 255),        # 设置方块颜色
       col.border = "navy",     # 设置边框线颜色
       xlim = c(-5, 15),
       text.common = "Random effect model",
       col.diamond = "darkblue",          # 汇总效应的菱形颜色
       col.diamond.lines = "darkblue")



#偏倚分析----
metaI <- metagen(
  TE = UMSARS_I_diff,
  lower = UMSARS_I_lci,
  upper = UMSARS_I_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE
)
metaII <- metagen(
  TE = UMSARS_II_diff,
  lower = UMSARS_II_lci,
  upper = UMSARS_II_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE
)
metatotal <- metagen(
  TE = UMSARS_total_diff,
  lower = UMSARS_total_lci,
  upper = UMSARS_total_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE
)
library(meta)
##egger检验
metabias(metaI, k.min = 7, method.bias = "Egger")
metabias(metaI,  method.bias = "Egger", plotit = T)
##Begg检验
metabias(metatotal, method.bias = "Begg", correct =T, plotit = T)#如果样本数少，设置k.min=7
metabias(metaI, method.bias = "linreg", correct = TRUE, plotit = TRUE)

##绘制漏斗图
funnel(metaII)
# 正确的函数名是update()，不是update.meta()
meta_sub <- update(metatotal, exclude = metatotal$studlab == "Fabbri M, 2025")

# 绘制过滤后的漏斗图
funnel(meta_sub)
funnel(meta_sub,studlab = T,pos.studlab = 1)
# meta$studlab <- ifelse(hrdata$study == "O'Sullivan SS, 2008", hrdata$study, "")
 funnel(metaI,studlab = T,pos.studlab = 1) # pos.studlab = 1 更改标注名称位置
#meta回归----
str(data)
meta_regression <- metareg(metaII, year+rigion+type+sample+disease_duration+follow_duration) #单变量meta回归
plot(meta_regression)
summary(meta_regression)
meta_reg_follow_I <- metareg(metaI, ~ disease_duration)
summary(meta_reg_follow_I)

meta_reg_follow_II <- metareg(metaII, ~ disease_duration)
summary(meta_reg_follow_II)

meta_reg_follow_total <- metareg(metatotal, ~ disease_duration)
summary(meta_reg_follow_total)


#功效分析----
library(pwr)
library(metafor)
Idata <- subset(data, !is.na(UMSARS_I_diff))
metaI <- metagen(
  TE = UMSARS_I_diff,
  lower = UMSARS_I_lci,
  upper = UMSARS_I_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE
)
metaII <- metagen(
  TE = UMSARS_II_diff,
  lower = UMSARS_II_lci,
  upper = UMSARS_II_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE
)
metatotal <- metagen(
  TE = UMSARS_total_diff,
  lower = UMSARS_total_lci,
  upper = UMSARS_total_uci,
  data = Idata,
  sm = "MD",
  studlab = study,
  comb.fixed = F,
  comb.random = TRUE
)


#  SD（标准差）----
# 1. 计算平均研究内方差
within_study_variances <- Idata$UMSARS_II_sd^2  # 各研究的方差
avg_within_variance <- mean(within_study_variances, na.rm = TRUE)

# 2. 获取研究间方差（来自随机效应meta分析）
between_study_variance <- metaII$tau2  # 这是τ²

# 3. 计算总方差和总标准差
total_variance <- between_study_variance + avg_within_variance
SD <- sqrt(total_variance)

#合并效应量----
meta_random <- metaII$TE.random     # 合并效应量

# 效应量百分比 20%~100%
effect_percent <- seq(0.2 , 1, by=0.1)  # 20%,25%,30%...100%
mean_diff <- effect_percent * meta_random        # 对应均值差

# 功效水平 80% 和 90%
power_levels <- c(0.8, 0.9)


# 创建空数据框存储结果
result <- data.frame()

# 循环计算样本量
for(pw in power_levels){
  for(i in 1:length(mean_diff)){
    md <- mean_diff[i]
    percent <- effect_percent[i] * 100  # 转换为百分比
    
    # 使用原始均值差和标准差计算样本量
    n <- power.t.test(
      delta = md,
      sd = SD,
      power = pw,
      sig.level = 0.05,
      type = "two.sample",
      alternative = "two.sided"
    )$n
    
    result <- rbind(result, data.frame(
      effect_percent = percent,  # 效应量百分比
      power = pw,
      n_per_group = ceiling(n),
      mean_diff = md,            # 实际均值差
      cohens_d = md / SD         # Cohen's d（仅供参考）
    ))
  }
}

# 查看结果
print(result)
cat("\n合并效应量:", round(meta_random, 2))
cat("\n总标准差:", round(SD, 2))

 ggplot(result, aes(x=effect_percent, y=n_per_group, color=factor(power))) +
  geom_line(size=1) +
  geom_point(size=2) +
  scale_color_manual(values=c("#3861ac","red"), labels=c("80% Power","90% Power")) +
  scale_x_continuous(
    breaks = seq(20, 100, by = 10)  # 从20到100，每隔20一个刻度
  )+
  scale_y_continuous(
    breaks = seq(0, 650, by = 100),  # 设置y轴刻度
    limits = c(0, 650)  # 设置y轴范围
  )+
  labs(
    x = "Effect size (%) of 1-year UMSARS II change",
    y = "Sample size (per group)",
    color = "Statistical power",
    title = ""
  ) +
  theme_prism(base_size = 14)+
  theme(legend.position = c(0.95, 0.95),  # 右上角位置 (x, y 坐标，0-1范围)
                                   legend.justification = c(1, 1))






