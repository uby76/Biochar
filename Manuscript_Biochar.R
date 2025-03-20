#===============================================
#alpha多样性计算
library(vegan)
library(ggpubr)
library(ggplot2)
library(patchwork)
data<-read.csv("asv.csv",row.names = 1) 
otu <- t(data)
shannon<-diversity(otu,index = "shannon")
pielou <- shannon/ log(specnumber(otu), exp(1))
simpson<-diversity(otu,index = "invsimpson")
chao1<-estimateR(otu)[2,]
ACE<-estimateR(otu)[4,]
richness<-estimateR(otu)[1,]
mdata<-rbind(shannon,pielou,simpson,chao1,ACE,richness)
write.csv(mdata,"index_16s.csv")
#===============================================
#混合效应模型的计算
library(dplyr)
library(ggplot2)
library(lme4)
library(boot)
library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(performance)
set.seed(123)  # 设置随机种子

data<-read.csv("alldata_index_16s.csv",row.names = 1)
# data <- data[data$Treatment == "zBiochar", ]
model <- lmer(ratio ~ addnum + (1 | project), data = data, REML = TRUE)
summary(model)  # 查看模型结果

r2_model <- r2(model)
print(r2_model)
#===============================================
#机器学习确定关键的微生物genus
library(caret)
library(randomForest)
library(e1071)
library(pROC)
setwd("")

# 读取数据并设置因子
df <- read.csv("merged_data_new.csv", row.names = 1, check.names = FALSE)  # 以第一列为行名
df$Treatment <- factor(df$Treatment)

# 提取自变量和因变量
otu_data <- (df[, -which(names(df) == "Treatment")])  # 转置数据并去除 Treatment 列
treatment <- df$Treatment  # 提取 Treatment 列
set.seed(365)

# 加载必要的包
library(pROC)
library(ROCR)
library(car)
library(randomForest)
library(e1071)
library(ggplot2)

# 导入数据
mydata <- read.csv("merged_data_new.csv", row.names = 1, check.names = FALSE)  
mydata$Treatment <- factor(mydata$Treatment)

# 检查 Treatment 是否为二分类
if (nlevels(mydata$Treatment) != 2) {
    stop("Error: Treatment 必须为二分类。")
}

# 替换列名中的斜杠
colnames(mydata) <- gsub("/", ".", colnames(mydata))

# 检查数据集中是否存在缺失值
if (any(is.na(mydata))) {
  # 可以使用均值填充缺失值
  mydata <- na.aggregate(mydata)
}

# 划分自变量和因变量
otu_data <- mydata[, -which(names(mydata) == "Treatment")]
treatment <- mydata$Treatment

# 数据标准化 (推荐)
otu_data <- scale(otu_data)

# 合并标准化后的自变量和因变量
mydata_scaled <- data.frame(otu_data, Treatment = treatment)

# 划分训练集和测试集
set.seed(123)  # 设置随机种子以确保结果可复现
ind <- sample(2, nrow(mydata_scaled), replace = TRUE, prob = c(0.8, 0.2))
train <- mydata_scaled[ind == 1, ]  # 训练集
test <- mydata_scaled[ind == 2, ]   # 测试集

# 确保 Treatment 的因子顺序一致
train$Treatment <- factor(train$Treatment, levels = levels(treatment))
test$Treatment <- factor(test$Treatment, levels = levels(treatment))

# 将 Treatment 转为 0/1 数值 (Biochar = 0, Control = 1)
train$Treatment_numeric <- ifelse(train$Treatment == levels(treatment)[1], 0, 1)
test$Treatment_numeric <- ifelse(test$Treatment == levels(treatment)[1], 0, 1)

### -------------------- Logistic 回归模型 -------------------- ###
# 构建 Logistic 回归模型
lr_model <- glm(Treatment_numeric ~ ., data = train, family = binomial)
summary(lr_model)

# 在测试集上进行预测，获取概率
pred_prob_lr <- predict(lr_model, newdata = test, type = "response")

# 确保 pred_prob_lr 是数值类型
pred_prob_lr <- as.numeric(pred_prob_lr)

# 检查 NA 并移除
if (sum(is.na(pred_prob_lr)) > 0) {
    warning("预测概率中包含 NA，已移除相关样本。")
    test <- test[!is.na(pred_prob_lr), ]
    pred_prob_lr <- pred_prob_lr[!is.na(pred_prob_lr)]
}

# 转换为预测类别（0.5 为阈值）
pred_class_lr <- ifelse(pred_prob_lr > 0.5, 1, 0)
pred_class_lr <- factor(pred_class_lr, levels = c(0, 1))
test$Treatment_numeric <- factor(test$Treatment_numeric, levels = c(0, 1))

# 绘制混淆矩阵
confusion_matrix_lr <- table(Predicted = pred_class_lr, Actual = test$Treatment_numeric)
print(confusion_matrix_lr)

# 计算 ROC 和 AUC
pred_lr <- prediction(pred_prob_lr, as.numeric(test$Treatment_numeric))
auc_lr <- performance(pred_lr, "auc")
auc_value_lr <- auc_lr@y.values[[1]]
print(paste("Logistic 回归 AUC:", auc_value_lr))

### -------------------- 随机森林模型 -------------------- ###
# 确保训练集和测试集的因子顺序一致
train$Treatment <- factor(train$Treatment, levels = levels(treatment))
test$Treatment <- factor(test$Treatment, levels = levels(treatment))

# 移除训练集和测试集中的 Treatment_numeric 列
train_no_numeric <- train[, -which(names(train) == "Treatment_numeric")]
test_no_numeric <- test[, -which(names(test) == "Treatment_numeric")]

# 构建随机森林模型的公式
formula_rf <- as.formula("Treatment ~ .")

# 构建随机森林模型
rf_model <- randomForest(formula_rf, data = train_no_numeric, importance = TRUE, proximity = TRUE, ntree = 500)
print(rf_model)

# 在测试集上进行预测，获取概率
pred_prob_rf <- predict(rf_model, newdata = test_no_numeric, type = "prob")[, 2]

# 确保 pred_prob_rf 是数值类型
pred_prob_rf <- as.numeric(pred_prob_rf)

# 计算 ROC 和 AUC
pred_rf <- prediction(pred_prob_rf, as.numeric(test$Treatment_numeric))
auc_rf <- performance(pred_rf, "auc")
auc_value_rf <- auc_rf@y.values[[1]]
print(paste("随机森林 AUC:", auc_value_rf))

### -------------------- SVM 模型 -------------------- ###
# 构建 SVM 模型 (径向基核函数)
svm_model <- svm(Treatment ~ ., data = train_no_numeric, type = "C-classification", kernel = "radial", probability = TRUE)
summary(svm_model)

# 在测试集上进行预测，获取概率
pred_prob_svm <- attr(predict(svm_model, newdata = test_no_numeric, probability = TRUE), "probabilities")[, 2]

# 确保 pred_prob_svm 是数值类型
pred_prob_svm <- as.numeric(pred_prob_svm)

# 计算 ROC 和 AUC
pred_svm <- prediction(pred_prob_svm, as.numeric(test$Treatment_numeric))
auc_svm <- performance(pred_svm, "auc")
auc_value_svm <- auc_svm@y.values[[1]]
print(paste("SVM AUC:", auc_value_svm))

### -------------------- 绘制 ROC 曲线 -------------------- ###
# 计算 ROC 曲线
roc_perf_lr <- performance(pred_lr, "tpr", "fpr")
roc_perf_rf <- performance(pred_rf, "tpr", "fpr")
roc_perf_svm <- performance(pred_svm, "tpr", "fpr")

# 绘制 ROC 曲线
plot(roc_perf_lr, main = "ROC Curve Comparison", col = "blue", lwd = 2, lty = 1)
plot(roc_perf_rf, col = "green", lwd = 2, lty = 2, add = TRUE)
plot(roc_perf_svm, col = "red", lwd = 2, lty = 3, add = TRUE)
abline(a = 0, b = 1, lwd = 2, lty = 2, col = "gray")

# 添加图例和 AUC 值
legend("bottomright", legend = c(
    paste("Logistic (AUC =", round(auc_value_lr, 3), ")"),
    paste("Random Forest (AUC =", round(auc_value_rf, 3), ")"),
    paste("SVM (AUC =", round(auc_value_svm, 3), ")")),
    col = c("blue", "green", "red"),
    lwd = 2, lty = c(1, 2, 3)
)

#============================================================
#============================================================
#============================================================

# 进行 5 次交叉验证
error.cv <- c()
for (i in 1:5) {
  print(i)
  set.seed(i)
  
  # 使用 rfcv 函数进行交叉验证
  fit <- rfcv(trainx = train_no_numeric[, -which(names(train_no_numeric) == "Treatment")],
              trainy = train_no_numeric$Treatment,
              cv.fold = 5,
              scale = "log",
              step = 0.9,
              ntree = 500)  # 使用之前随机森林模型的 ntree 值
  
  # 存储每次交叉验证的误差
  error.cv <- cbind(error.cv, fit$error.cv)
}

# 获取变量数量
n.var <- as.numeric(rownames(error.cv))

# 为误差矩阵的列命名
colnames(error.cv) <- paste('error', 1:5, sep = '.')

# 计算平均误差
err.mean <- apply(error.cv, 1, mean)

# 创建误差数据框
err.df <- data.frame(num = n.var,
                     err.mean = err.mean,
                     error.cv)

# 查看误差数据框的前几行
head(err.df[, 1:6])

# 打印误差数据框
print(err.df)

# 保存误差数据框到 CSV 文件
write.csv(err.df, "err.df.csv", row.names = TRUE)

err.df <- read.csv("err.df.csv", row.names = 1, check.names = FALSE)  # 以第一列为行名

optimal <- err.df$num[which(err.df$err.mean == min(err.df$err.mean))]
main_theme <-
 theme(
   panel.background = element_blank(),
   panel.grid = element_blank(),
   axis.line.x = element_line(linewidth = 0.5, color = "black"),
   axis.line.y = element_line(linewidth = 0.5, color = "black"),
   axis.ticks = element_line(color = "black"),
   axis.text = element_text(color = "black", size = 12),
   legend.position = "right",
   legend.background = element_blank(),
   legend.key = element_blank(),
   legend.text = element_text(size = 12),
   text = element_text(family = "sans", size = 12))

pl <-
 ggplot(data = err.df, aes(x = err.df$num)) +
   geom_line(aes(y = err.df$error.1), color = 'grey', linewidth = 0.5) +
   geom_line(aes(y = err.df$error.2), color = 'grey', linewidth = 0.5) +
   geom_line(aes(y = err.df$error.3), color = 'grey', linewidth = 0.5) +
   geom_line(aes(y = err.df$error.4), color = 'grey', linewidth = 0.5) +
   geom_line(aes(y = err.df$error.5), color = 'grey', linewidth = 0.5) +
   geom_line(aes(y = err.df$err.mean), color = 'black', linewidth = 0.5) +
   geom_vline(xintercept = optimal, color = 'red', lwd = 0.36, linetype = 2) +
   coord_trans(x = "log2") +
   scale_x_continuous(breaks = c(1, 10, 30, 50,80)) +
   labs(x = 'Number of Features ', y = 'Cross-validation error rate') +
   annotate("text",
            x = optimal,
            y = max(err.df$err.mean),
            label = paste("Optimal = ", optimal, sep = ""),
            color = "red") +
   main_theme
pl

ggsave('line_forest.pdf', pl, width = 5, height = 3.5)

# 获取特征重要性
importance <- randomForest::importance(rf_model)

# 转换为数据框并排序
importance_df <- data.frame(Feature = rownames(importance), Importance = importance[, "MeanDecreaseAccuracy"])
importance_df <- importance_df[order(-importance_df$Importance), ]

# 提取前138个特征
top_138_features <- head(importance_df, 138)

# 导出到 CSV 文件
write.csv(top_138_features, "top_138_features_forest.csv", row.names = FALSE)
#======================================================
#构建相关性网络
library(Hmisc)
 
#以属水平丰度为例，“genus_table.txt” 是一个属水平的微生物丰度表
genus <- read.csv("Biochar_asv_end_re_all.csv", row.names = 1, check.names = FALSE)  # 以第一列为行名

# genus1 <- genus
# genus1[genus1>0] <- 1
# genus <- genus[which(rowSums(genus1) >= 260), ]    #例如只保留在 5 个及以上样本中出现的属
# biochar:260 control:115

#计算两属之间是否存在丰度变化的相关性，以 spearman 相关系数为例
genus_corr <- rcorr(t(genus), type = 'spearman')
 
#阈值筛选
#将 spearman 相关系数低于 0.7 的关系剔除，即 r>=0.7
r <- genus_corr$r
r[abs(r) < 0.4] <- 0
 
#选取显著性 p 值小于 0.05 的相关系数，即 p<0.05
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
 
#根据上述筛选的 r 值和 p 值保留数据
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
 
#如此便得到了邻接矩阵格式的网络文件（微生物属的相关系数矩阵）
write.table(data.frame(z, check.names = FALSE), 'Biochar_genus_corr.matrix.txt', col.names = NA, sep = '\t', quote = FALSE)

##获得网络
library(igraph)
 
#将邻接矩阵转化为 igraph 网络的邻接列表
#构建含权的无向网络，权重代表了微生物属间丰度的 spearman 相关系数
g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected')
g
 
#自相关也可以通过该式去除
g <- simplify(g)
 
#孤立节点的删除（删除度为 0 的节点）
g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
 
#该模式下，边权重代表了相关系数
#由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
E(g)$correlation <- E(g)$weight
E(g)$weight <- abs(E(g)$weight)

plot(g)

#边列表
edge <- data.frame(as_edgelist(g))    #igraph 的邻接列表转为边列表
 
edge_list <- data.frame(
    source = edge[[1]],
    target = edge[[2]],
    weight = E(g)$weight,
    correlation = E(g)$correlation
)
head(edge_list)

#节点属性列表
node_list <- data.frame(
    label = names(V(g))
)
head(node_list)
write.table(edge_list, 'Biocharnetwork.edge_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)
write.table(node_list, 'Biocharnetwork.node_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)




######计算网络常用的几种拓扑系数#####
nodes_num = length(V(g))#节点数
edges_num = length(E(g))#边数
average_degree = mean(degree(g))#平均度
positive.cor_num = sum(E(g)$correlation>0)        #正相关的数量
negative.cor_num = sum(E(g)$correlation<0)        #负相关的数量



average_path_length = average.path.length(g, directed = FALSE)#平均路径长度
network_diameter = diameter(g, directed = FALSE)#网络直径
network_density = graph.density(g)#网络密度
clustering_coefficient = transitivity(g)#聚类系数


network_parameter = data.frame(nodes_num, 
                               edges_num,
                               positive.cor_num,
                               negative.cor_num,
                               average_degree,
                               average_path_length,
                               network_diameter, 
                               network_density,
                               clustering_coefficient                               
)

network_parameter
write.csv(network_parameter, 'network_parameter_Biochar.csv')
#======================================================
#随机抽取特定数量后构建相关性网络
library(Hmisc)
library(igraph)

# 读取数据
genus <- read.csv("Control_asv_end_re_all.csv", row.names = 1, check.names = FALSE)

# 设定重复次数和每次的样本数量
num_iterations <- 100
sample_size <- 300

# 创建一个列表用于存储每次计算的网络参数
network_results <- list()

set.seed(123)  # 确保随机性可复现

for (i in 1:num_iterations) {
    cat("Iteration:", i, "\n")
    
    # **随机抽取 300 个样本**
    sampled_indices <- sample(1:ncol(genus), sample_size, replace = FALSE)
    genus_sampled <- genus[, sampled_indices]
    
    # **计算 Spearman 相关系数**
    genus_corr <- rcorr(t(genus_sampled), type = 'spearman')
    
    # **筛选相关性**
    r <- genus_corr$r
    r[abs(r) < 0.4] <- 0  # 仅保留绝对值大于 0.4 的相关系数
    
    p <- genus_corr$P
    p <- p.adjust(p, method = 'BH')  # BH 方法校正 P 值
    p[p >= 0.05] <- -1
    p[p < 0.05 & p >= 0] <- 1
    p[p == -1] <- 0
    
    z <- r * p  # 只保留显著相关的相关系数
    diag(z) <- 0  # 去除自相关
    
    # **确保对称性**
    z <- (z + t(z)) / 2  # Make sure the matrix is symmetric
    z[is.na(z)] <- 0      # Replace NA values with 0 (in case there are any after p-value adjustments)
    
    
    # **构建网络**
    g <- graph_from_adjacency_matrix(z, weighted = TRUE, mode = 'undirected')
    g <- simplify(g)  # 去除自环
    g <- delete_vertices(g, names(degree(g)[degree(g) == 0]))  # 删除孤立节点
    
    # 计算边的属性
    E(g)$correlation <- E(g)$weight
    E(g)$weight <- abs(E(g)$weight)
    
    # **计算网络参数**
    nodes_num = length(V(g))  # 节点数
    edges_num = length(E(g))  # 边数
    average_degree = mean(degree(g))  # 平均度
    positive.cor_num = sum(E(g)$correlation > 0)  # 正相关数量
    negative.cor_num = sum(E(g)$correlation < 0)  # 负相关数量
    average_path_length = mean_distance(g, directed = FALSE)  # 平均路径长度
    network_diameter = diameter(g, directed = FALSE)  # 网络直径
    network_density = edge_density(g)  # 网络密度
    clustering_coefficient = transitivity(g)  # 聚类系数
    
    # **存储本次计算结果**
    network_results[[i]] <- data.frame(
        iteration = i,
        nodes_num = nodes_num,
        edges_num = edges_num,
        positive.cor_num = positive.cor_num,
        negative.cor_num = negative.cor_num,
        average_degree = average_degree,
        average_path_length = average_path_length,
        network_diameter = network_diameter,
        network_density = network_density,
        clustering_coefficient = clustering_coefficient
    )
}

# **合并所有计算结果并导出**
final_network_results <- do.call(rbind, network_results)
write.csv(final_network_results, "network_parameters_100_iterations.csv", row.names = FALSE)
cat("100 次随机抽样分析已完成，结果已保存至 'network_parameters_100_iterations.csv'")

#======================================================
#计算环境变量对多样性指数的影响
library(randomForest)
setwd("/forest/网络/subnetwork")
# 加载必要的包
library(randomForest)
library(caret)

# 读取数据
data_forest <- read.csv("forest.csv", row.names = 1, check.names = FALSE)

# 删除 shannon 和 avedtem 列中有 NA 的行
# data_forest <- data_forest[!is.na(data_forest$ratio) & !is.na(data_forest$avedtem), ]

# 提取目标变量和特征变量
response_var <- data_forest$AVD  # 响应变量
explanatory_vars <- data_forest[, 3:14]  # 解释变量

# 创建新的数据框
new_data <- data.frame(response_var = response_var, explanatory_vars)

# 删除含有缺失值的行
new_data_complete <- na.omit(new_data)

# 设置随机种子以保证结果可重复
set.seed(123)

# 调优 mtry 参数
optimal_mtry <- tuneRF(
  new_data_complete[, -1],  # 去除目标变量列
  new_data_complete$response_var,  # 响应变量
  ntreeTry = 500,            # 每次尝试的树数
  stepFactor = 1.5,          # mtry 调整步长因子
  improve = 0.01,            # 最小提升比例
  trace = TRUE,              # 打印调优过程
  plot = TRUE                # 显示调优图
)

# 获取最佳 mtry 值
best_mtry <- optimal_mtry[which.min(optimal_mtry[, 2]), 1]  # 最优 mtry 对应的值

# 使用最佳 mtry 和优化节点数 (nodesize) 训练模型
forest_model <- randomForest(
  response_var ~ ., 
  data = new_data_complete, 
  importance = TRUE, 
  ntree = 500,  # 增加树的数量以提高稳定性
  mtry = best_mtry, 
  nodesize = 5   # 优化叶子节点的最小样本数
)

# 查看变量重要性
importance_scores <- importance(forest_model, scale = TRUE)
varImpPlot(forest_model)  # 可视化变量重要性

# 根据重要性筛选特征
important_vars <- rownames(importance_scores[importance_scores[, 1] > 0, ])
new_data_reduced <- new_data_complete[, c("response_var", important_vars)]

# 用筛选后的重要特征重新训练随机森林模型
forest_model_optimized <- randomForest(
  response_var ~ ., 
  data = new_data_reduced, 
  importance = TRUE, 
  ntree = 500, 
  mtry = best_mtry, 
  nodesize = 5
)

# 输出最终模型
print(forest_model_optimized)

# 使用交叉验证评估模型性能
train_control <- trainControl(method = "cv", number = 5)  # 5 折交叉验证
cv_model <- train(
  response_var ~ ., 
  data = new_data_reduced, 
  method = "rf", 
  trControl = train_control, 
  ntree = 500
)
print(cv_model)

# 保存模型重要性得分为文件
write.csv(importance_scores, "variable_importance.csv")


# 加载必要的包
library(ggplot2)

# 从最终优化的随机森林模型中提取变量重要性
importance_scores <- importance(forest_model_optimized, scale = TRUE)
importance_df <- data.frame(
  Variable = rownames(importance_scores),
  X.IncMSE = importance_scores[, 1]  # 使用 %IncMSE 列
)

# 对变量按重要性排序

importance_df <- importance_df[order(importance_df$X.IncMSE, decreasing = TRUE), ]
importance_df$Variable <- factor(importance_df$Variable, levels = importance_df$Variable)

# 绘制变量重要性图
resp_imp <- ggplot(importance_df, aes(x = X.IncMSE, y = Variable)) + 
  geom_segment(aes(yend = Variable, xend = 0), linetype = 2) +
  geom_point(shape = 21, size = 5, fill = "grey") +
  ggtitle(bquote('Ratio:'~R^2~'= 0.64')) +  # 替换为实际的 R^2 值
  theme_classic() +
  ylab("") +
  theme(axis.title = element_text(size = 9)) +
  theme(text = element_text(size = 11)) +
  theme(legend.text = element_text(color = "black", size = 15)) +
  xlab("Variable Importance (% MSE Increase)") +
  theme(axis.line.x = element_line(colour = "black", size = 1)) + 
  theme(axis.line.y = element_line(colour = "black", size = 1)) + 
  theme(legend.title = element_blank()) +
  theme(legend.position = c(.8, .15)) +
  theme(legend.background = element_blank(), legend.box.background = element_rect(colour = "black"))+
  coord_flip()
# 绘制图表
plot(resp_imp)

ggsave('AVD_forest.pdf', resp_imp, width = 5, height = 3)
#======================================================
#群落构建
setwd("/betanti")

#请确保安装了这些 R 包
library(Hmisc)
library(minpack.lm)
library(stats4)

# 设置工作目录到文件所在位置

#spp: 物种或分类群的丰度表，行是分类群，列是样本
spp <- read.csv("Control_asv_end_re_all.csv", row.names = 1, check.names = FALSE)
spp<-t(spp)
 
##将 Sloan 等（2006）的中性模型拟合到一个物种或分类群的丰度表，并返回几个拟合统计数据。或者，将根据它们在元群落中的丰度返回每个分类群的预测出现频率
#用非线性最小二乘法（Non-linear least squares，NLS）拟合模型参数
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit  #获取 m 值
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr  #获取模型的 R2
 
#输出 3 个统计结果数据表，包括各物种或分类群的平均相对丰度（p.csv）、出现频率（freq.csv）和预测的出现频率（freq.pred.csv）
write.csv(p, file = "all_asv_p.csv")
write.csv(freq, file = "all_asv_freq.csv")
write.csv(freq.pred, file = "all_asv_freq.pred.csv")
 
#p 是平均相对丰度（mean relative abundance）
#freq 是出现频率（occurrence frequency）的观测值
#freq.pred 是出现频率（occurrence frequency）的预测值，即中性模型的拟合值
 
#绘制统计图
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#出现频率低于中性群落模型预测的部分
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#出现频率高于中性群落模型预测的部分
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')
 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

#6:5的比例输出

setwd("/betanti/NST")
#install.packages("NST")
library(NST)
library(dplyr)

 

 
comm <- read.csv("comm.csv", row.names = 1, check.names = FALSE)

 
group <- read.csv("group.csv", row.names = 1, check.names = FALSE) 
group


set.seed(123)
tnst <- tNST(comm = comm, group = group, dist.method = 'jaccard', null.model = 'PF', 
    rand = 1000, nworker = 1)
 
nst_group <- tnst$index.pair.grp
nst_group
 
#输出主要的统计结果
write.table(nst_group, 'nst_group.txt', sep = '\t', row.names = FALSE, quote = FALSE)



#为了更好地评估确定性和随机过程在两组群落中的相对重要性
#不妨简单绘制个箱线图，将两组群落内 MST 的分布可视化
#并通过统计检验（如这里以非参数的 wilcox test 为例），比较两组群落内 MST 数值是否存在显著差异
library(ggpubr)
 
ggboxplot(data = nst_group, x = 'group', y = 'MST.ij.ruzicka', color = 'group') +
stat_compare_means(method = 'wilcox.test', comparisons = list(c('Control', 'Treat'))) +
labs(y = 'Modified Stochasticity Ratio (MST)')


library(tidyverse)
library(gghalves)

ordercolors<-c( "#539ed5","#fdc187")

pend = ggplot(data = nst_group,
       aes(x=group, y=MST.ij.ruzicka, fill=group)) +
    geom_half_violin(side = "r", color=NA, alpha=0.35) +
    geom_half_boxplot(side = "r", errorbar.draw = FALSE, width=0.2, linewidth=0.5) +
    scale_fill_manual(values = ordercolors) +
    scale_x_discrete(labels = c('Biochar','Control')) +
    labs(y="Brain Proportion (%)",x=NULL) +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 16, color = "black"),
          axis.text = element_text(size=13, color = "black"))

ggsave('pend.pdf', pend, width = 6, height = 5)

#======================================================


