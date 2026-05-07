# 1. 加载必要的包
#if (!require("ropls")) install.packages("ropls")  # OPLS-DA核心包
#if (!require("ggplot2")) install.packages("ggplot2")  # 绘图
#if (!require("dplyr")) install.packages("dplyr")  # 数据处理
#if (!require("pheatmap")) install.packages("pheatmap")  # 热图
#if (!require("caret")) install.packages("caret")  # 交叉验证
# 加载必要的包
library(ropls)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(caret)

# 2. 数据准备与检查
# 假设您的数据框名为 all_mut_robust
# 检查数据结构
cat("数据维度:", dim(all_mut_robust), "\n")
cat("疾病类型:", unique(all_mut_robust$Disease), "\n")
cat("样本数量:", nrow(all_mut_robust), "\n")

# 分离特征矩阵X和响应变量Y
# 特征矩阵X: 第5到第19列（Co.evolution到Stiffness）
X <- as.matrix(all_mut_robust[, 5:19])
rownames(X) <- paste(all_mut_robust$Gene, all_mut_robust$Site, all_mut_robust$Mutation, sep="_")
# 【插入位置：设置行名之后，检查摘要之前】
# 1. 备份行列名
row_names_backup <- rownames(X)
col_names_backup <- colnames(X)

# 2. 将矩阵转换为数据框，便于逐列操作
X_df <- as.data.frame(X, stringsAsFactors = FALSE)

# 3. 遍历每一列，执行转换
X_numeric <- X_df # 创建副本用于存储结果
conversion_issues <- list() # 记录转换问题

for (col in colnames(X_df)) {
  # 3.1 将科学计数法中的大写E统一替换为小写e (R标准)
  temp_col <- gsub("E", "e", X_df[[col]], ignore.case = FALSE) # 仅替换大写E
  
  # 3.2 尝试转换为数值型
  converted <- suppressWarnings(as.numeric(temp_col))
  
  # 3.3 检查转换是否成功，记录问题
  na_count <- sum(is.na(converted))
  original_na_count <- sum(is.na(X_df[[col]]))
  
  if (na_count > original_na_count) {
    # 找出新增的NA（即转换失败的值）
    new_nas <- which(is.na(converted) & !is.na(X_df[[col]]))
    if (length(new_nas) > 0) {
      conversion_issues[[col]] <- list(
        count = length(new_nas),
        examples = head(X_df[[col]][new_nas], 3) # 记录前3个问题值
      )
    }
  }
  
  # 3.4 存储转换后的列
  X_numeric[[col]] <- converted
}

# 4. 报告转换问题（如果有）
if (length(conversion_issues) > 0) {
  cat("\n⚠️ 警告：以下列的部分值在转换时变成了NA（可能包含非数字字符）：\n")
  for (col_name in names(conversion_issues)) {
    cat(sprintf("  列 '%s': %d 个值转换失败，示例: %s\n", 
                col_name, 
                conversion_issues[[col_name]]$count,
                paste(conversion_issues[[col_name]]$examples, collapse = ", ")))
  }
  cat("\n")
} else {
  cat("\n✅ 所有列已成功转换为数值型。\n")
}

# 5. 将清理后的数据框转换回矩阵，并恢复行列名
X <- as.matrix(X_numeric)
rownames(X) <- row_names_backup
colnames(X) <- col_names_backup

# 6. 验证转换结果
cat("转换后数据类型:", class(X), "\n")
cat("转换后数据模式:", mode(X), "\n")
cat("检查是否还有非数值:", if(any(!is.numeric(X))) "是" else "否", "\n")
# 响应变量Y: 疾病类型（第二列）
Y <- factor(all_mut_robust$Disease)
names(Y) <- rownames(X)

# 检查特征矩阵的统计摘要
cat("\n特征矩阵统计摘要:\n")
print(summary(X))

# 检查是否有缺失值
if (any(is.na(X))) {
  cat("警告: 特征矩阵中存在缺失值!\n")
  # 可以选择删除有缺失值的行或进行插补
  na_rows <- which(apply(X, 1, function(x) any(is.na(x))))
  cat("含有缺失值的行:", na_rows, "\n")
}

# 3. 确定PLS-DA的预测成分数
# 对于多分类问题，预测成分数通常为min(变量数, 类别数-1, 样本数-1)
pred_components <- min(ncol(X), nlevels(Y)-1, nrow(X)-1)
cat("\n>>> 确定的预测成分数:", pred_components, "\n")

# 4. PLS-DA模型构建
cat("\n>>> 开始构建PLS-DA模型...\n")

# 使用ropls包进行PLS-DA分析（设置orthoI=0表示没有正交成分）
plsda_model <- opls(X, Y, 
                    predI = pred_components,  # 预测成分数
                    orthoI = 0,  # 正交成分数为0，表示这是PLS-DA而非OPLS-DA
                    crossvalI = min(7, nrow(X)),  # 交叉验证折数
                    log10L = FALSE,  # 不对数据取对数
                    scaleC = "none",  # 数据已经标准化，这里不再缩放
                    fig.pdfC = "none",  # 不输出PDF文件
                    info.txtC = "none")  # 不输出文本信息

# 5. 模型结果概览
cat("\n>>> PLS-DA模型结果概览:\n")
print(plsda_model)

# 检查模型质量指标
cat("\n模型质量指标:\n")
cat("R2X(cum):", plsda_model@modelDF$R2X[nrow(plsda_model@modelDF)], "\n")
cat("R2Y(cum):", plsda_model@modelDF$R2Y[nrow(plsda_model@modelDF)], "\n")
cat("Q2(cum):", plsda_model@modelDF$Q2[nrow(plsda_model@modelDF)], "\n")

# Q2是交叉验证的预测能力指标，通常>0.5表示模型良好
if (plsda_model@modelDF$Q2[nrow(plsda_model@modelDF)] < 0.5) {
  cat("警告: Q2值较低，模型预测能力可能不足。\n")
}

# 6. 模型诊断可视化
cat("\n>>> 生成模型诊断图...\n")

# 6.1 模型概述图
plot(plsda_model, 
     typeVc = "overview")

# 6.2 预测得分图（t1 vs t2）
plot(plsda_model, 
     typeVc = "x-score")

# 6.3 载荷图（p1 vs p2）
plot(plsda_model, 
     typeVc = "x-loading"
     )

# 6.4 重要特征分析
# 提取VIP值
if (!is.null(getVipVn(plsda_model))) {
  vip_values <- getVipVn(plsda_model)
  
  # 创建特征重要性数据框
  feature_importance <- data.frame(
    Feature = colnames(X),
    VIP = vip_values
  )
  
  # 按VIP值排序
  feature_importance <- feature_importance[order(-feature_importance$VIP), ]
  
  cat("\n>>> 特征重要性排序(VIP):\n")
  print(feature_importance)
  
  # 可视化VIP值
  vip_plot <- ggplot(feature_importance, aes(x = reorder(Feature, VIP), y = VIP)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
    coord_flip() +
    labs(x = "特征", y = "VIP值", 
         title = "PLS-DA特征重要性(VIP)",
         subtitle = paste("红色虚线表示VIP=1的阈值")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  
  print(vip_plot)
  
  # 6.5 载荷分析与权重提取
  # 对于多分类PLS-DA，通常关注第一个预测成分的载荷
  if (!is.null(plsda_model@loadingMN)) {
    # 获取第一个预测成分的载荷
    p1_values <- plsda_model@loadingMN[, 1]
    
    # 创建特征权重数据框
    feature_importance$Loading_p1 <- p1_values[match(feature_importance$Feature, names(p1_values))]
    feature_importance$abs_Loading <- abs(feature_importance$Loading_p1)
    
    cat("\n>>> 特征载荷排序(按绝对值):\n")
    print(feature_importance[order(-feature_importance$abs_Loading), ])
  }
}

# 7. 提取用于合成总指标的权重
cat("\n>>> 提取特征权重用于合成总指标...\n")

# 对于多分类PLS-DA，有几种权重提取策略

# 方法1：使用第一个预测成分的载荷作为权重
if (!is.null(plsda_model@loadingMN)) {
  weights_p1 <- plsda_model@loadingMN[, 1]
  cat("\n方法1: 基于第一个预测成分载荷的权重\n")
  print(weights_p1)
}

# 方法2：使用VIP值加权的第一个预测成分载荷
if (!is.null(getVipVn(plsda_model)) && !is.null(plsda_model@loadingMN)) {
  weights_vip_weighted <- weights_p1 * (vip_values / mean(vip_values))
  cat("\n方法2: VIP加权的预测载荷权重\n")
  print(weights_vip_weighted)
}

# 方法3：使用多个预测成分的加权组合
# 根据每个预测成分解释的Y方差(R2Y)进行加权
if (!is.null(plsda_model@loadingMN) && ncol(plsda_model@loadingMN) > 1) {
  # 提取各成分的R2Y贡献
  r2y_contrib <- c()
  for (i in 1:ncol(plsda_model@loadingMN)) {
    if (i == 1) {
      r2y_contrib[i] <- plsda_model@modelDF$R2Y[i]
    } else {
      r2y_contrib[i] <- plsda_model@modelDF$R2Y[i] - plsda_model@modelDF$R2Y[i-1]
    }
  }
  
  # 计算加权平均载荷
  weighted_loadings <- rep(0, nrow(plsda_model@loadingMN))
  for (i in 1:length(r2y_contrib)) {
    weighted_loadings <- weighted_loadings + plsda_model@loadingMN[, i] * r2y_contrib[i]
  }
  weighted_loadings <- weighted_loadings / sum(r2y_contrib)
  
  cat("\n方法3: 基于多个预测成分R2Y加权的权重\n")
  print(weighted_loadings)
}

# 创建权重比较数据框（使用方法2，VIP加权）
weights_comparison <- data.frame(
  Feature = colnames(X),
  Loading_p1 = weights_p1,
  VIP = vip_values,
  Weight_VIP_weighted = weights_vip_weighted,
  VIP_importance = ifelse(vip_values > 1, "重要", "一般")
)

cat("\n>>> 权重比较(按VIP加权权重排序):\n")
weights_comparison <- weights_comparison[order(-abs(weights_comparison$Weight_VIP_weighted)), ]
print(weights_comparison)

# 8. 计算样本的综合得分
cat("\n>>> 计算每个样本的综合得分...\n")

# 使用VIP加权的权重计算综合得分
composite_score <- X %*% weights_vip_weighted

# 创建包含综合得分的数据框
results_df <- data.frame(
  Sample = rownames(X),
  Gene = all_mut_robust$Gene,
  Disease = all_mut_robust$Disease,
  Site = all_mut_robust$Site,
  Mutation = all_mut_robust$Mutation,
  Composite_Score = as.numeric(composite_score)
)

# 按疾病类型和综合得分排序
results_df <- results_df[order(results_df$Disease, -results_df$Composite_Score), ]

cat("\n>>> 前10个样本的综合得分:\n")
print(head(results_df, 10))

# 9. 综合得分的可视化
# 按疾病类型分组的箱线图
score_boxplot <- ggplot(results_df, aes(x = Disease, y = Composite_Score, fill = Disease)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  labs(x = "疾病类型", y = "综合得分", 
       title = "不同疾病类型的综合得分分布",
       subtitle = "综合得分基于PLS-DA的VIP加权权重") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")

print(score_boxplot)

# 综合得分热图（前20个样本）
top_samples <- head(results_df[order(-abs(results_df$Composite_Score)), ], 20)

# 提取这些样本的原始特征值
top_X <- X[rownames(X) %in% top_samples$Sample, ]

# 创建注释信息
annotation_row <- data.frame(
  Disease = top_samples$Disease[match(rownames(top_X), top_samples$Sample)],
  Composite_Score = top_samples$Composite_Score[match(rownames(top_X), top_samples$Sample)]
)
rownames(annotation_row) <- rownames(top_X)

# 特征按权重排序
feature_order <- colnames(top_X)[order(-abs(weights_vip_weighted))]

pheatmap(top_X[, feature_order],
         scale = "none",  # 数据已经标准化
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "高综合得分样本的特征热图 (PLS-DA)",
         fontsize_row = 8,
         fontsize_col = 9)

# 10. 模型验证（permutation test）
cat("\n>>> 进行置换检验验证模型...\n")

set.seed(123)  # 设置随机种子以保证可重复性
permutation_results <- opls(X, Y,
                            predI = pred_components,
                            orthoI = 0,  # PLS-DA没有正交成分
                            permI = min(100, nrow(X)*2),  # 置换次数不超过样本数的2倍
                            crossvalI = min(7, nrow(X)),
                            scaleC = "none",
                            info.txtC = "none",
                            fig.pdfC = "none")

# 查看置换检验结果
if (!is.null(permutation_results@suppLs$permMN)) {
  perm_matrix <- permutation_results@suppLs$permMN
  cat("\n置换检验结果:\n")
  print(perm_matrix)
  
  # 提取原始模型的R2Y和Q2
  original_R2Y <- plsda_model@modelDF$R2Y[nrow(plsda_model@modelDF)]
  original_Q2 <- plsda_model@modelDF$Q2[nrow(plsda_model@modelDF)]
  
  # 计算置换检验的p值
  perm_R2Y <- perm_matrix[, "R2Y(cum)"]
  perm_Q2 <- perm_matrix[, "Q2(cum)"]
  
  p_R2Y <- sum(perm_R2Y >= original_R2Y) / length(perm_R2Y)
  p_Q2 <- sum(perm_Q2 >= original_Q2) / length(perm_Q2)
  
  cat("\n置换检验p值:\n")
  cat(sprintf("R2Y的p值: %.4f\n", p_R2Y))
  cat(sprintf("Q2的p值: %.4f\n", p_Q2))
  
  # 可视化置换检验结果
  perm_df <- data.frame(
    Iteration = 1:nrow(perm_matrix),
    R2Y = perm_R2Y,
    Q2 = perm_Q2
  )
  
  perm_plot <- ggplot(perm_df, aes(x = R2Y)) +
    geom_histogram(bins = 20, fill = "lightblue", alpha = 0.7) +
    geom_vline(xintercept = original_R2Y, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = original_R2Y, y = max(hist(perm_R2Y, plot = FALSE)$counts)/2,
             label = paste("原始R2Y =", round(original_R2Y, 3)),
             hjust = -0.1, color = "red") +
    labs(x = "置换检验的R2Y值", y = "频率",
         title = "置换检验: R2Y分布 (PLS-DA)",
         subtitle = paste("p值 =", p_R2Y)) +
    theme_minimal()
  
  print(perm_plot)
}

# 11. 计算每个类别与总体均值的偏最小二乘判别得分
# 这对于理解每个特征如何区分特定疾病类型很有帮助
cat("\n>>> 计算每个疾病类别的特征贡献...\n")

# 获取预测得分
if (!is.null(plsda_model@scoreMN)) {
  # 前两个预测成分的得分
  scores <- plsda_model@scoreMN[, 1:min(2, ncol(plsda_model@scoreMN))]
  
  # 计算每个疾病类别在预测成分上的平均得分
  mean_scores <- aggregate(scores, by = list(Disease = Y), FUN = mean)
  cat("\n各疾病类别的平均预测得分:\n")
  print(mean_scores)
  
  # 可视化预测得分
  scores_df <- data.frame(
    Sample = rownames(scores),
    Disease = Y,
    t1 = scores[, 1],
    t2 = if(ncol(scores) > 1) scores[, 2] else NA
  )
  
  if (ncol(scores) > 1) {
    score_scatter <- ggplot(scores_df, aes(x = t1, y = t2, color = Disease, shape = Disease)) +
      geom_point(size = 3, alpha = 0.7) +
      stat_ellipse(level = 0.95, alpha = 0.2) +  # 添加95%置信椭圆
      labs(x = "第一预测成分 (t1)", y = "第二预测成分 (t2)",
           title = "PLS-DA预测得分图",
           subtitle = "不同疾病类型的样本分布") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(score_scatter)
  }
}

# 12. 保存结果
cat("\n>>> 保存分析结果...\n")

# 保存权重到CSV文件
write.csv(feature_importance, "PLSDA_feature_importance.csv", row.names = FALSE)
write.csv(weights_comparison, "PLSDA_weights_comparison.csv", row.names = FALSE)
write.csv(results_df, "PLSDA_composite_scores.csv", row.names = FALSE)

# 保存重要结果到R数据文件
save(plsda_model, feature_importance, weights_comparison, results_df,
     file = "PLSDA_analysis_results.RData")

cat("\n>>> PLS-DA分析完成!\n")
cat("输出文件:\n")
cat("1. PLSDA_feature_importance.csv - 特征重要性排序\n")
cat("2. PLSDA_weights_comparison.csv - 权重比较\n")
cat("3. PLSDA_composite_scores.csv - 样本综合得分\n")
cat("4. PLSDA_analysis_results.RData - R数据文件\n")
cat("\n>>> 用于RWR网络的综合得分已计算完成。\n")
cat("您可以使用 results_df$Composite_Score 作为节点属性进行网络传播分析。\n")