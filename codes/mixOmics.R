# 1. 加载必要的包
if (!require("mixOmics")) install.packages("mixOmics")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("pheatmap")) install.packages("pheatmap")
library(devtools)
library(mixOmics)
library(ggplot2)
library(dplyr)
library(pheatmap)

# 2. 数据准备与检查
cat("数据维度:", dim(all_mut_robust), "\n")
cat("疾病类型:", unique(all_mut_robust$Disease), "\n")
cat("样本数量:", nrow(all_mut_robust), "\n")

# 3. 将15个指标按四个层面分组
# 根据您的描述，我将指标分组为：
# 序列层面: Co.evolution, Entropy, Consurf_Score
# 结构层面: RASA, ddG, Hydrophobicity
# 网络层面: Betweenness, Closeness, Eigenvector, CC
# 动力学层面: Effectiveness, Sensitivity, MSF, DFI, Stiffness

# 提取各层面的数据矩阵
X_seq <- as.matrix(all_mut_robust[, c("Co.evolution", "Entropy", "Consurf_Score")])
X_struct <- as.matrix(all_mut_robust[, c("RASA", "ddG", "Hydrophobicity")])
X_net <- as.matrix(all_mut_robust[, c("Betweenness", "Closeness", "Eigenvector", "CC")])
X_dyn <- as.matrix(all_mut_robust[, c("Effectiveness", "Sensitivity", "MSF", "DFI", "Stiffness")])

# 为每个矩阵设置行名
rownames(X_seq) <- rownames(X_struct) <- rownames(X_net) <- rownames(X_dyn) <- 
  paste(all_mut_robust$Gene, all_mut_robust$Site, all_mut_robust$Mutation, sep="_")

# 创建数据块列表
X_blocks <- list(
  Sequence = X_seq,
  Structure = X_struct,
  Network = X_net,
  Dynamics = X_dyn
)

# 响应变量
Y <- factor(all_mut_robust$Disease)
names(Y) <- rownames(X_seq)

# 4. 检查各数据块的基本信息
cat("\n>>> 各数据块维度:\n")
for (i in 1:length(X_blocks)) {
  cat(names(X_blocks)[i], ":", dim(X_blocks[[i]]), "\n")
}

# 5. DIABLO模型构建
cat("\n>>> 开始构建DIABLO模型...\n")

# 5.1 设置设计矩阵（定义数据块之间的连接强度）
# 对角线为0，非对角线表示连接强度（0-1之间，值越大表示连接越强）
design <- matrix(0.1, ncol = length(X_blocks), nrow = length(X_blocks),
                 dimnames = list(names(X_blocks), names(X_blocks)))
diag(design) <- 0  # 对角线设为0

cat("\n设计矩阵:\n")
print(design)

# 5.2 设置每个数据块保留的变量数（通过交叉验证选择）
# 这里我们先尝试一个简单的设置，实际分析中应通过交叉验证优化
list.keepX <- list(
  Sequence = rep(2, 2),      # 序列层面：2个成分，每个成分保留2个变量
  Structure = rep(2, 2),     # 结构层面：2个成分，每个成分保留2个变量
  Network = rep(2, 2),       # 网络层面：2个成分，每个成分保留2个变量
  Dynamics = rep(3, 2)       # 动力学层面：2个成分，每个成分保留3个变量
)

# 5.3 运行DIABLO
set.seed(123)  # 设置随机种子以保证可重复性
diablo_model <- block.splsda(
  X = X_blocks,
  Y = Y,
  ncomp = 2,          # 潜在成分数
  keepX = list.keepX,
  design = design,
  scale = FALSE       # 数据已经标准化
)

# 6. 模型评估
cat("\n>>> DIABLO模型结果概览:\n")
print(diablo_model)

# 6.1 评估模型性能（使用5折交叉验证）
cat("\n>>> 进行交叉验证评估模型性能...\n")
set.seed(123)
perf_diablo <- perf(diablo_model, validation = "Mfold",folds = 5, nrepeat = 10)

cat("\n交叉验证错误率:\n")
print(perf_diablo$error.rate)

# 绘制错误率图
plot(perf_diablo)

# 6.2 选择最佳成分数
# 基于交叉验证错误率选择成分数
if (!is.null(perf_diablo$choice.ncomp$ncomp)) {
  best_ncomp <- perf_diablo$choice.ncomp$ncomp
  cat("\n>>> 基于交叉验证的最佳成分数:", best_ncomp, "\n")
} else {
  best_ncomp <- 2
  cat("\n>>> 使用默认成分数:", best_ncomp, "\n")
}

# 7. 模型可视化
cat("\n>>> 生成模型诊断图...\n")

# 7.1 样本图（前两个成分）
plotIndiv(diablo_model,
          ind.names = FALSE,
          legend = TRUE,
          title = "DIABLO样本图",
          subtitle = "不同疾病类型的样本分布")

# 7.2 各层面样本图
plotIndiv(diablo_model,
          ind.names = FALSE,
          legend = TRUE,
          title = "DIABLO各层面样本图",
          subtitle = "按数据块展示",
          blocks = c("Sequence", "Structure", "Network", "Dynamics"),
          group = Y)

# 7.3 变量图
plotVar(diablo_model,
        var.names = TRUE,
        style = "ggplot2",
        title = "DIABLO变量相关图")

# 7.4 各层面变量图
plotVar(diablo_model,
        var.names = TRUE,
        style = "ggplot2",
        title = "DIABLO各层面变量图",
        blocks = c("Sequence", "Structure", "Network", "Dynamics"))

# 8. 特征重要性与权重提取
cat("\n>>> 提取特征权重...\n")

# 8.1 提取各数据块在第一个成分的载荷
weights_list <- list()

for (block_name in names(X_blocks)) {
  # 获取载荷矩阵
  loadings <- diablo_model$loadings[[block_name]]
  
  if (!is.null(loadings)) {
    # 提取第一个成分的载荷
    comp1_loadings <- loadings[, 1]
    
    # 创建权重数据框
    weights_df <- data.frame(
      Feature = rownames(loadings),
      Block = block_name,
      Loading_comp1 = comp1_loadings,
      Abs_Loading = abs(comp1_loadings)
    )
    
    # 按载荷绝对值排序
    weights_df <- weights_df[order(-weights_df$Abs_Loading), ]
    
    weights_list[[block_name]] <- weights_df
    
    cat(paste0("\n", block_name, "层面特征权重(按绝对值排序):\n"))
    print(weights_df)
  }
}

# 8.2 提取变量重要性投影（VIP）
if (!is.null(diablo_model$loadings)) {
  # 计算各变量的VIP值（基于所有成分）
  vip_values <- vip(diablo_model)
  
  # 创建VIP数据框
  vip_df <- data.frame(
    Feature = rownames(vip_values),
    VIP = vip_values[, 1]  # 第一个成分的VIP
  )
  
  # 按VIP值排序
  vip_df <- vip_df[order(-vip_df$VIP), ]
  
  cat("\n>>> 变量重要性投影(VIP)排序:\n")
  print(head(vip_df, 15))
}

# 8.3 提取各数据块的贡献权重
cat("\n>>> 各数据块对模型的贡献:\n")
block_contributions <- diablo_model$prop_expl_var
print(block_contributions)

# 可视化各数据块贡献
if (!is.null(block_contributions)) {
  contrib_df <- data.frame(
    Block = colnames(block_contributions),
    Contribution = block_contributions[1, ]  # 第一个成分的贡献
  )
  
  contrib_plot <- ggplot(contrib_df, aes(x = reorder(Block, -Contribution), y = Contribution)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    labs(x = "数据块", y = "解释方差比例",
         title = "各数据块对DIABLO模型的贡献",
         subtitle = "第一个潜在成分") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(contrib_plot)
}

# 9. 计算综合得分（DIABLO方式）
cat("\n>>> 计算每个样本的综合得分...\n")

# 9.1 获取DIABLO的潜在变量得分
diablo_scores <- predict(diablo_model, newdata = X_blocks)$variates

# 使用第一个成分的得分作为综合指标（或者可以加权组合多个成分）
if (!is.null(diablo_scores$Sequence)) {
  # 提取每个数据块第一个成分的得分
  seq_score <- diablo_scores$Sequence[, 1]
  struct_score <- diablo_scores$Structure[, 1]
  net_score <- diablo_scores$Network[, 1]
  dyn_score <- diablo_scores$Dynamics[, 1]
  
  # 使用各数据块的贡献权重进行加权
  if (!is.null(block_contributions)) {
    weights <- block_contributions[1, ]
    
    # 标准化权重
    weights_norm <- weights / sum(weights)
    
    # 计算加权综合得分
    composite_score_diablo <- 
      seq_score * weights_norm["Sequence"] +
      struct_score * weights_norm["Structure"] +
      net_score * weights_norm["Network"] +
      dyn_score * weights_norm["Dynamics"]
  } else {
    # 如果无法获取贡献权重，使用简单平均
    composite_score_diablo <- (seq_score + struct_score + net_score + dyn_score) / 4
  }
  
  # 9.2 创建结果数据框
  results_diablo <- data.frame(
    Sample = names(composite_score_diablo),
    Gene = all_mut_robust$Gene,
    Disease = all_mut_robust$Disease,
    Site = all_mut_robust$Site,
    Mutation = all_mut_robust$Mutation,
    Composite_Score = as.numeric(composite_score_diablo),
    Seq_Score = seq_score,
    Struct_Score = struct_score,
    Net_Score = net_score,
    Dyn_Score = dyn_score
  )
  
  # 按疾病类型和综合得分排序
  results_diablo <- results_diablo[order(results_diablo$Disease, -results_diablo$Composite_Score), ]
  
  cat("\n>>> DIABLO综合得分前10个样本:\n")
  print(head(results_diablo, 10))
  
  # 9.3 可视化综合得分
  diablo_boxplot <- ggplot(results_diablo, aes(x = Disease, y = Composite_Score, fill = Disease)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    labs(x = "疾病类型", y = "综合得分", 
         title = "不同疾病类型的综合得分分布(DIABLO)",
         subtitle = "基于四个层面整合的得分") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  print(diablo_boxplot)
  
  # 9.4 各层面得分热图
  score_matrix <- as.matrix(results_diablo[, c("Seq_Score", "Struct_Score", "Net_Score", "Dyn_Score")])
  rownames(score_matrix) <- results_diablo$Sample
  
  # 选择得分差异最大的前20个样本
  top_samples_diablo <- head(results_diablo[order(-abs(results_diablo$Composite_Score)), ], 20)
  top_score_matrix <- score_matrix[rownames(score_matrix) %in% top_samples_diablo$Sample, ]
  
  # 添加注释
  annotation_row_diablo <- data.frame(
    Disease = top_samples_diablo$Disease[match(rownames(top_score_matrix), top_samples_diablo$Sample)],
    Composite_Score = top_samples_diablo$Composite_Score[match(rownames(top_score_matrix), top_samples_diablo$Sample)]
  )
  rownames(annotation_row_diablo) <- rownames(top_score_matrix)
  
  pheatmap(top_score_matrix,
           scale = "column",  # 按列标准化以比较不同层面
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           annotation_row = annotation_row_diablo,
           show_rownames = TRUE,
           show_colnames = TRUE,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           main = "高综合得分样本的各层面得分热图(DIABLO)",
           fontsize_row = 8,
           fontsize_col = 10)
}

# 10. 网络图展示DIABLO结果
cat("\n>>> 生成DIABLO网络图...\n")

# 10.1 变量相关网络
network_result <- network(diablo_model,
                          blocks = c(1, 2, 3, 4),  # 所有数据块
                          color.node = c("Sequence" = "lightgreen", 
                                         "Structure" = "orange", 
                                         "Network" = "lightblue", 
                                         "Dynamics" = "pink"),
                          cutoff = 0.7,  # 相关性阈值
                          save = "png",
                          name.save = "DIABLO_network")

# 10.2 圆形布局网络图
plotArrow(diablo_model,
          ind.names = FALSE,
          legend = TRUE,
          title = "DIABLO箭头图")

# 11. 保存结果
cat("\n>>> 保存DIABLO分析结果...\n")

# 保存权重
for (block_name in names(weights_list)) {
  filename <- paste0("DIABLO_weights_", block_name, ".csv")
  write.csv(weights_list[[block_name]], filename, row.names = FALSE)
  cat("保存文件:", filename, "\n")
}

# 保存综合得分
if (exists("results_diablo")) {
  write.csv(results_diablo, "DIABLO_composite_scores.csv", row.names = FALSE)
  cat("保存文件: DIABLO_composite_scores.csv\n")
}

# 保存VIP值
if (exists("vip_df")) {
  write.csv(vip_df, "DIABLO_VIP_values.csv", row.names = FALSE)
  cat("保存文件: DIABLO_VIP_values.csv\n")
}

# 保存模型对象
save(diablo_model, perf_diablo, weights_list, results_diablo,
     file = "DIABLO_analysis_results.RData")
cat("保存文件: DIABLO_analysis_results.RData\n")

cat("\n>>> DIABLO分析完成!\n")
cat(">>> 综合得分已计算完成，可用于RWR网络传播分析。\n")
cat(">>> 建议使用 results_diablo$Composite_Score 作为节点属性。\n")