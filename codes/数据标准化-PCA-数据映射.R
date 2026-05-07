library(readxl)
library(writexl)
library(ggplot2)
library(tidyr)
library(dplyr)
all_mut<-read_xlsx('C:\\毕业论文\\数据\\ALL_DATA.xlsx',sheet = 1)
# 先计算需要提取的列索引
start_col <- ncol(all_mut) - 14  # 后15列的开始位置
end_col <- ncol(all_mut)         # 最后一列

# 1. 提取后15列数据
data_all <- all_mut[, start_col:end_col]

# 2. 将所有列从字符型转换为数值型，同时处理科学计数法
data_all <- data_all %>%
  mutate(across(everything(), ~ as.numeric(.)))

data_all<-cbind(all_mut[1:4],data_all)
data_all_copy<-data_all

# 对最后15列进行z-score标准化
for (i in (ncol(data_all) - 14):ncol(data_all)) {
  # 计算每列的均值和标准差
  col_mean <- mean(data_all[[i]], na.rm = TRUE)
  col_sd <- sd(data_all[[i]], na.rm = TRUE)
  
  # 计算z-score
  data_all[[i]] <- (data_all[[i]] - col_mean) / col_sd
}

######################PCA########################
# 创建分组并分别进行PCA分析
# 提取各组数据
sequence_data <- data_all[, c(5, 6, 7)]
structure_data <- data_all[, c(8, 9, 14)]
network_data <- data_all[, c(10, 11, 12, 13)]  # 仅提取，不进行PCA
dynamic_data <- data_all[, c(15, 16, 17, 18, 19)]

# 对sequence组进行PCA
if (ncol(sequence_data) >= 2) {
  pca_sequence <- prcomp(sequence_data, scale. = TRUE, center = TRUE)
  cat("=== Sequence组PCA结果 ===\n")
  print(summary(pca_sequence))
  
  # 输出方差解释
  cat("\nSequence组主成分方差解释：\n")
  print(pca_sequence$sdev^2)
} else {
  cat("Sequence组数据不足，无法进行PCA\n")
}

# 对structure组进行PCA
if (ncol(structure_data) >= 2) {
  pca_structure <- prcomp(structure_data, scale. = TRUE, center = TRUE)
  cat("\n=== Structure组PCA结果 ===\n")
  print(summary(pca_structure))
  
  # 输出方差解释
  cat("\nStructure组主成分方差解释：\n")
  print(pca_structure$sdev^2)
} else {
  cat("\nStructure组数据不足，无法进行PCA\n")
}

# 对network组进行PCA
if (ncol(network_data) >= 2) {
  pca_network <- prcomp(network_data, scale. = TRUE, center = TRUE)
  cat("\n===Network组PCA结果 ===\n")
  print(summary(pca_network))
  
  # 输出方差解释
  cat("\nNetwork组主成分方差解释：\n")
  print(pca_network$sdev^2)
} else {
  cat("\nNetwork组数据不足，无法进行PCA\n")
}

# 对dynamic组进行PCA
if (ncol(dynamic_data) >= 2) {
  pca_dynamic <- prcomp(dynamic_data, scale. = TRUE, center = TRUE)
  cat("\n=== Dynamic组PCA结果 ===\n")
  print(summary(pca_dynamic))
  
  # 输出方差解释
  cat("\nDynamic组主成分方差解释：\n")
  print(pca_dynamic$sdev^2)
} else {
  cat("\nDynamic组数据不足，无法进行PCA\n")
}

# 如果需要将PCA结果保存到变量中，可以创建列表存储
pca_results <- list(
  sequence = if (exists("pca_sequence")) pca_sequence else NULL,
  structure = if (exists("pca_structure")) pca_structure else NULL,
  dynamic = if (exists("pca_dynamic")) pca_dynamic else NULL
)

# 可选：查看各组数据的基本信息
cat("\n=== 各组数据维度 ===\n")
cat("Sequence组:", dim(sequence_data), "\n")
cat("Structure组:", dim(structure_data), "\n")
cat("Network组:", dim(network_data), "\n")
cat("Dynamic组:", dim(dynamic_data), "\n")

#############合并##########################
# 继续你的代码，进行各组数据压缩和合并

# 1. 对sequence_data：提取PC1（1维）
sequence_compressed <- data.frame(seq_PC1 = pca_sequence$x[, 1])

# 2. 对structure_data：提取PC1和PC2（2维）
structure_compressed <- data.frame(
  struct_PC1 = pca_structure$x[, 1],
  struct_PC2 = pca_structure$x[, 2]
)

# 3. 对network_data：保持原状（4维），但需要重命名列以避免混淆
# 获取原始列名，并添加前缀
network_cols <- paste0("network_", colnames(network_data))
network_compressed <- network_data
colnames(network_compressed) <- network_cols

# 4. 对dynamic_data：提取PC1, PC2, PC3（3维）
# 先查看dynamic组的PCA结果，确定PC1, PC2, PC3的累计解释比例
cat("Dynamic组累计解释比例：\n")
print(summary(pca_dynamic)$importance[3, 1:3])

# 提取前三个主成分
dynamic_compressed <- data.frame(
  dyn_PC1 = pca_dynamic$x[, 1],
  dyn_PC2 = pca_dynamic$x[, 2],
  dyn_PC3 = pca_dynamic$x[, 3]
)

# 5. 将所有压缩后的数据合并为10列
compressed_data <- cbind(
  sequence_compressed,    # 1列
  structure_compressed,   # 2列
  network_compressed,     # 4列
  dynamic_compressed      # 3列
)

# 检查总列数
cat("\n压缩后的数据维度：", dim(compressed_data), "\n")
cat("应为：行数 x 10列 (1+2+4+3=10)\n")

# 6. 提取data_all的前4列
first_four_cols <- data_all[, 1:4]

# 7. 合并：前4列在左侧，压缩后的10列在右侧
final_data <- cbind(first_four_cols, compressed_data)

# 8. 查看最终数据结构
cat("\n最终数据框维度：", dim(final_data), "\n")
cat("列名：\n")
print(colnames(final_data))

# 9. 验证合并结果
cat("\n=== 合并结果验证 ===\n")
cat("前4列列名：", colnames(final_data)[1:4], "\n")
cat("第5-14列列名（压缩后的10列）：", colnames(final_data)[5:14], "\n")

# 10. 保存最终数据（可选）
# 保存为Excel文件
write_xlsx(final_data, "C:\\毕业论文\\数据\\PCA_compressed_data.xlsx")

# 保存为CSV文件（可选）
#write.csv(final_data, "C:\\毕业论文\\数据\\PCA_compressed_data.csv", row.names = FALSE)

# 11. 创建压缩信息摘要
compression_summary <- data.frame(
  组别 = c("Sequence", "Structure", "Network", "Dynamic", "总计"),
  原始维度 = c(3, 3, 4, 5, 15),
  压缩后维度 = c(1, 2, 4, 3, 10),
  信息保留比例 = c(
    paste0(round(summary(pca_sequence)$importance[2, 1]*100, 1), "%"),
    paste0(round(summary(pca_structure)$importance[3, 2]*100, 1), "%"),
    "100%",
    paste0(round(summary(pca_dynamic)$importance[3, 3]*100, 1), "%"),
    "-"
  ),
  使用的PCs = c("PC1", "PC1+PC2", "原始4列", "PC1+PC2+PC3", "-")
)

print(compression_summary)

# 13. 检查压缩后数据的基本统计信息
cat("\n=== 压缩后数据统计摘要 ===\n")
cat("数据框维度：", dim(final_data), "\n")
cat("\n前6行数据：\n")
print(head(final_data))

cat("\n各列数据类型：\n")
print(sapply(final_data, class))

# 14. 检查是否有缺失值
cat("\n缺失值统计：\n")
missing_counts <- sapply(final_data, function(x) sum(is.na(x)))
print(missing_counts[missing_counts > 0])

if (all(missing_counts == 0)) {
  cat("没有缺失值。\n")
}
# 18. 最终结果输出
cat("\n=== 处理完成 ===\n")
cat("1. 原始数据：", ncol(data_all), "列\n")
cat("2. 压缩后数据：", ncol(final_data), "列\n")
cat("3. 维度减少：", ncol(data_all) - ncol(final_data), "列\n")
cat("4. 最终数据已保存为：PCA_compressed_data.xlsx\n")
cat("5. PCA模型已保存为：pca_models.rds\n")
cat("6. 压缩函数已创建：apply_pca_compression()\n")

##################################PPIN加权#############################################
all_ppi<-read_xlsx("C:\\毕业论文\\数据处理及绘图\\ppi_build\\PPIN_score700_SYMBOL_success.xlsx")
all_ppi<-all_ppi[1:3]
PPIN<-read_xlsx("C:\\毕业论文\\数据处理及绘图\\ppi_build\\top10%_CoreNet.xlsx")

start_col_idx <- 5
end_col_idx <- 14

# 计算每个样本的L2范数

total_score <- sqrt(rowSums(as.matrix(final_data[, start_col_idx:end_col_idx])^2, na.rm = TRUE))


# 将total_score列添加到final_data的第一列后面
final_data <- cbind(final_data[, 1:4], 
                    total_score = total_score,
                    final_data[, start_col_idx:end_col_idx])

#####为PPI边加权
add_ppi_edge_weights <- function(PPIN, all_ppi) {
  # 高效方法：使用规范化基因对
  cat("开始为PPIN添加边权重...\n")
  
  # 确保基因列为字符型
  PPIN$gene1 <- as.character(PPIN$gene1)
  PPIN$gene2 <- as.character(PPIN$gene2)
  all_ppi$gene1 <- as.character(all_ppi$gene1)
  all_ppi$gene2 <- as.character(all_ppi$gene2)
  
  # 创建规范化基因对
  create_norm_pair <- function(g1, g2) {
    # 按字母顺序排序，确保(g1,g2)和(g2,g1)相同
    pair <- ifelse(g1 < g2, paste(g1, g2, sep = "||"), paste(g2, g1, sep = "||"))
    return(pair)
  }
  
  # 为PPIN和all_ppi创建规范化基因对
  PPIN$norm_pair <- create_norm_pair(PPIN$gene1, PPIN$gene2)
  all_ppi$norm_pair <- create_norm_pair(all_ppi$gene1, all_ppi$gene2)
  
  # 创建all_ppi的查找表
  # 如果有重复的规范化基因对，取combined_score的平均值
  library(dplyr)
  lookup_table <- all_ppi %>%
    group_by(norm_pair) %>%
    summarise(combined_score = mean(combined_score, na.rm = TRUE)) %>%
    ungroup()
  
  # 合并PPIN和查找表
  PPIN_complete <- PPIN %>%
    left_join(lookup_table, by = "norm_pair") %>%
    select(gene1, gene2, combined_score)  # 移除norm_pair列
  
  # 统计
  total_edges <- nrow(PPIN_complete)
  matched_edges <- sum(!is.na(PPIN_complete$combined_score))
  match_rate <- matched_edges / total_edges * 100
  
  cat("匹配完成：\n")
  cat(sprintf("  总边数: %d\n", total_edges))
  cat(sprintf("  匹配成功: %d\n", matched_edges))
  cat(sprintf("  匹配率: %.2f%%\n", match_rate))
  
  # 检查是否有完全未匹配的情况
  if (matched_edges == 0) {
    cat("警告：没有成功匹配任何边！\n")
    cat("可能原因：\n")
    cat("  1. PPIN和all_ppi使用不同的基因命名方式\n")
    cat("  2. 数据格式有问题\n")
    cat("  3. 样本没有重叠\n")
  }
  
  return(PPIN_complete)
}

# 使用函数
PPIN_complete <- add_ppi_edge_weights(PPIN, all_ppi)

# 查看结果
cat("\nPPIN_complete的前10行：\n")
print(head(PPIN_complete, 10))






