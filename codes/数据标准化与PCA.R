library(readxl)
library(writexl)
library(dplyr)
all_mut<-read_xlsx('C:\\毕业论文\\数据\\ALL_DATA.xlsx',sheet = 1)
# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

############################检查数据分布
# 提取后15列数据
# 先计算需要提取的列索引
start_col <- ncol(all_mut) - 14  # 后15列的开始位置
end_col <- ncol(all_mut)         # 最后一列

# 1. 提取后15列数据
data_to_plot <- all_mut[, start_col:end_col]

# 2. 将所有列从字符型转换为数值型，同时处理科学计数法
data_to_plot <- data_to_plot %>%
  mutate(across(everything(), ~ as.numeric(.)))

# 3. 将宽数据转换为长数据（便于使用ggplot2分面绘制）
all_mut_long <- gather(
  data_to_plot,
  key = "variable",
  value = "value"
)

# 4. 移除NA值（如果转换过程中产生了NA）
all_mut_long <- all_mut_long %>%
  filter(!is.na(value))

# 5. 绘制分面直方图
ggplot(all_mut_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title = "后15列数值分布直方图",
    x = "数值",
    y = "频数"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    panel.spacing = unit(1, "lines")
  )
##########################################################z-score

feature_cols <- c("Co.evolution","Entropy","Consurf_Score","RASA","ddG",
                  "Betweenness","Closeness","Eigenvector","CC","Hydrophobicity",
                  "Effectiveness","Sensitivity","MSF","DFI","Stiffness")
pop_sd <- function(x) sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))

z_pop <- as.data.frame(
  lapply(data_to_plot[feature_cols], function(x) {
    sd0 <- pop_sd(x)
    if (is.na(sd0) || sd0 == 0) {
      # 如果方差为0，全部设为0（或按你需要处理）
      return(rep(0, length(x)))
    } else {
      (x - mean(x, na.rm = TRUE)) / sd0
    }
  })
)
colnames(z_pop) <- paste0(feature_cols, "_zpop")
all_mut_zpop <- cbind(data_to_plot, z_pop)

# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 提取后15列数据
# 先计算需要提取的列索引
start_col <- ncol(all_mut_zpop) - 14  # 后15列的开始位置
end_col <- ncol(all_mut_zpop)         # 最后一列

# 1. 提取后15列数据
data_to_plot <- all_mut_zpop[, start_col:end_col]

# 2. 将所有列从字符型转换为数值型，同时处理科学计数法
data_to_plot <- data_to_plot %>%
  mutate(across(everything(), ~ as.numeric(.)))

# 3. 将宽数据转换为长数据（便于使用ggplot2分面绘制）
all_mut_long <- gather(
  data_to_plot,
  key = "variable",
  value = "value"
)

# 4. 移除NA值（如果转换过程中产生了NA）
all_mut_long <- all_mut_long %>%
  filter(!is.na(value))

# 5. 绘制分面直方图
ggplot(all_mut_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  labs(
    title = "后15列数值分布直方图",
    x = "数值",
    y = "频数"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    panel.spacing = unit(1, "lines")
  )


##########################################################Robust Scaling
library(ggplot2)
library(tidyr)
library(dplyr)

# 提取后15列数据
start_col <- ncol(all_mut) - 14
end_col <- ncol(all_mut)

# 1. 提取并转换数据
data_to_plot <- all_mut[, start_col:end_col] %>%
  mutate(across(everything(), ~ as.numeric(.)))

# 2. 转换为长格式并移除NA值
all_mut_long <- data_to_plot %>%
  pivot_longer(cols = everything(), 
               names_to = "variable", 
               values_to = "value") %>%
  filter(!is.na(value))

# 3. 计算每个变量的统计信息
variable_stats <- all_mut_long %>%
  group_by(variable) %>%
  summarise(
    mean_val = mean(value),
    median_val = median(value),
    sd_val = sd(value),
    min_val = min(value),
    max_val = max(value),
    n = n()
  )

# 4. 将统计信息合并到长格式数据中
all_mut_long <- all_mut_long %>%
  left_join(variable_stats, by = "variable")

# 5. 绘制带统计信息的密度曲线图
ggplot(all_mut_long, aes(x = value)) +
  # 密度曲线
  geom_density(fill = "steelblue", alpha = 0.4, color = "steelblue", size = 0.8) +
  
  # 添加均值线
  geom_vline(aes(xintercept = mean_val), 
             color = "red", linetype = "dashed", size = 0.7) +
  
  # 添加中位数线
  geom_vline(aes(xintercept = median_val), 
             color = "darkgreen", linetype = "dashed", size = 0.7) +
  
  # 分面显示
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  
  # 添加统计信息标签
  geom_text(data = distinct(all_mut_long, variable, .keep_all = TRUE),
            aes(x = min_val + 0.1*(max_val - min_val),
                y = Inf,
                label = paste0("n=", n, "\nmean=", round(mean_val, 2), 
                               "\nsd=", round(sd_val, 2))),
            hjust = 0, vjust = 1.2, size = 2.5, color = "darkblue") +
  
  labs(
    title = "特征数值分布密度曲线图",
    subtitle = "红色虚线为均值，绿色虚线为中位数",
    x = "数值",
    y = "密度"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 9, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

##################################################robust scaling
# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

# 提取后15列数据（特征列）
start_col <- ncol(all_mut) - 14  # 后15列的开始位置
end_col <- ncol(all_mut)         # 最后一列

# 获取特征列名
feature_names <- names(all_mut)[start_col:end_col]
cat("要处理的特征列（后15列）：\n")
print(feature_names)

# 1. 提取特征数据并进行字符到数值的转换
# 创建专门的函数处理科学计数法和字符转换
clean_and_convert <- function(x) {
  # 如果已经是数值型，直接返回
  if (is.numeric(x)) return(x)
  
  # 确保是字符型
  x_char <- as.character(x)
  
  # 处理科学计数法：统一将E/e替换为e，R可以识别e为科学计数法
  x_char <- tolower(x_char)  # 将E转换为e
  
  # 移除可能存在的空格、换行符等
  x_char <- trimws(x_char)
  
  # 将空字符串或无法识别的值转换为NA
  x_char[x_char == ""] <- NA
  x_char[x_char == "NA"] <- NA
  x_char[x_char == "NaN"] <- NA
  
  # 转换为数值型
  x_num <- suppressWarnings(as.numeric(x_char))
  
  # 检查转换失败的情况
  failed_indices <- which(is.na(x_num) & !is.na(x_char))
  if (length(failed_indices) > 0) {
    cat(sprintf("警告: %d 个值转换失败，示例如下:\n", length(failed_indices)))
    for (i in head(failed_indices, 3)) {
      cat(sprintf("  索引 %d: '%s'\n", i, x_char[i]))
    }
  }
  
  return(x_num)
}

# 应用清理和转换
feature_data <- all_mut[, feature_names]
feature_data <- feature_data %>%
  mutate(across(everything(), clean_and_convert))

# 检查转换后的数据类型和统计摘要
cat("\n=== 数据类型转换结果 ===\n")
type_summary <- sapply(feature_data, function(x) {
  c(
    "数据类型" = class(x)[1],
    "NA数量" = sum(is.na(x)),
    "NA比例" = round(sum(is.na(x)) / length(x) * 100, 2),
    "唯一值数量" = length(unique(x[!is.na(x)])),
    "范围" = ifelse(sum(!is.na(x)) > 0, 
                  paste(round(min(x, na.rm = TRUE), 4), "到", 
                        round(max(x, na.rm = TRUE), 4)),
                  "无数据")
  )
})
print(t(type_summary))

# 2. 应用Robust Scaling函数（保留16位有效数字）
robust_scale <- function(x) {
  # 确保输入是数值型
  x_numeric <- as.numeric(x)
  
  # 移除NA值进行计算
  x_clean <- x_numeric[!is.na(x_numeric)]
  
  if (length(x_clean) == 0) {
    return(rep(NA, length(x_numeric)))
  }
  
  # 计算中位数和IQR（保留16位有效数字）
  median_val <- signif(median(x_clean, na.rm = TRUE), 16)
  iqr_val <- signif(IQR(x_clean, na.rm = TRUE), 16)
  
  # 处理IQR为0的情况
  if (iqr_val == 0) {
    # 如果IQR为0但数据不是常数，使用标准差
    sd_val <- signif(sd(x_clean, na.rm = TRUE), 16)
    if (sd_val > 0) {
      # 应用标准化，结果保留16位有效数字
      result <- signif((x_numeric - median_val) / sd_val, 16)
      # 对NA值保持为NA
      result[is.na(x_numeric)] <- NA
      return(result)
    } else {
      # 如果所有值都相同，则返回0
      result <- rep(0, length(x_numeric))
      result[is.na(x_numeric)] <- NA
      return(result)
    }
  }
  
  # 应用Robust Scaling公式：(x - median) / IQR
  # 结果保留16位有效数字
  result <- signif((x_numeric - median_val) / iqr_val, 16)
  result[is.na(x_numeric)] <- NA
  return(result)
}

# 对每一列应用Robust Scaling
scaled_data <- feature_data %>%
  mutate(across(everything(), robust_scale))

# 3. 检查标准化后的统计特性
cat("\n=== Robust Scaling后各列统计摘要 ===\n")
# 创建统计摘要数据框
stats_after <- data.frame(
  特征列 = feature_names,
  中位数 = sapply(scaled_data, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) signif(median(x_clean), 6) else NA
  }),
  均值 = sapply(scaled_data, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) signif(mean(x_clean), 6) else NA
  }),
  Q1 = sapply(scaled_data, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) signif(quantile(x_clean, 0.25), 6) else NA
  }),
  Q3 = sapply(scaled_data, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) signif(quantile(x_clean, 0.75), 6) else NA
  }),
  IQR = sapply(scaled_data, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) signif(IQR(x_clean), 6) else NA
  }),
  标准差 = sapply(scaled_data, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) signif(sd(x_clean), 6) else NA
  }),
  最小值 = sapply(scaled_data, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) signif(min(x_clean), 6) else NA
  }),
  最大值 = sapply(scaled_data, function(x) {
    x_clean <- x[!is.na(x)]
    if (length(x_clean) > 0) signif(max(x_clean), 6) else NA
  }),
  有效样本数 = sapply(scaled_data, function(x) sum(!is.na(x))),
  stringsAsFactors = FALSE
)

# 打印统计摘要（只对数值列进行格式化）
numeric_cols <- names(stats_after)[-1]  # 排除"特征列"
for (col in numeric_cols) {
  if (col != "有效样本数") {
    stats_after[[col]] <- sapply(stats_after[[col]], function(x) {
      if (!is.na(x)) sprintf("%.6f", x) else "NA"
    })
  }
}

print(stats_after)

# 4. 计算有多少列的IQR为0（常数特征）
constant_features <- feature_names[stats_after$IQR == "0.000000" | 
                                     stats_after$IQR == "0" |
                                     as.numeric(stats_after$IQR) == 0]
if (length(constant_features) > 0) {
  cat(sprintf("\n有 %d 个特征的IQR为0（常数特征）:\n", length(constant_features)))
  print(constant_features)
}

# 5. 将数据转换为长格式以便绘图
scaled_long <- scaled_data %>%
  pivot_longer(cols = everything(), 
               names_to = "variable", 
               values_to = "value") %>%
  filter(!is.na(value))  # 移除NA值

# 检查绘图数据
cat(sprintf("\n可用于绘图的有效数据点: %d 个\n", nrow(scaled_long)))

# 6. 绘制Robust Scaling后的密度分布图
if (nrow(scaled_long) > 0) {
  p1 <- ggplot(scaled_long, aes(x = value)) +
    geom_density(fill = "steelblue", alpha = 0.6, color = "steelblue", size = 0.8) +
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), color = "darkgreen", linetype = "dashed", size = 0.5, alpha = 0.7) +
    labs(
      title = "Robust Scaling后特征分布密度曲线图",
      subtitle = "红色虚线：中位数=0，绿色虚线：Q1=-0.5，Q3=+0.5",
      x = "Robust标准化值",
      y = "密度"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, margin = margin(b = 10)),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 9, face = "bold"),
      panel.spacing = unit(1, "lines")
    )
  
  print(p1)
} else {
  cat("警告: 没有足够的数据绘制密度分布图\n")
}

# 7. 可选：绘制标准化前后对比图
# 创建标准化前后的长格式数据
original_long <- feature_data %>%
  pivot_longer(cols = everything(), 
               names_to = "variable", 
               values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(type = "原始数据")

scaled_long_with_type <- scaled_long %>%
  mutate(type = "Robust标准化后")

# 合并数据
combined_long <- rbind(original_long, scaled_long_with_type)

# 绘制对比图
if (nrow(combined_long) > 0) {
  p2 <- ggplot(combined_long, aes(x = value, fill = type)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    scale_fill_manual(values = c("原始数据" = "darkorange", "Robust标准化后" = "steelblue")) +
    labs(
      title = "特征标准化前后分布对比",
      x = "数值",
      y = "密度",
      fill = "数据类型"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 9, face = "bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "bottom"
    )
  
  print(p2)
} else {
  cat("警告: 没有足够的数据绘制对比图\n")
}

# 8. 创建统计摘要图
# 计算每个特征的统计量
if (nrow(scaled_long) > 0) {
  summary_stats <- scaled_long %>%
    group_by(variable) %>%
    summarise(
      median_val = median(value, na.rm = TRUE),
      mean_val = mean(value, na.rm = TRUE),
      q1_val = quantile(value, 0.25, na.rm = TRUE),
      q3_val = quantile(value, 0.75, na.rm = TRUE),
      iqr_val = IQR(value, na.rm = TRUE),
      min_val = min(value, na.rm = TRUE),
      max_val = max(value, na.rm = TRUE),
      n = n()
    ) %>%
    ungroup()
  
  # 绘制统计摘要图
  p3 <- ggplot(summary_stats, aes(x = reorder(variable, median_val))) +
    geom_boxplot(aes(ymin = min_val, lower = q1_val, middle = median_val, 
                     upper = q3_val, ymax = max_val),
                 stat = "identity", fill = "lightblue", alpha = 0.7) +
    geom_point(aes(y = mean_val), color = "red", size = 2) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = c(-0.5, 0.5), color = "darkgreen", linetype = "dashed", alpha = 0.5) +
    labs(
      title = "Robust标准化后特征统计摘要",
      subtitle = "箱线图显示分布范围，红点表示均值",
      x = "特征",
      y = "Robust标准化值"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9)
    )
  
  print(p3)
} else {
  cat("警告: 没有足够的数据绘制统计摘要图\n")
}

# 9. 保存标准化后的数据
# 将标准化后的特征数据合并回原始数据框
all_mut_robust <- all_mut
all_mut_robust[, feature_names] <- scaled_data

# 保存为CSV文件（保留16位有效数字）
write.csv(all_mut_robust, "C:\\毕业论文\\数据\\all_mut_robust_scaled.csv", row.names = FALSE)
cat("\nRobust标准化后的数据已保存为 'C:\\毕业论文\\数据\\all_mut_robust_scaled.csv'\n")

# 10. 生成Robust Scaling报告
cat("\n=== Robust Scaling处理报告 ===\n")
cat(sprintf("处理的特征数量: %d\n", length(feature_names)))
cat(sprintf("常数特征数量 (IQR=0): %d\n", length(constant_features)))

# 检查标准化后中位数接近0的特征数量
median_check <- sum(abs(as.numeric(stats_after$中位数)) < 0.01, na.rm = TRUE)
cat(sprintf("标准化后中位数接近0的特征: %d\n", median_check))

# 检查标准化后IQR接近1的特征数量（非常数特征）
iqr_vals <- as.numeric(stats_after$IQR)
iqr_check <- sum(abs(iqr_vals - 1) < 0.01 & iqr_vals > 0, na.rm = TRUE)
cat(sprintf("标准化后IQR接近1的特征: %d\n", iqr_check))

# 检查哪些特征的IQR不为1（非常数特征）
non_constant_features <- feature_names[iqr_vals != 0 & !is.na(iqr_vals)]
if (length(non_constant_features) > 0) {
  cat("\n非常数特征的IQR统计:\n")
  iqr_stats <- stats_after %>%
    filter(特征列 %in% non_constant_features) %>%
    select(特征列, IQR) %>%
    mutate(IQR_num = as.numeric(IQR)) %>%
    arrange(desc(IQR_num))
  
  # 只显示数值部分
  print(iqr_stats[, c("特征列", "IQR")])
  
  # 绘制IQR分布图
  if (nrow(iqr_stats) > 0) {
    p4 <- ggplot(iqr_stats, aes(x = reorder(特征列, IQR_num), y = IQR_num)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
      labs(
        title = "非常数特征的IQR值",
        x = "特征",
        y = "IQR"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(hjust = 0.5, size = 14)
      )
    
    print(p4)
  }
}

# 11. 额外：精度检查
cat("\n=== 精度检查 ===\n")
cat("原始数据类型: 字符型\n")
cat("转换后数据类型: 数值型\n")
cat("有效数字保留: 16位\n")

# 检查科学计数法转换示例
if (nrow(feature_data) > 0) {
  # 查找包含e/E的原始数据
  sci_notation_examples <- character(0)
  for (col in feature_names) {
    # 获取原始数据（字符型）
    orig_col <- all_mut[[col]]
    # 查找包含e或E的值
    sci_idx <- grep("e|E", orig_col, ignore.case = TRUE)
    if (length(sci_idx) > 0) {
      sci_notation_examples <- c(sci_notation_examples, 
                                 paste(col, ":", orig_col[sci_idx[1]]))
    }
  }
  
  if (length(sci_notation_examples) > 0) {
    cat("\n科学计数法转换示例:\n")
    for (example in head(sci_notation_examples, 3)) {
      cat(example, "\n")
    }
  } else {
    cat("\n未发现科学计数法数据\n")
  }
}

# 12. 验证转换精度
cat("\n=== 转换精度验证 ===\n")
cat("随机选择3个特征进行精度验证:\n")
set.seed(123)  # 确保可重复性
sample_features <- sample(feature_names, min(3, length(feature_names)))

for (feat in sample_features) {
  # 获取原始数据（字符型）
  orig_vals <- as.character(all_mut[[feat]])
  # 获取转换后的数值
  conv_vals <- feature_data[[feat]]
  
  # 随机选择3个非NA值进行验证
  non_na_idx <- which(!is.na(conv_vals))
  if (length(non_na_idx) >= 3) {
    test_idx <- sample(non_na_idx, min(3, length(non_na_idx)))
    
    cat(sprintf("\n特征: %s\n", feat))
    for (idx in test_idx) {
      orig <- orig_vals[idx]
      conv <- conv_vals[idx]
      cat(sprintf("  原始值: %s -> 转换后: %.10e\n", 
                  orig, conv))
    }
  }
}


