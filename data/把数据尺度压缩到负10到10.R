

library(readxl)
library(dplyr)

# 1. 导入Excel数据
file_path <- "C:/毕业论文/数据/all_data_numeric.xlsx"
data <- read_xlsx(file_path)

# 2. 确定后15列
n_cols <- ncol(data)
last_15_cols <- (n_cols - 14):n_cols
col_names <- names(data)[last_15_cols]

# 3. 自定义标准化函数：缩放到[-10, 10]，以0为中心
scale_to_range <- function(x) {
  if (sd(x, na.rm = TRUE) == 0) {
    # 如果标准差为0（常数列），直接返回0
    return(rep(0, length(x)))
  }
  
  # 计算当前列的最大绝对值
  max_abs <- max(abs(x), na.rm = TRUE)
  
  # 如果最大绝对值为0，返回0向量
  if (max_abs == 0) {
    return(rep(0, length(x)))
  }
  
  # 标准化到[-1, 1]范围，保留正负号
  scaled <- x / max_abs
  
  # 缩放到[-10, 10]范围
  scaled <- scaled * 10
  
  return(scaled)
}

# 4. 对后15列应用标准化
data_scaled <- data
data_scaled[col_names] <- lapply(data[col_names], scale_to_range)

# 5. 验证标准化结果
# 查看后15列的范围
range_summary <- sapply(data_scaled[col_names], function(x) {
  c(min = min(x, na.rm = TRUE), 
    max = max(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE))
})

print("标准化后各列的范围：")
print(t(range_summary))

# 6. 可选：保存结果到新文件
write.csv(data_scaled, "C:/毕业论文/数据/all_data_scaled.csv", row.names = FALSE)

# 7. 查看数据基本信息
cat("\n数据基本信息：\n")
cat("原始数据维度：", dim(data), "\n")
cat("标准化后数据维度：", dim(data_scaled), "\n")
cat("后15列列名：", paste(col_names, collapse = ", "), "\n")