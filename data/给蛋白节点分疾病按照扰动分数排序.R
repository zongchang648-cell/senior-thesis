# 加载必要的包
library(dplyr)
library(tidyr)
library(VennDiagram)

# 假设 result_matrix 已存在，行名是基因名，列名是字符串
# 例如：result_matrix <- read.csv("your_data.csv", row.names = 1)

# 解析列名，提取疾病和基因信息
col_info <- data.frame(
  full_name = colnames(result_matrix),
  stringsAsFactors = FALSE
) %>%
  separate(full_name, into = c("Disease", "Gene", "Site", "Mutation"), sep = "_", remove = FALSE)

# 将疾病和基因信息添加到列名中便于后续筛选
# 这里我们直接用 col_info 来分组

# 任务1：按疾病分组，计算每行均值，取前50节点
diseases <- unique(col_info$Disease)
top_50_list <- list()

for (d in diseases) {
  # 选择当前疾病的所有列
  cols <- col_info$full_name[col_info$Disease == d]
  # 计算行均值
  mean_vec <- rowMeans(result_matrix[, cols, drop = FALSE], na.rm = TRUE)
  # 降序排序并取前50个节点名
  sorted_names <- names(sort(mean_vec, decreasing = TRUE))
  top_50 <- sorted_names[1:50]  # 若有并列，取前50个
  top_50_list[[d]] <- top_50
}

# 构建数据框 top_50_ns，50行3列
top_50_ns <- data.frame(
  Cancer = top_50_list[["Cancer"]],
  NDD = top_50_list[["NDD"]],
  OR = top_50_list[["OR"]],
  stringsAsFactors = FALSE
)

# 查看结果
#head(top_50_ns)

# 画韦恩图
venn.plot <- venn.diagram(
  x = list(Cancer = top_50_list[["Cancer"]],
           NDD = top_50_list[["NDD"]],
           OR = top_50_list[["OR"]]),
  category.names = c("Cancer", "NDD", "OR"),
  filename = NULL,  # 直接显示在R中，可保存为文件
  output = TRUE,
  imagetype = "png",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = c("red", "green", "blue"),
  fill = c(alpha("red", 0.3), alpha("green", 0.3), alpha("blue", 0.3)),
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("red", "green", "blue")
)
grid.draw(venn.plot)

# 保存韦恩图到文件
# tiff("venn_diseases.tiff", width = 480, height = 480)
# grid.draw(venn.plot)
# dev.off()

# 计算交集并保存
intersect_all <- Reduce(intersect, top_50_list)
intersect_cancer_ndd <- intersect(top_50_list[["Cancer"]], top_50_list[["NDD"]])
intersect_cancer_or <- intersect(top_50_list[["Cancer"]], top_50_list[["OR"]])
intersect_ndd_or <- intersect(top_50_list[["NDD"]], top_50_list[["OR"]])

# 将交集保存为变量
intersections <- list(
  all = intersect_all,
  cancer_ndd = intersect_cancer_ndd,
  cancer_or = intersect_cancer_or,
  ndd_or = intersect_ndd_or
)

write_xlsx(top_50_ns,"C:\\毕业论文\\数据\\top_50_ns.xlsx")

# 任务2：按疾病和基因分组，共27组
# 获取所有基因名
genes <- unique(col_info$Gene)

# 创建组合标签
combinations <- expand.grid(Disease = diseases, Gene = genes, stringsAsFactors = FALSE)
combinations$label <- paste(combinations$Disease, combinations$Gene, sep = "_")

# 初始化列表存储每个组合的前50节点
top_50_gene_disease <- list()

for (i in 1:nrow(combinations)) {
  d <- combinations$Disease[i]
  g <- combinations$Gene[i]
  # 选择同时满足疾病和基因的列
  cols <- col_info$full_name[col_info$Disease == d & col_info$Gene == g]
  if (length(cols) == 0) next  # 若没有突变，跳过
  mean_vec <- rowMeans(result_matrix[, cols, drop = FALSE], na.rm = TRUE)
  sorted_names <- names(sort(mean_vec, decreasing = TRUE))
  top_50 <- sorted_names[1:50]
  top_50_gene_disease[[combinations$label[i]]] <- top_50
}

# 构建数据框，50行27列
# 由于各列可能长度不一致（如果某组突变数少导致排名不足50？但假设每组都有足够突变），这里用数据框合并
top_50_by_gene <- as.data.frame(lapply(top_50_gene_disease, function(x) {
  # 确保每个向量长度为50，不足时补NA
  if(length(x) < 50) {
    c(x, rep(NA, 50 - length(x)))
  } else {
    x[1:50]
  }
}), stringsAsFactors = FALSE)

# 确保列名顺序与combinations$label一致
top_50_by_gene <- top_50_by_gene[, combinations$label, drop = FALSE]

# 查看结果
dim(top_50_by_gene)  # 应为 50 27
head(top_50_by_gene)

# 可选：将结果保存为CSV文件
# write.csv(top_50_ns, "top_50_ns.csv", row.names = FALSE)
# write.csv(top_50_by_gene, "top_50_by_gene_disease.csv", row.names = FALSE)