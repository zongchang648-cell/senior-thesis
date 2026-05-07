#给PRS结果附上列名，以便区分是什么疾病什么蛋白的什么突变产生的扰动
# 需要的包（如果未安装请先安装：install.packages(c("readxl","openxlsx"))）
library(readxl)
library(openxlsx)  # 仅用于可选保存

# 文件路径（按你的示例）
responses_path <- "C:\\\\毕业论文\\\\数据处理及绘图\\\\ppi_build\\\\python_prs_output\\\\responses.xlsx"
pca_path       <- "C:\\\\毕业论文\\\\数据\\\\PCA_compressed_data.xlsx"

# 1) 读入 responses.xlsx，第一列为行名（基因名）
prs_raw <- readxl::read_excel(responses_path, col_names = TRUE)
# 把 tibble 转为 data.frame
prs_raw <- as.data.frame(prs_raw, stringsAsFactors = FALSE)
rowname<-prs_raw[,1]
rownames(prs_raw)<-rowname
prs_raw<-prs_raw[,-1]


# 1) 读入 PCA 文件（全部读为字符，避免类型问题）
pca_path <- "C:/毕业论文/数据/PCA_compressed_data.xlsx"
pca_df <- read_excel(pca_path, col_types = "text")

# 2) 检查必要列是否存在
required_cols <- c("Disease", "Gene", "Site", "Mutation")
missing_cols <- setdiff(required_cols, names(pca_df))
if(length(missing_cols) > 0){
  stop("PCA 文件缺少必要列: ", paste(missing_cols, collapse = ", "))
}

# 3) 把每列强制为字符并处理 NA（用 "NA" 占位）
disease <- ifelse(is.na(pca_df[["Disease"]]), "NA", as.character(pca_df[["Disease"]]))
gene    <- ifelse(is.na(pca_df[["Gene"]]),    "NA", as.character(pca_df[["Gene"]]))
site    <- ifelse(is.na(pca_df[["Site"]]),    "NA", as.character(pca_df[["Site"]]))
mut     <- ifelse(is.na(pca_df[["Mutation"]]),"NA", as.character(pca_df[["Mutation"]]))

# 4) 生成新列名（顺序：Disease_Gene_Site_Mutation）
new_colnames <- paste(disease, gene, site, mut, sep = "_")

# 5) 可选：清理字符（去空格、把文件名/列名中可能导致问题的字符替换为短横）
new_colnames <- gsub("\\s+", "", new_colnames)                     # 去掉所有空白
new_colnames <- gsub("[/\\\\:*?\"<>|]", "-", new_colnames)         # 把一些特殊符号替为 -

# 6) 检查 prs_raw 是否存在并分配列名
if(!exists("prs_raw")){
  stop("对象 prs_raw 在当前环境中不存在。请先把 prs_raw 读入 R 并确保它是一个 data.frame 或 matrix。")
}

# 确保 prs_raw 可设置列名（matrix/data.frame）
if(!is.data.frame(prs_raw) && !is.matrix(prs_raw)){
  stop("prs_raw 必须是 data.frame 或 matrix。")
}

# 7) 按可对应的最小长度进行赋值（避免行/列数不一致导致错误）
n_pca <- length(new_colnames)
n_prs <- ncol(prs_raw)
n_assign <- min(n_pca, n_prs)

if(n_pca != n_prs){
  warning(sprintf("PCA 行数 (%d) 与 prs_raw 列数 (%d) 不相等。将为前 %d 列赋新列名，剩余列保持不变。", 
                  n_pca, n_prs, n_assign))
}

colnames(prs_raw)[seq_len(n_assign)] <- new_colnames[seq_len(n_assign)]

# 8) 快速查看结果
#message("已为 prs_raw 的前 ", n_assign, " 列赋新列名（格式：Disease_Gene_Site_Mutation）。")
#message("示例新列名（前10个）：")
#print(head(colnames(prs_raw), 10))
prs_raw_origin<-prs_raw
prs_raw<-as.data.frame(prs_raw)
prs_raw[] <- lapply(prs_raw, function(col) {
  # 把字符转为数值（抑制转换时的警告），并检查是否有非数字被转为 NA
  num <- suppressWarnings(as.numeric(col))
  if (any(is.na(num) & !is.na(col))) {
    warning("注意：存在无法转为数值的元素，它们已被设置为 NA（见具体位置）。")
  }
  abs(num)
})
###################################################################################################
#根据扰动得分与疾病对节点进行排序
# 需要包
#if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# ---------- 1. 把 prs_raw 转为数值矩阵并取绝对值 ----------
# 假定 prs_raw 已在工作区，且有 rownames/prs_raw 和 colnames(prs_raw)
# 先备份（以防万一）
prs_backup <- prs_raw

# 转换为矩阵并把每个元素转为 numeric，再取绝对值
mat_char <- as.matrix(prs_raw)
mat_num <- apply(mat_char, 2, function(col) suppressWarnings(as.numeric(col)))
# apply 可能会在单列上返回向量，确保是矩阵（当只有一列时）
mat_num <- as.matrix(mat_num)
rownames(mat_num) <- rownames(mat_char)
colnames(mat_num) <- colnames(mat_char)

# 检查是否有非数值被转为 NA
na_count <- sum(is.na(mat_num) & !is.na(mat_char))
if (na_count > 0) warning("转换过程中发现 ", na_count, " 个无法转换为数值的元素，被设为 NA。")

# 取绝对值
mat_num <- abs(mat_num)

# ---------- 2. 从列名解析 Disease（列名格式假设为 Disease_...） ----------
if (is.null(colnames(mat_num))) stop("prs_raw (或 mat_num) 没有列名，无法从列名解析 Disease。")
disease_vec <- sapply(strsplit(colnames(mat_num), "_"), `[`, 1)
unique_diseases <- unique(disease_vec)
message("解析到的 Disease 类型：", paste(unique_diseases, collapse = ", "))

# 期望三类
expected <- c("Cancer", "NDD", "OR")
missing_expected <- setdiff(expected, unique_diseases)
if (length(missing_expected) > 0) {
  warning("未找到以下期望的 Disease 分组: ", paste(missing_expected, collapse = ", "), 
          "。将对存在的分组进行处理。")
}

# ---------- 3. 按 Disease 划分列，计算每行均值并排序 ----------
# 辅助函数：给定列索引，返回排序后的子矩阵与对应均值
make_group_sorted <- function(mat, col_idx) {
  subm <- mat[, col_idx, drop = FALSE]
  # 计算每行均值（忽略 NA）
  row_means <- rowMeans(subm, na.rm = TRUE)
  # 按均值降序排列（从高到低）
  ord <- order(row_means, decreasing = TRUE)
  list(
    mat_sorted = subm[ord, , drop = FALSE],
    means_sorted = row_means[ord],
    order_idx = ord
  )
}

# 取得每组的列索引（若某组不存在，返回 integer(0)）
idx_cancer <- which(disease_vec == "Cancer")
idx_ndd    <- which(disease_vec == "NDD")
idx_or     <- which(disease_vec == "OR")

res_cancer <- if (length(idx_cancer) > 0) make_group_sorted(mat_num, idx_cancer) else NULL
res_ndd    <- if (length(idx_ndd) > 0)    make_group_sorted(mat_num, idx_ndd)    else NULL
res_or     <- if (length(idx_or) > 0)     make_group_sorted(mat_num, idx_or)     else NULL

# 确认行数（应该一致）
n_rows <- nrow(mat_num)
message("原始行数（基因数）: ", n_rows, ". 你期望每个子数据框为 1266 行。")
if (n_rows != 1266) warning("当前行数不是 1266（实际为 ", n_rows, "）。将按实际行数绘图/排序。")

# 若某组不存在，创建长度为 n_rows 的 NA 向量以便绘图（保持 x 轴对齐）
means_cancer <- if (!is.null(res_cancer)) res_cancer$means_sorted else rep(NA_real_, n_rows)
means_ndd    <- if (!is.null(res_ndd))    res_ndd$means_sorted    else rep(NA_real_, n_rows)
means_or     <- if (!is.null(res_or))     res_or$means_sorted     else rep(NA_real_, n_rows)

# 若某组的排序结果行数与 n_rows 不一致（理论上应该一致，因为我们对所有行都计算均值），提醒
if (length(means_cancer) != n_rows) warning("Cancer 组排序后的行数与原始行数不一致。")
if (length(means_ndd)    != n_rows) warning("NDD 组排序后的行数与原始行数不一致。")
if (length(means_or)     != n_rows) warning("OR 组排序后的行数不一致。")

# ---------- 4. 准备用于绘图的数据框 ----------
# x 轴为 1:n_rows（代表名次 1..n）
df_plot <- data.frame(
  x = rep(seq_len(n_rows), times = 3),
  y = c(means_cancer, means_ndd, means_or),
  group = factor(rep(c("Cancer", "NDD", "OR"), each = n_rows), levels = c("Cancer", "NDD", "OR"))
)

# 去掉所有 y 为 NA 的点（如果某组不存在，会产生 NA）
df_plot <- df_plot[!is.na(df_plot$y), ]

# ---------- 5. 绘图：散点图，去除网格线 ----------
p <- ggplot(df_plot, aes(x = x, y = y, color = group)) +
  geom_point(alpha = 0.8, size = 1) +
  labs(x = "Rank ", y = "Row mean (sorted descending)", color = "Disease") +
  theme_minimal() +
  theme(
    #panel.grid.major = element_blank(),   # 去除主网格线
    #panel.grid.minor = element_blank(),   # 去除次网格线
    panel.border = element_blank()
  )

print(p)

# （可选）保存图片
# ggsave("prs_group_means_scatter.png", p, width = 8, height = 5, dpi = 300)


###############################################################################################
#计算出按疾病和得分排序的蛋白质基因，并用于后续的药物重定位研究
the_rownames<-rownames(prs_raw)

cancer_order<-res_cancer$order_idx
cancer_score<-res_cancer$means_sorted
cancer_store<-c()
for (i in cancer_order){
  cancer_store<-c(cancer_store,the_rownames[i])
}
cancer<-cbind(cancer_store,cancer_score)
cancer<-as.data.frame(cancer)


ndd_order<-res_ndd$order_idx
ndd_score<-res_ndd$means_sorted
ndd_store<-c()
for (i in ndd_order){
  ndd_store<-c(ndd_store,the_rownames[i])
}
ndd<-cbind(ndd_store,ndd_score)
ndd<-as.data.frame(ndd)


or_order<-res_or$order_idx
or_score<-res_or$means_sorted
or_store<-c()
for (i in or_order){
  or_store<-c(or_store,the_rownames[i])
}
or<-cbind(or_store,or_score)
or<-as.data.frame(or)

write_xlsx(cancer,'C:\\毕业论文\\数据处理及绘图\\ppi_build\\cancer_prs.xlsx')
write_xlsx(ndd,'C:\\毕业论文\\数据处理及绘图\\ppi_build\\ndd_prs.xlsx')
write_xlsx(or,'C:\\毕业论文\\数据处理及绘图\\ppi_build\\or_prs.xlsx')


##################################################################################
#画出prs_raw对应的热图

# 安装/加载包
# install.packages("pheatmap")
library(pheatmap)

# 1) 将数据框转成矩阵
prs_mat <- as.matrix(prs_raw)

# 2) 转成数值型（防止读入时是字符）
mode(prs_mat) <- "numeric"

# 3) 处理 0 值，避免 -log10(0) 变成 Inf
pseudo <- 1e-10
prs_log <- -log10(prs_mat + pseudo)

# 4) 画热图
pheatmap(
  prs_log,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  fontsize = 8,
  main = "-log10 transformed heatmap"
)




