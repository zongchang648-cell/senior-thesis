
#step_1.构建加权PPI网络并量化初始扰动分数

#这一步骤在前一个脚本中已完成
PPIN_weighed <- PPIN_complete
mutation_score <-final_data
module_type<-read_xlsx("C:\\毕业论文\\modules_fast_greedy.xlsx",sheet = 1)
module_1<-module_type[module_type$module==1,]

#过滤出模块1的PPI加权网络
library(dplyr)

# 简单检查
stopifnot(all(c("gene1","gene2","combined_score") %in% colnames(PPIN_weighed)))
stopifnot(all(c("name","module") %in% colnames(module_1)))

# 去除首尾空白并确保为字符（防止因因子/空格匹配失败）
trim_ws <- function(x) gsub("^\\s+|\\s+$", "", as.character(x))
PPIN_weighed <- PPIN_weighed %>%
  mutate(gene1 = trim_ws(gene1),
         gene2 = trim_ws(gene2))
module_1 <- module_1 %>% mutate(name = trim_ws(name))

# 构建模块基因集合（去重）
genes_in_module <- unique(module_1$name)

# 过滤：只有当两端基因都在 module 列表中才保留
PPIN_module_1 <- PPIN_weighed %>%
  filter(!is.na(gene1) & !is.na(gene2)) %>%
  filter(gene1 %in% genes_in_module & gene2 %in% genes_in_module)

# 保存结果
#saveRDS(PPIN_module_1, file = "PPIN_module_1.rds")
#write.csv(PPIN_module_1, file = "PPIN_module_1.csv", row.names = FALSE)

# 简单信息输出(5794条边)
message("保留的边数：", nrow(PPIN_module_1))

#———————————————————————————————————————————————————————————————————————————————
#step_2. 构建拉普拉斯矩阵
#利用原本用于阐明蛋白质结构内在性质的GNM来分析 DTN 中的 PPI 层。
#这涉及将蛋白质相互作用组件想象成一个由相互连接的弹性珠子和弹簧精心构建的网络结构。
#全局网络连通性通过拉普拉斯矩阵来描绘，其中对角元素表示节点度数，非零
#的非对角元素表示节点之间的成对连接。在此背景下，PRS 计算作为一种使我们能够评估系统内节点如何对外部扰动响应的方法而出现。

# 依赖包
#if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
library(Matrix)

# 假设 PPIN_weighed 已经在环境中，列名为 gene1, gene2, combined_score
# 简单检查
stopifnot(all(c("gene1","gene2","combined_score") %in% colnames(PPIN_module_1)))

# 1) 去重（如果存在重复边，按平均权重合并；如果你偏好其它聚合方法，可改 FUN）
PPIN_edges <- aggregate(combined_score ~ gene1 + gene2, data = PPIN_module_1, FUN = mean)

# 2) 确保是无向图的一致表示：把 (u,v) 和 (v,u) 视为同一条边
#    统一边的顺序（小字母顺序或按因子顺序），防止 (A,B) 和 (B,A) 被当作不同条目
ord1 <- pmin(as.character(PPIN_edges$gene1), as.character(PPIN_edges$gene2))
ord2 <- pmax(as.character(PPIN_edges$gene1), as.character(PPIN_edges$gene2))
PPIN_edges$gA <- ord1
PPIN_edges$gB <- ord2
PPIN_edges2 <- aggregate(combined_score ~ gA + gB, data = PPIN_edges, FUN = mean)
colnames(PPIN_edges2) <- c("gene1","gene2","combined_score")

# 3) 构造节点列表（确定矩阵的行列顺序）
nodes <- sort(unique(c(PPIN_edges2$gene1, PPIN_edges2$gene2)))
n <- length(nodes)
message("节点数：", n, "；边数（无向，去重后）：", nrow(PPIN_edges2))

# 4) 索引映射
i_idx <- match(PPIN_edges2$gene1, nodes)
j_idx <- match(PPIN_edges2$gene2, nodes)
x_w   <- PPIN_edges2$combined_score

# 5) 构造对称的稀疏邻接矩阵 A （无自环）
#    把每条边放入 (i,j) 和 (j,i)
A <- sparseMatrix(i = c(i_idx, j_idx),
                  j = c(j_idx, i_idx),
                  x = c(x_w, x_w),
                  dims = c(n, n),
                  dimnames = list(nodes, nodes),
                  giveCsparse = TRUE)

# （可选）确保对角为0（无自环）
#diag(A) <- 0

# 6) 计算度向量与度矩阵 D（稀疏对角）
deg_vec <- rowSums(A)    # Matrix::rowSums 适用于稀疏矩阵
D <- Diagonal(x = deg_vec)

# 7) 计算拉普拉斯矩阵 L = D - A（依然是稀疏矩阵）
L <- D - A

# 8) 输出/检查
#print(L)                # 概览（稀疏矩阵信息）
message("非零元素数（L）：", length(L@x))
message("L 是否对称？ ", isSymmetric(L))

# 9) 保存（可选）
# saveRDS(L, file = "PPIN_laplacian_sparse.rds")
# writeMM(L, file = "PPIN_laplacian.mtx") # 如果需要 Matrix Market 格式导出

# 返回 L（如果在脚本/函数内使用，可以直接把 L 返回）
#L

#————————————————————————————————————————————————————————————————————————————
#对拉普拉斯矩阵性质的检查
# 1) 检查基本性质（稀疏矩阵下的高效检查）
message("---- 基本检查（稀疏形式） ----")
# 对称性（用范数衡量 L - t(L)）
sym_norm <- norm(L - Matrix::t(L), type = "F")
message("对称性（Frobenius范数 of L - t(L)）：", signif(sym_norm, 6))
if (sym_norm < 1e-8) message("=> L 近似对称（通过数值容差判断）") else message("=> 注意：L 可能不完全对称，请检查边/合并步骤或数值精度")

# 对角元素是否等于加权度
diag_diff_max <- max(abs(diag(L) - deg_vec))
message("对角与度差的最大绝对值：", signif(diag_diff_max, 6))
if (diag_diff_max < 1e-8) message("=> diag(L) 与 deg_vec 一致") else message("=> 注意：diag(L) 与 deg_vec 存在差异")

# 每行之和是否为0（理论上 Laplacian 行和为0）
row_sum_max <- max(abs(rowSums(L)))
message("行和的最大绝对值：", signif(row_sum_max, 6))
if (row_sum_max < 1e-8) message("=> 每行和接近 0（数值精度内）") else message("=> 注意：行和不为0，需检查 A 或 D 的构建")

# 非对角项是否非正（加权拉普拉斯 off-diagonal = -weight <= 0）
# 用稀疏数据直接检查：建立副本去掉对角，再看最大值
L_offdiag <- L
diag(L_offdiag) <- 0
max_offdiag <- max(L_offdiag)   # 如果存在正数，这会显示
message("非对角元素的最大值：", signif(max_offdiag, 6))
if (max_offdiag <= 1e-12) message("=> 非对角项均为非正（符合 L = D - A）") else message("=> 注意：存在正的非对角元素（不符合标准定义）")

# 2) 保存为 dense matrix（matrix 格式）并写入 RDS 文件
message("---- 转换为 dense matrix 并保存 ----")
# 安全转换：先估算内存
dense_bytes_est <- nrow(L) * ncol(L) * 8  # 8 bytes per double
message("估算保存为 dense matrix 需要内存（bytes）：", format(dense_bytes_est, scientific = FALSE))
# 尝试转换并保存，捕获可能的内存错误
saved <- FALSE
try({
  L_mat <- as.matrix(L)   # Matrix::as.matrix 对 1266x1266 是安全的
  # 额外检查尺寸
  stopifnot(dim(L_mat)[1] == dim(L_mat)[2])
  saveRDS(L_mat, file = "PPIN_laplacian_matrix_1266x1266.rds")
  message("已保存为 'PPIN_laplacian_matrix_1266x1266.rds' （RDS，二进制 R 矩阵）")
  saved <- TRUE
}, silent = FALSE)

if (!saved) {
  message("保存 dense 矩阵失败（内存或其他错误）。你可以改为保存稀疏格式：")
  message("  Matrix::writeMM(L, file='PPIN_laplacian_sparse.mtx') 或 saveRDS(L, file='PPIN_laplacian_sparse.rds')")
}

# 3) 近似检验：正半定性（PSD）——计算最小的几个特征值（稀疏方式，推荐 RSpectra）
message("---- PSD 检验（计算若干最小特征值） ----")
if (requireNamespace("RSpectra", quietly = TRUE)) {
  # 计算 6 个最小特征值（SM: smallest magnitude）
  # 注意：对于 Laplacian 最小特征值应接近 0（可能存在若干个零，取决于连通分量数）
  eigs_small <- try(RSpectra::eigs_sym(L, k = 6, which = "SM"), silent = TRUE)
  if (!inherits(eigs_small, "try-error")) {
    vals <- sort(Re(eigs_small$values))
    message("计算到的最小特征值（近似，按升序列出）：", paste(signif(vals,6), collapse = ", "))
    if (min(vals) >= -1e-8) message("=> 最小特征值近似 >= 0，L 似乎为正半定（数值容差内）")
    else message("=> 存在明显负特征值，需检查计算或数据")
  } else {
    message("使用 RSpectra 计算特征值失败：", eigs_small)
  }
} else {
  message("未安装包 'RSpectra'，若要做稀疏特征值计算请安装：install.packages('RSpectra').")
  message("备用（但可能较慢）：若已成功生成 L_mat，可以用 eigen(as.matrix(L_mat), symmetric=TRUE) 获取全部特征值。")
}

message("---- 检查完毕 ----")

#———————————————————————————————————————————————————————————————————————————————
#计算拉普拉斯矩阵的伪逆
# 计算拉普拉斯的 Moore-Penrose 伪逆（可选精确或近似）
# 输入：L (dgCMatrix) 已在环境中
# 输出：L_pinv (dense matrix) 并保存到 "L_pinv.rds"
# 并同时计算 sigma = diag(L_pinv) 并保存 "sigma.rds"

# 依赖
#if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
#if (!requireNamespace("RSpectra", quietly = TRUE)) install.packages("RSpectra")  # 仅近似模式需要
library(Matrix)

compute_pseudoinverse <- function(L,
                                  mode = c("auto", "exact", "approx"),
                                  n_threshold = 6000,   # 若 n <= n_threshold 则尝试 exact
                                  approx_k = 500,      # 近似模式下的特征个数（可调）
                                  tol_rel = 1e-8       # 特征值相对阈值（用于判断零模）
) {
  mode <- match.arg(mode)
  n <- nrow(L)
  if (is.null(n) || ncol(L) != n) stop("L 必须是方阵 dgCMatrix")
  message("矩阵尺寸：", n, " x ", n)
  
  # 决定模式
  if (mode == "auto") {
    mode_use <- if (n <= n_threshold) "exact" else "approx"
  } else {
    mode_use <- mode
  }
  message("计算模式：", mode_use)
  
  if (mode_use == "exact") {
    # 精确谱法（需要把 L 转成 dense，并做 eigen；内存/时间开销大）
    message("将 L 转为密集矩阵，开始做特征分解（可能耗内存/时间） ...")
    L_dense <- as.matrix(L)   # 注意：这会分配 O(n^2) 内存
    gc()
    eig <- eigen(L_dense, symmetric = TRUE)
    vals <- eig$values
    vecs <- eig$vectors
    # 设置阈值：相对最大特征值
    tol <- max(abs(vals)) * tol_rel
    message("特征值最大值：", max(vals), "， 相对阈值 tol = ", tol)
    keep_mask <- vals > tol
    num_keep <- sum(keep_mask)
    message("保留的非零特征个数：", num_keep, "（剔除 ", sum(!keep_mask), " 个小/零特征）")
    inv_vals <- numeric(length(vals))
    inv_vals[keep_mask] <- 1/vals[keep_mask]
    # 构造伪逆
    message("构造伪逆矩阵（dense） ...")
    L_pinv <- vecs %*% (inv_vals * t(vecs))   # 高效写法： vecs %*% diag(inv_vals) %*% t(vecs)
    # 清理
    rm(L_dense, eig); gc()
  } else {
    # 近似模式：用 RSpectra 求最小 k 个特征（which = "SM"）
    # 我们想得到最小的非零特征值/向量，跳过精确零模（若存在）
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      stop("近似模式需要安装 RSpectra 包：install.packages('RSpectra')")
    }
    k <- approx_k
    if (k >= n) {
      stop("approx_k 必须小于 n。若想做精确计算请把 mode 设置为 'exact'。")
    }
    message("近似模式：用 RSpectra 计算最小 ", k, " 个特征（可能需要一些时间） ...")
    # 如果 L 是稀疏，RSpectra 接受稀疏矩阵
    res <- RSpectra::eigs_sym(L, k = k, which = "SM")
    vals <- res$values
    vecs <- res$vectors
    # 忽略接近 0 的特征（零模）
    tol <- max(abs(vals)) * tol_rel
    keep_mask <- vals > tol
    if (sum(keep_mask) == 0) {
      stop("RSpectra 取得的最小特征中全部接近 0；请增大 approx_k 或使用 exact 模式。")
    }
    vals_k <- vals[keep_mask]
    vecs_k <- vecs[, keep_mask, drop = FALSE]
    message("在近似特征中保留 ", length(vals_k), " 个非零特征用于构造低秩伪逆（k_used=", length(vals_k), "）")
    inv_vals <- 1 / vals_k
    # 低秩近似的伪逆： V_k * diag(1/lambda_k) * V_k^T
    message("构造低秩近似的伪逆 ...")
    L_pinv <- vecs_k %*% (inv_vals * t(vecs_k))
    # 说明：L_pinv 是一个 dense 矩阵（n x n），但秩被限制为 k_used，可用作近似
  }
  
  # 校验对称性并返回
  if (!isSymmetric(L_pinv, tol = 1e-6)) {
    warning("计算得到的 L_pinv 非完全对称（数值误差），将进行对称化处理。")
    L_pinv <- (L_pinv + t(L_pinv)) / 2
  }
  return(L_pinv)
}

# -----------------------------
# 调用示例（默认 auto，会根据 n 选择 exact/approx）
# -----------------------------
# 假设 L 已经在全局环境中
L_pinv <- compute_pseudoinverse(L, mode = "exact", n_threshold = 6000, approx_k = 800, tol_rel = 1e-8)

# 保存结果（可能很大）
#saveRDS(L_pinv, file = "L_pinv.rds")
message("伪逆已保存为 L_pinv.rds")

# 可选：提取对角作为节点特异方差 sigma
sigma <- diag(L_pinv)
#saveRDS(sigma, file = "sigma_L_pinv_diag.rds")
message("sigma (diag of L_pinv) 已保存为 sigma_L_pinv_diag.rds")

# 赋予基因名称
rownames(L_pinv) <- rownames(L)
colnames(L_pinv) <- colnames(L)
names(sigma) <- rownames(L)
sigma_diag <- diag(sigma)

# 检查 sigma 中是否存在零值（求逆会导致无穷大）
if (any(sigma == 0)) {
  warning("sigma 中存在零值，求逆将产生无穷大，结果可能无效。")
}

# 构造 sigma_diag 的逆矩阵（对角矩阵，对角线元素为 1/sigma）
sigma_diag_inv <- diag(1 / sigma)   # 自动处理 Inf

# 矩阵乘法：L_pinv 左乘 sigma_diag_inv
PRS <- L_pinv %*% sigma_diag_inv

# 确保行名列名与 L_pinv 一致（理论上乘法后会自动保留，此处显式设置以防万一）
rownames(PRS) <- rownames(L_pinv)
colnames(PRS) <- colnames(L_pinv)


#———————————————————————————————————————————————————————————————————————————————
#计算扰动结果矩阵与突变扰动分数

# 假设 mutation_score 和 L_pinv 已存在于环境中

# 1. 提取前五列并调整列顺序为：Disease, Gene, Site, Mutation, total_score
# 假设原始列名为 "Gene", "Disease", "Site", "Mutation", "total_score"
mutation_data <- mutation_score[, 1:5]
# 按指定顺序重排列
mutation_data <- mutation_data[, c("Disease", "Gene", "Site", "Mutation", "total_score")]

# 2. 多列排序：Disease -> Gene -> Site -> Mutation -> total_score (全部升序)
mutation_data <- mutation_data[order(mutation_data$Disease,
                                     mutation_data$Gene,
                                     mutation_data$Site,
                                     mutation_data$Mutation,
                                     mutation_data$total_score), ]
# 确保 Gene 列为字符型（如果是因子则转换）
if (is.factor(mutation_data$Gene)) {
  mutation_data$Gene <- as.character(mutation_data$Gene)
}

# 替换基因名
mutation_data$Gene[mutation_data$Gene == "MEK1"] <- "MAP2K1"
mutation_data$Gene[mutation_data$Gene == "MEK2"] <- "MAP2K2"
# 3. 准备 PRS 的行名（基因名）作为结果数据框的行名
gene_names <- rownames(PRS)          # 或 colnames(PRS)，两者一致
n_genes <- length(gene_names)
n_muts  <- nrow(mutation_data)

# 4. 预分配结果矩阵（1266行 × 4267列）
result_matrix <- matrix(0, nrow = n_genes, ncol = n_muts)
rownames(result_matrix) <- gene_names

# 为列设置描述性名称（使用疾病、基因、位点和突变组合，并确保唯一性）
col_names <- with(mutation_data, paste(Disease, Gene, Site, Mutation, sep = "_"))
col_names <- make.unique(col_names, sep = "_")   # 若重复则自动添加后缀
colnames(result_matrix) <- col_names

# 5. 遍历每一行，计算并填充矩阵
for (i in seq_len(n_muts)) {
  score <- mutation_data$total_score[i]
  gene  <- mutation_data$Gene[i]
  
  # 取出 PRS 中对应基因的列（向量），乘以 score，存入结果矩阵
  result_matrix[, i] <- PRS[, gene] * score
}

# 6. 转换为数据框（便于后续分析）
result_df <- as.data.frame(result_matrix)

# 可选：保存结果
# saveRDS(result_df, file = "mutation_effect_columns.rds")


result_df_origin<-result_df

result_df<-abs(result_df)










