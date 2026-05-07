# 加载必要的包
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(car)       # leveneTest
library(stats)     # 基础统计函数

# 假设 result_df 已存在，行名是蛋白标识，列名包含疾病类型
# 示例：result_df <- read.csv("your_file.csv", row.names=1)

# 1. 按疾病分组列 ------------------------------------------------------------
# 从列名中提取疾病类型（第一个下划线前的部分）
disease_labels <- sapply(strsplit(colnames(result_df), "_"), `[`, 1)
# 检查是否只有 Cancer, NDD, OR
unique(disease_labels)

# 获取各组的列索引
cancer_cols <- which(disease_labels == "Cancer")
ndd_cols    <- which(disease_labels == "NDD")
or_cols     <- which(disease_labels == "OR")

# 检查各组突变数
cat("Cancer mutations:", length(cancer_cols), "\n")
cat("NDD mutations:", length(ndd_cols), "\n")
cat("OR mutations:", length(or_cols), "\n")

# 2. 定义函数：对单个蛋白进行两组比较，返回结果 ---------------------------------
# 输入：蛋白对应的数值向量（已按组分好），组名，样本数
# 输出：包含均值、log2FC、p值等的数据框一行
compare_two_groups <- function(values_group1, values_group2, name_group1, name_group2) {
  # 计算均值
  mean1 <- mean(values_group1, na.rm = TRUE)
  mean2 <- mean(values_group2, na.rm = TRUE)
  
  # 计算log2FC（假设值非负，若出现0则加一个小常数）
  if (all(values_group1 >= 0, na.rm = TRUE) & all(values_group2 >= 0, na.rm = TRUE)) {
    fc <- mean1 / mean2
    if (fc == 0) fc <- 1e-10  # 避免log2(0)
    log2fc <- log2(fc)
  } else {
    # 若有负值，改用均值差
    log2fc <- mean1 - mean2
    warning("Negative values detected, using mean difference instead of log2FC.")
  }
  
  # 检查样本量是否足够进行正态性检验
  n1 <- length(values_group1)
  n2 <- length(values_group2)
  
  # 默认使用Wilcoxon
  use_t_test <- FALSE
  
  # 尝试进行正态性和方差齐性检验
  if (n1 >= 3 & n2 >= 3) {
    # 正态性检验 (Shapiro-Wilk)
    norm1 <- tryCatch(shapiro.test(values_group1)$p.value > 0.05, error = function(e) FALSE)
    norm2 <- tryCatch(shapiro.test(values_group2)$p.value > 0.05, error = function(e) FALSE)
    
    # 方差齐性检验 (Levene检验，基于中位数)
    combined <- c(values_group1, values_group2)
    groups <- factor(rep(c(name_group1, name_group2), times = c(n1, n2)))
    var_equal <- tryCatch(leveneTest(combined, groups)$`Pr(>F)`[1] > 0.05, error = function(e) FALSE)
    
    if (norm1 & norm2 & var_equal) {
      use_t_test <- TRUE
    }
  }
  
  # 执行检验
  if (use_t_test) {
    test_result <- t.test(values_group1, values_group2, var.equal = TRUE)
    p_value <- test_result$p.value
    method <- "t-test"
  } else {
    test_result <- wilcox.test(values_group1, values_group2, exact = FALSE)
    p_value <- test_result$p.value
    method <- "Wilcoxon"
  }
  
  # 返回结果
  data.frame(
    protein = rownames(result_df)[i],  # 将在循环中传入i
    group1 = name_group1,
    group2 = name_group2,
    mean1 = mean1,
    mean2 = mean2,
    log2fc = log2fc,
    p_value = p_value,
    method = method,
    stringsAsFactors = FALSE
  )
}

# 3. 循环所有蛋白，进行三组两两比较 --------------------------------------------
# 初始化空数据框存储结果
comparisons <- c("Cancer_vs_NDD", "Cancer_vs_OR", "NDD_vs_OR")
results <- list()

# 获取蛋白总数
n_proteins <- nrow(result_df)

# 预分配列表以加快速度
res_cv_ndd <- vector("list", n_proteins)
res_cv_or  <- vector("list", n_proteins)
res_ndd_or <- vector("list", n_proteins)

# 将数据框转为矩阵以便快速索引（可选，但保留行名）
result_mat <- as.matrix(result_df)

# 循环每个蛋白
for (i in 1:n_proteins) {
  # 提取该蛋白在各组的值
  vals_cancer <- result_mat[i, cancer_cols]
  vals_ndd    <- result_mat[i, ndd_cols]
  vals_or     <- result_mat[i, or_cols]
  
  # 比较 Cancer vs NDD
  res_cv_ndd[[i]] <- compare_two_groups(vals_cancer, vals_ndd, "Cancer", "NDD")
  res_cv_ndd[[i]]$protein <- rownames(result_df)[i]  # 覆盖蛋白名
  
  # 比较 Cancer vs OR
  res_cv_or[[i]] <- compare_two_groups(vals_cancer, vals_or, "Cancer", "OR")
  res_cv_or[[i]]$protein <- rownames(result_df)[i]
  
  # 比较 NDD vs OR
  res_ndd_or[[i]] <- compare_two_groups(vals_ndd, vals_or, "NDD", "OR")
  res_ndd_or[[i]]$protein <- rownames(result_df)[i]
}

# 合并为数据框
df_cv_ndd <- do.call(rbind, res_cv_ndd)
df_cv_or  <- do.call(rbind, res_cv_or)
df_ndd_or <- do.call(rbind, res_ndd_or)

# 4. 对每个比较进行FDR校正 ----------------------------------------------------
df_cv_ndd$padj <- p.adjust(df_cv_ndd$p_value, method = "BH")
df_cv_or$padj  <- p.adjust(df_cv_or$p_value, method = "BH")
df_ndd_or$padj <- p.adjust(df_ndd_or$p_value, method = "BH")

# 添加显著性标签（用于火山图）
threshold_fc <- 0.5   # log2FC阈值，可调整
threshold_p  <- 0.05

df_cv_ndd$sig <- with(df_cv_ndd, ifelse(padj < threshold_p & abs(log2fc) > threshold_fc, 
                                        ifelse(log2fc > 0, "Up in Cancer", "Down in Cancer"), "Not Sig"))
df_cv_or$sig  <- with(df_cv_or, ifelse(padj < threshold_p & abs(log2fc) > threshold_fc,
                                       ifelse(log2fc > 0, "Up in Cancer", "Down in Cancer"), "Not Sig"))
df_ndd_or$sig <- with(df_ndd_or, ifelse(padj < threshold_p & abs(log2fc) > threshold_fc,
                                        ifelse(log2fc > 0, "Up in NDD", "Down in NDD"), "Not Sig"))

# 5. 绘制火山图 --------------------------------------------------------------
# 自定义火山图函数
plot_volcano <- function(df, title) {
  ggplot(df, aes(x = log2fc, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up in Cancer" = "red", "Down in Cancer" = "blue",
                                  "Up in NDD" = "red", "Down in NDD" = "blue",
                                  "Not Sig" = "gray")) +
    geom_vline(xintercept = c(-threshold_fc, threshold_fc), linetype = "dashed", color = "darkgray") +
    geom_hline(yintercept = -log10(threshold_p), linetype = "dashed", color = "darkgray") +
    labs(title = title, x = "log2(Fold Change)", y = "-log10(adjusted p-value)") +
    theme_minimal() +
    theme(legend.title = element_blank())
}

# 绘制并保存（可调整）
p1 <- plot_volcano(df_cv_ndd, "Cancer vs NDD")
p2 <- plot_volcano(df_cv_or,  "Cancer vs OR")
p3 <- plot_volcano(df_ndd_or, "NDD vs OR")

print(p1)
print(p2)
print(p3)

# 保存为文件（可选）
# ggsave("volcano_CvN.pdf", p1, width = 8, height = 6)
# ggsave("volcano_CvO.pdf", p2, width = 8, height = 6)
# ggsave("volcano_NvO.pdf", p3, width = 8, height = 6)

# 6. 提取显著蛋白用于热图 ----------------------------------------------------
# 定义显著阈值（可调整）
sig_padj <- 0.05
sig_log2fc <- 0.5

# 从每个比较中提取显著蛋白（上调或下调）
sig_cv_ndd <- df_cv_ndd$protein[df_cv_ndd$padj < sig_padj & abs(df_cv_ndd$log2fc) > sig_log2fc]
sig_cv_or  <- df_cv_or$protein[df_cv_or$padj < sig_padj & abs(df_cv_or$log2fc) > sig_log2fc]
sig_ndd_or <- df_ndd_or$protein[df_ndd_or$padj < sig_padj & abs(df_ndd_or$log2fc) > sig_log2fc]

# 取并集
sig_proteins <- unique(c(sig_cv_ndd, sig_cv_or, sig_ndd_or))
cat("Number of significant proteins:", length(sig_proteins), "\n")

# 如果没有显著蛋白，提示并退出
if (length(sig_proteins) == 0) {
  stop("No significant proteins found. Try relaxing thresholds.")
}

# 7. 计算这些蛋白在各疾病组的平均扰动得分 --------------------------------------
# 提取显著蛋白的行数据
sig_data <- result_df[sig_proteins, ]

# 按疾病组计算每行的均值
mean_cancer <- rowMeans(sig_data[, cancer_cols], na.rm = TRUE)
mean_ndd    <- rowMeans(sig_data[, ndd_cols], na.rm = TRUE)
mean_or     <- rowMeans(sig_data[, or_cols], na.rm = TRUE)

# 构建热图矩阵（行：蛋白，列：疾病组）
heatmap_mat <- cbind(Cancer = mean_cancer, NDD = mean_ndd, OR = mean_or)
rownames(heatmap_mat) <- sig_proteins

# 可选：对每行进行z-score标准化（使不同蛋白可比）
heatmap_mat_scaled <- t(scale(t(heatmap_mat)))  # 按行标准化

# 8. 绘制热图 -----------------------------------------------------------------
# 使用pheatmap，可添加注释等
pheatmap(heatmap_mat_scaled,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",  # 已手动标准化
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Average perturbation scores of significant proteins",
         fontsize_row = 8,
         display_numbers = FALSE)

# 若需保存，可指定文件名
# pheatmap(heatmap_mat_scaled, filename = "heatmap_sig_proteins.pdf", ...)

# 9. 可选：输出差异分析结果到CSV ----------------------------------------------
write.csv(df_cv_ndd, "C:\\毕业论文\\数据\\diff_Cancer_vs_NDD.csv", row.names = FALSE)
write.csv(df_cv_or,  "C:\\毕业论文\\数据\\diff_Cancer_vs_OR.csv", row.names = FALSE)
write.csv(df_ndd_or, "C:\\毕业论文\\数据\\diff_NDD_vs_OR.csv", row.names = FALSE)


#富集分析
# -------------------------
library(clusterProfiler)
library(org.Hs.eg.db)   # 如果不是人类请换相应 OrgDb
library(enrichplot)
library(ggplot2)
library(dplyr)

# -------------------------
# 通用函数：对单个数据框做富集并画图
# -------------------------
run_go_and_plot <- function(df, sig_label, out_prefix, species_orgdb = org.Hs.eg.db, 
                            symbol_col = "protein", sig_col = "sig",
                            p_cutoff = 0.05, top_n = 20) {
  message("=== 处理：", out_prefix, " | 筛选 sig = ", sig_label, " ===")
  # 1) 提取基因（symbol）
  genes_sym <- df %>% filter(.data[[sig_col]] == sig_label) %>% pull(.data[[symbol_col]]) %>% unique() %>% na.omit()
  bg_sym <- df %>% pull(.data[[symbol_col]]) %>% unique() %>% na.omit()
  
  if (length(genes_sym) == 0) {
    message("未找到满足条件的基因：", out_prefix, " - 跳过。")
    return(NULL)
  }
  
  # 2) SYMBOL -> ENTREZID 映射
  gene_map <- bitr(genes_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_orgdb)
  bg_map   <- bitr(bg_sym,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = species_orgdb)
  if (nrow(gene_map) == 0) {
    message("筛选出的基因没有成功映射到 ENTREZID（请检查是否为人类 SYMBOL 或更换 OrgDb）。")
    return(NULL)
  }
  genes_entrez <- unique(gene_map$ENTREZID)
  universe_entrez <- unique(bg_map$ENTREZID)
  
  # 3) 对三个本体做富集并绘图
  ontologies <- c("MF", "BP", "CC")
  res_list <- list()
  for (ont in ontologies) {
    message("  -> 进行富集：", ont)
    ego <- tryCatch({
      enrichGO(gene         = genes_entrez,
               universe     = universe_entrez,
               OrgDb        = species_orgdb,
               keyType      = "ENTREZID",
               ont          = ont,
               pvalueCutoff = p_cutoff,
               qvalueCutoff = p_cutoff,
               readable     = TRUE)
    }, error = function(e) {
      message("    enrichGO 错误：", e$message)
      return(NULL)
    })
    
    if (is.null(ego) || (is.data.frame(as.data.frame(ego)) && nrow(as.data.frame(ego)) == 0)) {
      message("    没有显著富集项（p < ", p_cutoff, "），跳过绘图。")
      res_list[[ont]] <- NULL
      next
    }
    
    # 取前 top_n 项（按 p.adjust 或 pvalue 排序）
    ego_df <- as.data.frame(ego)
    show_k <- min(top_n, nrow(ego_df))
    
    # dotplot（气泡图）
    p <- dotplot(ego, showCategory = show_k, orderBy = "GeneRatio") +
      ggtitle(paste0(out_prefix, " | GO-", ont, " | sig=", sig_label)) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5))
    print(p)
    # 保存图片
    fname <- paste0(out_prefix, "_GO", ont, "_sig-", gsub(" ", "_", sig_label), ".png")
    ggsave(filename = fname, plot = p, width = 8, height = 6, dpi = 300)
    message("    已保存：", fname)
    
    # 也保存结果表到变量
    res_list[[ont]] <- ego
  }
  
  return(res_list)
}

# -------------------------
# 对你指定的三组分别运行
# -------------------------
# df_cv_ndd: 筛选 sig == "Down in Cancer"
res_cv_ndd <- run_go_and_plot(df = df_cv_ndd,
                              sig_label = "Down in Cancer",
                              out_prefix = "df_cv_ndd")

# df_cv_or: 筛选 sig == "Down in Cancer"
res_cv_or  <- run_go_and_plot(df = df_cv_or,
                              sig_label = "Down in Cancer",
                              out_prefix = "df_cv_or")

# df_ndd_or: 筛选 sig == "Down in NDD"
res_ndd_or <- run_go_and_plot(df = df_ndd_or,
                              sig_label = "Down in NDD",
                              out_prefix = "df_ndd_or")

# 现在 res_cv_ndd / res_cv_or / res_ndd_or 各自是一个 list(MF=..., BP=..., CC=...)
# 可以查看某个富集结果，例如：
# as.data.frame(res_cv_ndd$BP)