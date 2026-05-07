# =========================================================
# 节点筛选：只做三组两两比较，不再做Kruskal-Wallis
# =========================================================

library(dplyr)
library(ggplot2)

out_dir <- "C:/毕业论文/数据处理及绘图/ppi_build/python_prs_output/disease_difference"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05   # 显著性阈值，按FDR理解
eps <- 1e-8     # 防止log2FC计算时出现0除错

# ---------------------------------------------------------
# 1. 提取疾病分组信息
# ---------------------------------------------------------
mut_names <- colnames(prs_raw)

mut_info <- data.frame(
  mutation = mut_names,
  disease = stringr::str_match(mut_names, "^(Cancer|NDD|OR)_")[, 2],
  protein = stringr::str_match(mut_names, "^(Cancer|NDD|OR)_([^_]+)_")[, 3],
  stringsAsFactors = FALSE
)

if (any(is.na(mut_info$disease))) {
  stop("部分列名未能解析出 Disease type，请检查列名格式是否为 Cancer/NDD/OR_蛋白_位点_变体。")
}

mut_info$disease <- factor(mut_info$disease, levels = c("Cancer", "NDD", "OR"))

idx_C <- which(mut_info$disease == "Cancer")
idx_N <- which(mut_info$disease == "NDD")
idx_O <- which(mut_info$disease == "OR")

# ---------------------------------------------------------
# 2. 定义两两比较函数
#    每一行 = 一个节点
#    每个节点比较两个疾病组中的全部突变列
# ---------------------------------------------------------
pairwise_node_test <- function(mat, idx1, idx2, g1, g2, eps = 1e-6) {
  med1 <- apply(mat[, idx1, drop = FALSE], 1, median, na.rm = TRUE)
  med2 <- apply(mat[, idx2, drop = FALSE], 1, median, na.rm = TRUE)
  
  # 这里用 Wilcoxon rank-sum test（非配对）
  pval <- vapply(seq_len(nrow(mat)), function(i) {
    suppressWarnings(
      wilcox.test(
        mat[i, idx1],
        mat[i, idx2],
        exact = FALSE
      )$p.value
    )
  }, numeric(1))
  
  fdr <- p.adjust(pval, method = "BH")
  
  # log2FC 用中位数计算，适合非正态、偏态的扰动值
  log2fc <- log2((med1 + eps) / (med2 + eps))
  
  data.frame(
    node = rownames(mat),
    group1 = g1,
    group2 = g2,
    med1 = med1,
    med2 = med2,
    log2FC = log2fc,
    diff_median = med1 - med2,
    p_value = pval,
    fdr = fdr,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------
# 3. 运行三组两两比较
# ---------------------------------------------------------
mat <- as.matrix(prs_raw)

res_CN <- pairwise_node_test(mat, idx_C, idx_N, "Cancer", "NDD", eps = eps)
res_CO <- pairwise_node_test(mat, idx_C, idx_O, "Cancer", "OR",  eps = eps)
res_NO <- pairwise_node_test(mat, idx_N, idx_O, "NDD",    "OR",  eps = eps)

write.csv(res_CN, file.path(out_dir, "Cancer_vs_NDD_node_results.csv"), row.names = FALSE)
write.csv(res_CO, file.path(out_dir, "Cancer_vs_OR_node_results.csv"), row.names = FALSE)
write.csv(res_NO, file.path(out_dir, "NDD_vs_OR_node_results.csv"), row.names = FALSE)

# ---------------------------------------------------------
# 4. 火山图函数
#    x轴：log2FC
#    y轴：-log10(FDR)
# ---------------------------------------------------------
make_volcano <- function(df, title_text, group1, group2, alpha = 0.05) {
  df <- df %>%
    mutate(
      direction = case_when(
        fdr < alpha & log2FC > 0.5 ~ paste0(group1, " higher"),
        fdr < alpha & log2FC < -0.5 ~ paste0(group2, " higher"),
        TRUE ~ "Not significant"
      ),
      neglog10_fdr = -log10(fdr + 1e-300)
    )
  
  cols <- c("#D73027", "#4575B4", "grey70")
  names(cols) <- c(
    paste0(group1, " higher"),
    paste0(group2, " higher"),
    "Not significant"
  )
  
  ggplot(df, aes(x = log2FC, y = neglog10_fdr)) +
    geom_point(aes(color = direction), size = 1.3, alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.6) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, linewidth = 0.6, color = "grey50") +
    geom_hline(yintercept = -log10(alpha), linetype = 2, linewidth = 0.6) +
    scale_color_manual(values = cols) +
    theme_classic(base_size = 12) +
    labs(
      title = title_text,
      x = "log2FC_(median difference)",
      y = expression(-log[10](FDR)),
      color = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

p_CN <- make_volcano(res_CN, "Cancer vs NDD", "Cancer", "NDD", alpha)
p_CO <- make_volcano(res_CO, "Cancer vs OR", "Cancer", "OR", alpha)
p_NO <- make_volcano(res_NO, "NDD vs OR", "NDD", "OR", alpha)

ggsave(file.path(out_dir, "volcano_Cancer_vs_NDD.png"), p_CN, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "volcano_Cancer_vs_NDD.pdf"), p_CN, width = 8, height = 6)

ggsave(file.path(out_dir, "volcano_Cancer_vs_OR.png"), p_CO, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "volcano_Cancer_vs_OR.pdf"), p_CO, width = 8, height = 6)

ggsave(file.path(out_dir, "volcano_NDD_vs_OR.png"), p_NO, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "volcano_NDD_vs_OR.pdf"), p_NO, width = 8, height = 6)

# ---------------------------------------------------------
# 5. 提取疾病特异扰动节点
#    定义：
#    - Cancer特异：Cancer在 Cancer vs NDD 和 Cancer vs OR 中都显著更高
#    - NDD特异：NDD在 Cancer vs NDD 和 NDD vs OR 中都显著更高
#    - OR特异：OR在 Cancer vs OR 和 NDD vs OR 中都显著更高
# ---------------------------------------------------------
CN_Cancer_high <- res_CN %>% filter(fdr < alpha, log2FC > 0.5) %>% pull(node)
CN_NDD_high    <- res_CN %>% filter(fdr < alpha, log2FC < -0.5) %>% pull(node)

CO_Cancer_high <- res_CO %>% filter(fdr < alpha, log2FC > 0.5) %>% pull(node)
CO_OR_high     <- res_CO %>% filter(fdr < alpha, log2FC < -0.5) %>% pull(node)

NO_NDD_high    <- res_NO %>% filter(fdr < alpha, log2FC > 0.5) %>% pull(node)
NO_OR_high     <- res_NO %>% filter(fdr < alpha, log2FC < -0.5) %>% pull(node)

Cancer_specific <- Reduce(intersect, list(CN_Cancer_high, CO_Cancer_high))
NDD_specific    <- Reduce(intersect, list(CN_NDD_high, NO_NDD_high))
OR_specific     <- Reduce(intersect, list(CO_OR_high, NO_OR_high))

write.table(Cancer_specific, file.path(out_dir, "Cancer_specific_nodes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(NDD_specific, file.path(out_dir, "NDD_specific_nodes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(OR_specific, file.path(out_dir, "OR_specific_nodes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(CN_Cancer_high, file.path(out_dir, "CN_Cancer_high.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(CN_NDD_high, file.path(out_dir, "CN_NDD_high.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(CO_Cancer_high, file.path(out_dir, "CO_Cancer_high.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(CO_OR_high, file.path(out_dir, "CO_OR_high.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(NO_NDD_high, file.path(out_dir, "NO_NDD_high.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(NO_OR_high, file.path(out_dir, "NO_OR_high.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---------------------------------------------------------
# 6. 两组共享的无显著性差异节点
#    例如 Cancer vs NDD 中不显著
# ---------------------------------------------------------
CN_shared_stable <- res_CN %>% filter(fdr >= alpha) %>% pull(node)
CO_shared_stable <- res_CO %>% filter(fdr >= alpha) %>% pull(node)
NO_shared_stable <- res_NO %>% filter(fdr >= alpha) %>% pull(node)

write.table(CN_shared_stable, file.path(out_dir, "shared_stable_Cancer_vs_NDD.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(CO_shared_stable, file.path(out_dir, "shared_stable_Cancer_vs_OR.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(NO_shared_stable, file.path(out_dir, "shared_stable_NDD_vs_OR.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---------------------------------------------------------
# 7. 三组共同稳定节点
#    三个两两比较都不显著
# ---------------------------------------------------------
common_stable <- Reduce(intersect, list(
  CN_shared_stable,
  CO_shared_stable,
  NO_shared_stable
))

write.table(common_stable, file.path(out_dir, "common_stable_nodes.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---------------------------------------------------------
# 8. 汇总输出
# ---------------------------------------------------------
summary_df <- data.frame(
  category = c(
    "Cancer-specific nodes",
    "NDD-specific nodes",
    "OR-specific nodes",
    "Cancer vs NDD shared stable",
    "Cancer vs OR shared stable",
    "NDD vs OR shared stable",
    "Common stable nodes"
  ),
  count = c(
    length(Cancer_specific),
    length(NDD_specific),
    length(OR_specific),
    length(CN_shared_stable),
    length(CO_shared_stable),
    length(NO_shared_stable),
    length(common_stable)
  )
)

write.csv(summary_df, file.path(out_dir, "summary_counts.csv"), row.names = FALSE)

cat("节点筛选与火山图分析完成！结果已输出到：\n", out_dir, "\n")

############################################################################
###富集分析###
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# =========================
# 1. 生成 OR_high
# =========================
OR_high <- intersect(CO_OR_high, NO_OR_high)

# 去掉NA和重复
CN_NDD_high   <- unique(na.omit(CN_NDD_high))
CN_Cancer_high <- unique(na.omit(CN_Cancer_high))
OR_high       <- unique(na.omit(OR_high))

# =========================
# 2. 设置输出目录
# =========================
out_dir <- "C:/毕业论文/数据处理及绘图/ppi_build/python_prs_output/disease_difference/GO_enrichment"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# 3. 定义富集 + 作图函数
# =========================
run_go_and_plot <- function(gene_vec, group_name, ont, out_dir, showCategory = 20) {
  # GO富集
  ego <- enrichGO(
    gene          = gene_vec,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,          # "BP" / "MF" / "CC"
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  # 保存富集结果表
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    message(group_name, " - ", ont, ": no enriched terms.")
    return(NULL)
  }
  
  res_file <- file.path(out_dir, paste0(group_name, "_GO_", ont, "_enrichment.csv"))
  write.csv(as.data.frame(ego), res_file, row.names = FALSE)
  
  # 气泡图
  # 气泡图
  p <- dotplot(ego, showCategory = showCategory) +
    ggtitle(paste0(group_name, " GO_", ont, " enrichment")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # 保存图片
  png_file <- file.path(out_dir, paste0(group_name, "_GO_", ont, "_dotplot.png"))
  pdf_file <- file.path(out_dir, paste0(group_name, "_GO_", ont, "_dotplot.pdf"))
  
  ggsave(png_file, plot = p, width = 9, height = 10, dpi = 300)
  ggsave(pdf_file, plot = p, width = 9, height = 10)
  
  return(ego)
}

# =========================
# 4. 分别对三组做 GO_MF / GO_BP / GO_CC
# =========================
groups <- list(
  CN_NDD_high    = CN_NDD_high,
  CN_Cancer_high = CN_Cancer_high,
  OR_high        = OR_high
)

onts <- c("MF", "BP", "CC")

# 保存结果对象
go_results <- list()

for (gname in names(groups)) {
  gene_vec <- groups[[gname]]
  
  for (ont in onts) {
    key_name <- paste0(gname, "_GO_", ont)
    go_results[[key_name]] <- run_go_and_plot(
      gene_vec   = gene_vec,
      group_name = gname,
      ont        = ont,
      out_dir    = out_dir,
      showCategory = 20
    )
  }
}

cat("GO富集分析完成，结果已保存到：\n", out_dir, "\n")











