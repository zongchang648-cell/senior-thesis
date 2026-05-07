
library(rrvgo)
library(GOSemSim)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(wordcloud)
library(RColorBrewer)

# 通用处理函数：输入富集表（含 ID, Description, qvalue 或 p.adjust, Count, geneID 等），和 ontology ('BP','MF','CC')
process_rrvgo_module <- function(enrich_df, ont = c("BP","MF","CC"), orgdb = "org.Hs.eg.db", max_clusters = 20) {
  ont <- match.arg(ont)
  # 检查必须列
  stopifnot("ID" %in% colnames(enrich_df))
  if(!("Description" %in% colnames(enrich_df))) enrich_df$Description <- NA_character_
  
  go_ids <- unique(as.character(enrich_df$ID))
  if(length(go_ids) < 2) stop("GO term 数量 < 2，无法计算语义相似性。")
  
  message("计算语义相似性矩阵（this can take some seconds/minutes depending on number of terms)...")
  # 使用 rrvgo 的 calculateSimMatrix（内部调用 GOSemSim）
  simMatrix <- calculateSimMatrix(go_ids, orgdb = orgdb, ont = ont, method = "Rel")
  
  # 构建 score：优先 qvalue（若存在），否则 p.adjust；用 -log10 转换（越大越好）
  if("qvalue" %in% colnames(enrich_df)) {
    sco_raw <- enrich_df$qvalue
  } else if("p.adjust" %in% colnames(enrich_df)) {
    sco_raw <- enrich_df$p.adjust
  } else if("pvalue" %in% colnames(enrich_df)) {
    sco_raw <- enrich_df$pvalue
  } else {
    sco_raw <- rep(NA_real_, nrow(enrich_df))
  }
  scores <- -log10(sco_raw)
  names(scores) <- enrich_df$ID
  # 保留只有 simMatrix 包含的那些 IDs 的 scores，并按照 simMatrix 的行列顺序设置
  common_ids <- intersect(rownames(simMatrix), names(scores))
  scores <- scores[common_ids]
  # 若 scores 都为 NA，则设置 NULL（让 reduceSimMatrix 使用 term size/uniqueness）
  if(all(is.na(scores))) scores2 <- NULL else scores2 <- scores
  
  # 自动搜索 threshold（从高到低），以得到 <= max_clusters 个簇
  thresholds <- seq(0.95, 0.40, by = -0.01)
  chosen_th <- NA_real_
  reduced <- NULL
  for(th in thresholds){
    rt <- tryCatch({
      reduceSimMatrix(simMatrix, scores2, threshold = th, orgdb = orgdb)
    }, error = function(e) {
      warning("reduceSimMatrix 出错（threshold=", th, "）：", conditionMessage(e))
      NULL
    })
    if(is.null(rt)) next
    # reduced 返回的是每个原始 term 的 mapping（列名通常含 go, term, parent, parentTerm, score, size, uniq, disp）
    if("parent" %in% colnames(rt)) {
      nclus <- length(unique(rt$parent))
    } else if("go" %in% colnames(rt) && "parent" %in% colnames(rt)) {
      nclus <- length(unique(rt$parent))
    } else {
      # fallback: 以 parentTerm 计数
      nclus <- length(unique(rt$parentTerm))
    }
    if(nclus <= max_clusters) {
      chosen_th <- th
      reduced <- rt
      message(sprintf("选定 threshold = %.2f -> %d clusters", th, nclus))
      break
    }
  }
  
  # 如果循环中没有找到合适 threshold（簇数仍 > max_clusters），做后备处理：
  # 用较低 threshold 做一次 reduce，然后取得分最高的 max_clusters 个代表 parent，
  # 把其它 terms 映射到这 20 个代表里相似度最高者
  if(is.null(reduced)) {
    message("未能在搜索范围内得到 <= ", max_clusters, " 簇，退回到 threshold = 0.40 并做 top-N 代表映射。")
    reduced0 <- reduceSimMatrix(simMatrix, scores2, threshold = 0.40, orgdb = orgdb)
    # 先确定所有 parent（簇代表）
    parents_all <- unique(reduced0$parent)
    # 计算每个 parent 的代表性分数（取 reduced0 中每个 parent 的最大 score）
    if("score" %in% colnames(reduced0)) {
      parent_scores <- tapply(reduced0$score, reduced0$parent, max, na.rm = TRUE)
    } else {
      # 若没有 score（极少），尽量用 term size拟代
      parent_scores <- setNames(rep(1, length(parents_all)), parents_all)
    }
    top_parents <- names(sort(parent_scores, decreasing = TRUE))[seq_len(min(length(parent_scores), max_clusters))]
    # 对每个 go term，找到与 top_parents 中 similarity 最大的 parent（用 simMatrix）
    sim_sub <- simMatrix
    # 保证有行列
    all_terms <- rownames(simMatrix)
    assign_parent_by_sim <- sapply(all_terms, function(g){
      # 对 top_parents 找相似度（若 g 本身是 parent 之一则可能为 1）
      sims <- sim_sub[g, top_parents, drop = TRUE]
      top <- top_parents[which.max(sims)]
      return(top)
    }, USE.NAMES = TRUE)
    # 构造 reduced dataframe with mapping (one row per term)
    reduced <- data.frame(
      go = all_terms,
      term = ifelse(all_terms %in% enrich_df$ID,
                    enrich_df$Description[match(all_terms, enrich_df$ID)],
                    NA_character_),
      parent = as.character(assign_parent_by_sim[all_terms]),
      parentTerm = NA_character_,
      score = scores[all_terms],
      stringsAsFactors = FALSE
    )
    # fetch parentTerm descriptions
    reduced$parentTerm <- getGoTerm(reduced$parent)
  }
  
  # 现在 reduced 应该是一个 data.frame，列含 go, term, parent, parentTerm, score
  # 为安全起见，标准化列名到我们需要的字段
  if(!"go" %in% colnames(reduced)) {
    if("GO" %in% colnames(reduced)) reduced$go <- reduced$GO
    else stop("reduceSimMatrix 返回结果中找不到 'go' 列。")
  }
  if(!"term" %in% colnames(reduced)) {
    if("term" %in% colnames(reduced)) NULL else reduced$term <- getGoTerm(reduced$go)
  }
  if(!"parent" %in% colnames(reduced)) reduced$parent <- reduced$go  # 兜底
  
  # 为每个 parent 分配一个 numeric cluster id（按 parent 排序）
  parent_levels <- unique(reduced$parent)
  cluster_id_map <- setNames(seq_along(parent_levels), parent_levels)
  reduced$cluster_id <- cluster_id_map[as.character(reduced$parent)]
  
  # 构造输出 mapping dataframe（每个原始条目属于哪个聚类）
  mapping_df <- reduced %>%
    dplyr::select(ID = go, Description = term, parent = parent, parentTerm = parentTerm, score = score, cluster_id) %>%
    mutate(ID = as.character(ID), parent = as.character(parent), parentTerm = as.character(parentTerm))
  
  # 如果原始 enrich_df 中有更多列（Count, geneID 等），把它们合并回 mapping_df（按 ID）
  mapping_df <- mapping_df %>%
    left_join(enrich_df, by = c("ID" = "ID")) %>%
    # 保持必要列顺序（把 cluster info 前置）
    dplyr::select(ID, Description = Description.x, parent, parentTerm, cluster_id, everything())
  
  # --------- MDS（cmdscale）散点图准备并绘图 ----------
  # 使用 (1 - simMatrix) 作为距离（与 rrvgo 文档一致）
  dist_mat <- 1 - simMatrix
  # 保证对称、非负，转为 dist
  dist_mat[dist_mat < 0] <- 0
  d <- as.dist(dist_mat)
  mds_coords <- cmdscale(d, k = 2, add = FALSE)
  colnames(mds_coords) <- c("MDS1","MDS2")
  mds_df <- data.frame(ID = rownames(mds_coords),
                       MDS1 = mds_coords[,1],
                       MDS2 = mds_coords[,2],
                       stringsAsFactors = FALSE)
  # 将 cluster 信息并入
  mds_df <- mds_df %>% left_join(mapping_df %>% dplyr::select(ID, cluster_id, parent, parentTerm, score), by = "ID")
  # 颜色
  ncol <- length(unique(mds_df$cluster_id))
  palette <- colorRampPalette(brewer.pal(min(8, max(3, ncol)), "Set2"))(ncol)
  names(palette) <- sort(unique(mds_df$cluster_id))
  
  # 绘图：点按 cluster 上色，点大小按 score（若缺则按 Count 或 1）
  mds_df$size_for_plot <- ifelse(!is.na(mds_df$score), scales::rescale(mds_df$score, to = c(2,6), na.rm = TRUE),
                                 ifelse(!is.na(mds_df$Count), scales::rescale(mds_df$Count, to = c(2,6), na.rm = TRUE), 2))
  
  p_mds <- ggplot(mds_df, aes(x = MDS1, y = MDS2, color = factor(cluster_id), size = size_for_plot)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = palette, name = "cluster_id") +
    guides(size = FALSE) +
    theme_minimal() +
    ggtitle(paste0("MDS (cmdscale) plot - ", ont, " (clusters=", length(unique(mds_df$cluster_id)), ")")) +
    theme(plot.title = element_text(hjust = 0.5))
  # 在图上标注 representative parent（每个 cluster 取 parentTerm 的坐标即 parent 的点位置）
  reps <- mds_df %>% filter(ID == parent) %>% distinct(parent, parentTerm, cluster_id, MDS1, MDS2)
  # 若 representative 的点未包含（有时 parent 可能不在列表），尝试从 mapping_df 找对应 ID 再合并
  if(nrow(reps) < length(unique(mds_df$cluster_id))) {
    # 尝试提取 parent coords by matching
    rep_coords <- mds_df %>% group_by(cluster_id) %>%
      filter(ID == parent | row_number()==1) %>% slice(1) %>% ungroup()
    reps <- rep_coords %>% distinct(cluster_id, parent, parentTerm, MDS1, MDS2)
  }
  p_mds <- p_mds + geom_text_repel(data = reps, aes(x = MDS1, y = MDS2, label = parentTerm),
                                   size = 3, max.overlaps = 30, show.legend = FALSE)
  
  # 保存图片
  png_filename <- paste0("rrvgo_mds_", ont, ".png")
  ggsave(png_filename, p_mds, width = 8, height = 6)
  message("MDS 图保存为: ", png_filename)
  
  # --------- 词云（代表性 parent） ----------
  # 代表性 terms = 所有 parent 的 parentTerm（并以 parent 的 score 或 cluster 内最大 score 作为权重）
  parents_info <- mapping_df %>% group_by(parent, parentTerm) %>%
    summarize(weight = max(score, na.rm = TRUE), n_terms = n(), .groups = "drop")
  # 若 weight 全为 -Inf（说明 score 全为 NA），改用 n_terms 作为权重
  if(all(is.infinite(parents_info$weight))) {
    parents_info$weight <- parents_info$n_terms
  }
  # 若 parentTerm 有 NA，就用 getGoTerm 来取得
  missing_parentTerm <- is.na(parents_info$parentTerm) | parents_info$parentTerm == ""
  if(any(missing_parentTerm)) {
    parents_info$parentTerm[missing_parentTerm] <- getGoTerm(parents_info$parent[missing_parentTerm])
  }
  # 画词云
  wc_file <- paste0("rrvgo_wordcloud_", ont, ".png")
  png(wc_file, width = 1200, height = 800, res = 150)
  par(mar = c(0,0,0,0))
  set.seed(123)
  wordcloud(words = parents_info$parentTerm,
            freq = parents_info$weight,
            min.freq = 1,
            random.order = FALSE,
            scale = c(4, 0.6),
            rot.per = 0.1)
  dev.off()
  message("词云保存为: ", wc_file)
  
  # 返回结果
  list(
    mapping_df = mapping_df,
    simMatrix = simMatrix,
    reduced = reduced,
    mds_plot = p_mds,
    wordcloud_file = wc_file,
    mds_file = png_filename
  )
}

# -------------------- 对三个数据框分别运行 --------------------
# 假设你的数据框变量名就是 ego_bp_df, ego_mf_df, ego_cc_df
res_bp <- process_rrvgo_module(ego_bp_df[1:500,], ont = "BP", orgdb = "org.Hs.eg.db", max_clusters = 20)
res_mf <- process_rrvgo_module(ego_mf_df, ont = "MF", orgdb = "org.Hs.eg.db", max_clusters = 20)
res_cc <- process_rrvgo_module(ego_cc_df, ont = "CC", orgdb = "org.Hs.eg.db", max_clusters = 20)

# 最终输出的数据框（每个条目属于哪个聚类）
bp_clusters_df <- res_bp$mapping_df
mf_clusters_df <- res_mf$mapping_df
cc_clusters_df <- res_cc$mapping_df

# 保存到 CSV（可选）
write.csv(bp_clusters_df, "bp_clusters_mapping.csv", row.names = FALSE)
write.csv(mf_clusters_df, "mf_clusters_mapping.csv", row.names = FALSE)
write.csv(cc_clusters_df, "cc_clusters_mapping.csv", row.names = FALSE)
message("三份 mapping csv 已保存：bp_clusters_mapping.csv, mf_clusters_mapping.csv, cc_clusters_mapping.csv")

# 返回这些对象供后续分析 / 交互式查看
list(bp = res_bp, mf = res_mf, cc = res_cc)


#画气泡图#############################################################################################
# plot_enrichment_bubbles.R
# 依赖包：readxl, dplyr, stringr, ggplot2, viridis (可选)
# install.packages(c("readxl","dplyr","stringr","ggplot2","viridis"))

library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)   # color scale
library(scales)    # for pretty numbers

# --- 配置 -------------------------------------------------------------
base_path <- "C:/毕业论文/数据处理及绘图/富集分析"
files <- c("mf_clusters_mapping.xlsx", "bp_clusters_mapping.xlsx", "cc_clusters_mapping.xlsx")
out_dir <- file.path(base_path, "bubble_plots_output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 辅助函数：把 "87/1250" 之类的 GeneRatio 转为数值 87/1250
parse_ratio <- function(ratio_str) {
  # 处理 NA 或空字符
  if (is.na(ratio_str) || trimws(ratio_str) == "") return(NA_real_)
  # 允许形如 "87/1250" 或 " 87 / 1250"
  parts <- str_split(ratio_str, "/", simplify = TRUE)
  if (ncol(parts) < 2) return(NA_real_)
  num <- as.numeric(str_trim(parts[1]))
  den <- as.numeric(str_trim(parts[2]))
  if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
  return(num / den)
}

# 主处理与绘图函数
process_and_plot <- function(xlsx_path, out_dir) {
  message("Processing: ", xlsx_path)
  df <- read_excel(xlsx_path)
  
  # 确保列存在并转成合适类型（不抛错）
  required_cols <- c("cluster_id","score","GeneRatio","p.adjust","Count","Description")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("缺少必要列: ", paste(missing_cols, collapse = ", "))
  }
  
  df <- df %>%
    mutate(
      score = as.numeric(score),
      p.adjust = as.numeric(p.adjust),
      Count = as.numeric(Count),
      GeneRatio_raw = as.character(GeneRatio),
      GeneRatio_num = vapply(GeneRatio_raw, parse_ratio, numeric(1))
    )
  
  # 按 cluster_id 分组，取 score 最大的一行（with_ties = FALSE -> 若并列只取第一个）
  top_per_cluster <- df %>%
    group_by(cluster_id) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # 为绘图准备字段（处理可能的 NA / 0）
  top_per_cluster <- top_per_cluster %>%
    mutate(
      # 防止 p.adjust 为 0 导致 -log10(0) Inf
      p_adj_safe = ifelse(is.na(p.adjust), NA_real_, p.adjust + 1e-300),
      neglog10_padj = -log10(p_adj_safe),
      # 若 Count 非数值则设置为 1
      Count_plot = ifelse(is.na(Count) | Count <= 0, 1, Count)
    )
  
  # 按 GeneRatio_num 对 Description 排序（y 轴显示 Description）
  # 也可以改为按 score 或 cluster_id 排序，按需修改
  top_per_cluster <- top_per_cluster %>%
    mutate(Description_plot = ifelse(is.na(Description), paste0("ID:", ID), Description))
  
  # 绘图
  p <- ggplot(top_per_cluster, 
              aes(x = GeneRatio_num, y = reorder(Description_plot, GeneRatio_num),
                  size = Count_plot, color = neglog10_padj)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(range = c(3, 10), name = "Gene count (Count)") +
    scale_color_viridis(option = "D", direction = -1, na.value = "grey80",
                        name = expression(-log[10](p.adjust))) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
      title = paste0("Enrichment bubble plot — ", basename(xlsx_path)),
      x = "GeneRatio (fraction)",
      y = NULL,
      caption = "each cluster_id chose highest score ；bubble_size = Count，color = -log10(p.adjust)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(size = 11),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  print(p)
  # 保存图像和 csv（输出选项：PNG & PDF）
  base_name <- tools::file_path_sans_ext(basename(xlsx_path))
  png_file <- file.path(out_dir, paste0(base_name, "_bubble.png"))
  pdf_file <- file.path(out_dir, paste0(base_name, "_bubble.pdf"))
  csv_file <- file.path(out_dir, paste0(base_name, "_top_per_cluster.csv"))
  
  ggsave(png_file, p, width = 10, height = max(6, 0.25 * nrow(top_per_cluster)), dpi = 300)
  ggsave(pdf_file, p, width = 10, height = max(6, 0.25 * nrow(top_per_cluster)))
  write.csv(top_per_cluster, csv_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  message("Saved: ", png_file)
  message("Saved: ", pdf_file)
  message("Saved top-per-cluster csv: ", csv_file)
  
  return(list(plot = p, top_df = top_per_cluster))
}

# 批量处理三个文件
results <- list()
for (f in files) {
  fp <- file.path(base_path, f)
  if (!file.exists(fp)) {
    warning("文件不存在，跳过: ", fp)
    next
  }
  results[[f]] <- process_and_plot(fp, out_dir)
}

message("全部完成。输出文件夹：", out_dir)



















