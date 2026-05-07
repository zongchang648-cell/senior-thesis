
library(readxl)
library(dplyr)
library(tibble)
library(igraph)
library(devtools)
library(ggplot2)
# 安装 ne.PCA（如果尚未安装）
# README 中建议使用 devtools::install_github("wf-frank2019/ne.PCA")
#if (!"ne.PCA" %in% installed.packages()[,1]) {
#  message("Installing ne.PCA from GitHub...")
#  devtools::install_github("wf-frank2019/ne.PCA")
#}
library(ne.PCA)


#file<- "C:/毕业论文/数据/STRING_output_confidence0.7_about400nodes.xlsx"
file<- "C:/毕业论文/数据/PPIN/PPIN_with_genes.xlsx"
nodes<-read_xlsx(file,sheet = 3)
edges_raw<-read_xlsx(file,sheet = 2)
edges_raw<-edges_raw[,1:2]

# 假设 edges_raw 的前两列是端点
a <- as.character(edges_raw[[1]])
b <- as.character(edges_raw[[2]])

u <- pmin(a, b)   # 较小的放第一列（字典序）
v <- pmax(a, b)   # 较大的放第二列

edges_unique <- unique(data.frame(from = u, to = v, stringsAsFactors = FALSE))

# 可选：去掉自环（如果不想保留 A-A）
edges<- subset(edges_unique, from != to)
colnames(edges)<-c('gene symbol A',"gene symbol B")
nodes<-as.data.frame(nodes)
edges<-as.data.frame(edges)

#基于原始PPI获得RWR分数排序以及核心子网
RWRs.CN(nodes, edges, r = 0.8, core = 10)
#基于新的PPI获得每条边的权重
edges<-read_xlsx('filtered_ppi_edges.xlsx',sheet = 1)
edges<-as.data.frame(edges)
TFCs(edges)
#对核心子网进行模块划分，与文章方案一致
edges_core<-read_xlsx('top10%_CoreNet.xlsx')
edges_core<-as.data.frame(edges_core)
#edges_core<-as.data.frame(edges_core)
edges<-as.data.frame(read_xlsx("TFC_8525interactions.xlsx",sheet = 1))
edges<-edges[,1:2]
MODULEs(edges_core, heatmap = TRUE, tarGet = TRUE)

###########################################################

nodes_core<-read_xlsx('top10%_CoreNet.xlsx',sheet=2)
nodes_core<-as.data.frame(nodes_core)

library(readxl)
library(writexl)
library(dplyr)

# ---- 参数：文件名（或改成完整路径） ----
input_xlsx  <- "TFC_3822interactions.xlsx"          # 输入文件（可改为完整路径）
output_xlsx <- "TFC_3822interactions_filtered.xlsx" # 输出文件名

# ---- 检查 nodes_core 是否存在并提取节点向量 ----
if(!exists("nodes_core")) stop("未找到 R 对象 'nodes_core'。请先在环境中准备好该数据框。")

# 假设 nodes_core 是只有一列的 data.frame；提取为字符向量并去重
nodes_vec <- unique(as.character(nodes_core[[1]]))
nodes_vec <- nodes_vec[!is.na(nodes_vec)]

if(length(nodes_vec) == 0) stop("'nodes_core' 中没有可用节点名。")

# ---- 读取 TFC xlsx ----
tfc_all <- readxl::read_excel(input_xlsx)

# 确保列名为 name1, name2（如果列名不同，请相应修改下面的列名）
if(!all(c("name1","name2") %in% colnames(tfc_all))) {
  stop("输入文件中没有找到 'name1' 和 'name2' 两列，请检查列名。")
}

# ---- 过滤：只保留 name1 和 name2 都在 nodes_core 的行 ----
tfc_filtered <- tfc_all %>%
  filter((name1 %in% nodes_vec) & (name2 %in% nodes_vec))

# ---- 保存为新的 xlsx ----
writexl::write_xlsx(tfc_filtered, path = output_xlsx)

# ---- 信息反馈 ----
message("过滤完成：")
message("原始行数： ", nrow(tfc_all))
message("保留行数：   ", nrow(tfc_filtered))
message("已保存为：   ", normalizePath(output_xlsx))

##############富集分析####################
nodes_module<-read_xlsx('node_Module.xlsx')
nodes_module<-as.data.frame(nodes_module)



# ======================
# 模块富集分析完整代码
# ======================

# 加载必要的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(readr)

# === 基本参数 ===
modules_to_analyze <- c(1, 2, 3, 6, 7, 8, 9)
out_prefix <- "MODULE_enrichment"  # 输出文件前缀

# === 检查 nodes_module ===
if(!exists("nodes_module")) {
  stop("请先在 R 环境中准备好 'nodes_module' 数据框，且包含 'name' 和 'module' 两列。")
}

if(!all(c("name", "module") %in% colnames(nodes_module))) {
  stop("nodes_module 必须包含 name 和 module 两列。")
}

# 去除没有模块的基因（NA）
nodes_df <- nodes_module %>%
  dplyr::filter(!is.na(module)) %>%
  mutate(module = as.numeric(module))  # 确保 module 是数值

# 过滤只保留要分析的模块（若模块不存在会跳过）
nodes_df <- nodes_df %>% filter(module %in% modules_to_analyze)

if(nrow(nodes_df) == 0) {
  stop("在 nodes_module 中未找到待分析模块内的基因。请检查 modules_to_analyze 与你的数据。")
}

# === 预映射 SYMBOL -> ENTREZID（全体出现过的基因） ===
all_genes <- unique(as.character(nodes_df$name))
# 进行基因 ID 转换（可能有部分无法映射）
map_df <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
# map_df: SYMBOL, ENTREZID
# 若有未映射基因，请在后续注释中注意

# 用于汇总所有模块显著条目的列表（后面画合并热图）
GO_terms_all <- list()
KEGG_terms_all <- list()

# === 针对每个模块做富集并输出（csv + dotplot png） ===
for(mod in modules_to_analyze) {
  genes_mod <- nodes_df %>% filter(module == mod) %>% pull(name) %>% unique()
  
  if(length(genes_mod) < 5) {
    message(sprintf("模块 %s 基因数 < 5，跳过富集（建议至少 5 个以获得稳健结果）。", mod))
    next
  }
  
  # SYMBOL -> ENTREZ
  mapped <- map_df %>% filter(SYMBOL %in% genes_mod)
  entrez_ids <- unique(mapped$ENTREZID)
  
  if(length(entrez_ids) == 0) {
    message(sprintf("模块 %s 无法映射到 ENTREZ ID，跳过。", mod))
    next
  }
  
  # --- GO Biological Process 富集（不直接用 pvalue cutoff，后面用 p.adjust 筛选） ---
  ego <- tryCatch({
    enrichGO(gene = entrez_ids,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "BP",
             pvalueCutoff = 1,
             qvalueCutoff = 1,
             readable = TRUE)
  }, error = function(e) NULL)
  
  if(!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego_df <- as.data.frame(ego)
    ego_sig <- ego_df %>% filter(p.adjust <= 0.05)
    
    # 保存结果
    write_csv(ego_df, paste0(out_prefix, "_module", mod, "_GO_BP_all.csv"))
    write_csv(ego_sig, paste0(out_prefix, "_module", mod, "_GO_BP_sig.csv"))
    
    # dotplot（气泡图）
    if(nrow(ego_sig) > 0) {
      png(filename = paste0(out_prefix, "_module", mod, "_GO_BP_dotplot.png"),
          width = 1500, height = 1800, res = 150)
      print(dotplot(ego, showCategory = 20) + 
              ggtitle(paste0("Module ", mod, " GO:BP (top sig categories)")))
      dev.off()
      
      # 保存显著条目名字，用于合并热图
      GO_terms_all[[as.character(mod)]] <- ego_sig$Description
    } else {
      message(sprintf("模块 %s 的 GO:BP 没有 p.adjust <= 0.05 的条目。", mod))
    }
  } else {
    message(sprintf("模块 %s 的 GO:BP 富集失败或没有结果。", mod))
  }
  
  # --- KEGG 富集（注意需要 ENTREZ ID 且 organism='hsa'） ---
  ekegg <- tryCatch({
    enrichKEGG(gene = entrez_ids,
               organism = "hsa",
               pvalueCutoff = 1,
               qvalueCutoff = 1)
  }, error = function(e) NULL)
  
  if(!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg_df <- as.data.frame(ekegg)
    # 将 ENTREZ id 转为 readable gene symbols（clusterProfiler 提供 setReadable）
    ekegg_readable <- tryCatch(
      setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"), 
      error = function(e) ekegg
    )
    ekegg_df_r <- as.data.frame(ekegg_readable)
    ekegg_sig <- ekegg_df_r %>% filter(p.adjust <= 0.05)
    
    write_csv(ekegg_df_r, paste0(out_prefix, "_module", mod, "_KEGG_all.csv"))
    write_csv(ekegg_sig, paste0(out_prefix, "_module", mod, "_KEGG_sig.csv"))
    
    if(nrow(ekegg_sig) > 0) {
      png(filename = paste0(out_prefix, "_module", mod, "_KEGG_dotplot.png"),
          width = 1500, height = 1800, res = 150)
      print(dotplot(ekegg_readable, showCategory = 20) + 
              ggtitle(paste0("Module ", mod, " KEGG (top sig pathways)")))
      dev.off()
      
      KEGG_terms_all[[as.character(mod)]] <- ekegg_sig$Description
    } else {
      message(sprintf("模块 %s 的 KEGG 没有 p.adjust <= 0.05 的条目。", mod))
    }
  } else {
    message(sprintf("模块 %s 的 KEGG 富集失败或没有结果。", mod))
  }
} # for modules

# === 合并热图：把每个富集类型在所有模块的显著条目合并，构建 -log10(p.adjust) 矩阵并画热图 ===

# Helper：读取单模块结果 CSV 并抽取 p.adjust （若没有则跳过）
collect_adjust <- function(pattern, modules, prefix) {
  # pattern 应该是 like "_module{mod}_GO_BP_all.csv" etc.
  res_list <- list()
  for(mod in modules) {
    fn_all <- paste0(prefix, "_module", mod, pattern)
    if(!file.exists(fn_all)) next
    
    # 关键修正：使用 as.data.frame() 转换 spec_tbl_df 为普通 data.frame
    df <- as.data.frame(read_csv(fn_all, show_col_types = FALSE))
    
    if(!all(c("Description", "p.adjust") %in% colnames(df))) next
    
    res_list[[as.character(mod)]] <- df %>% dplyr::select(Description, p.adjust)
  }
  return(res_list)
}

# GO heatmap（使用 GO_all files）
go_res <- collect_adjust("_GO_BP_all.csv", modules_to_analyze, out_prefix)

if(length(go_res) > 0) {
  # 合并所有显著条目（p.adjust <= 0.05）作为行集合
  all_terms <- unique(unlist(lapply(go_res, function(df) {
    df$Description[df$p.adjust <= 0.05]
  })))
  
  if(length(all_terms) > 0) {
    mat <- matrix(0, 
                  nrow = length(all_terms), 
                  ncol = length(modules_to_analyze),
                  dimnames = list(all_terms, paste0("mod", modules_to_analyze)))
    
    for(i in seq_along(modules_to_analyze)) {
      mod <- modules_to_analyze[i]
      key <- as.character(mod)
      
      if(!is.null(go_res[[key]])) {
        df <- go_res[[key]]
        
        # for terms in all_terms, fill -log10(p.adjust) if p.adjust<=0.05 else 0
        for(term in all_terms) {
          row <- df %>% filter(Description == term)
          
          if(nrow(row) > 0 && row$p.adjust[1] <= 0.05) {
            mat[term, i] <- -log10(row$p.adjust[1])
          }
        }
      }
    }
    
    # 限制行数以保证热图可视（若行太多，先取 -log10 总和排序前 80）
    if(nrow(mat) > 120) {
      row_order_score <- rowSums(mat)
      top_rows <- names(sort(row_order_score, decreasing = TRUE))[1:120]
      mat <- mat[top_rows, , drop = FALSE]
    }
    
    png(filename = paste0(out_prefix, "_GO_BP_heatmap_allModules.png"), 
        width = 1600, height = 1800, res = 150)
    
    pheatmap::pheatmap(
      mat, 
      cluster_rows = TRUE, 
      cluster_cols = TRUE, 
      main = "GO:BP -log10(p.adjust) across modules", 
      color = colorRampPalette(c("white", "yellow", "red"))(50),
      fontsize_row = 8,
      fontsize_col = 12
    )
    dev.off()
  } else {
    message("所有模块的 GO:BP 中均无 p.adjust <= 0.05 的条目，未生成 GO 热图。")
  }
} else {
  message("未找到任何模块 GO 结果文件，跳过 GO 合并热图。")
}

# KEGG heatmap
kegg_res <- collect_adjust("_KEGG_all.csv", modules_to_analyze, out_prefix)

if(length(kegg_res) > 0) {
  all_terms_kegg <- unique(unlist(lapply(kegg_res, function(df) {
    df$Description[df$p.adjust <= 0.05]
  })))
  
  if(length(all_terms_kegg) > 0) {
    matk <- matrix(0, 
                   nrow = length(all_terms_kegg), 
                   ncol = length(modules_to_analyze),
                   dimnames = list(all_terms_kegg, paste0("mod", modules_to_analyze)))
    
    for(i in seq_along(modules_to_analyze)) {
      mod <- modules_to_analyze[i]
      key <- as.character(mod)
      
      if(!is.null(kegg_res[[key]])) {
        df <- kegg_res[[key]]
        
        for(term in all_terms_kegg) {
          row <- df %>% filter(Description == term)
          
          if(nrow(row) > 0 && row$p.adjust[1] <= 0.05) {
            matk[term, i] <- -log10(row$p.adjust[1])
          }
        }
      }
    }
    
    if(nrow(matk) > 120) {
      row_order_score <- rowSums(matk)
      top_rows <- names(sort(row_order_score, decreasing = TRUE))[1:120]
      matk <- matk[top_rows, , drop = FALSE]
    }
    
    png(filename = paste0(out_prefix, "_KEGG_heatmap_allModules.png"), 
        width = 1600, height = 1800, res = 150)
    
    pheatmap::pheatmap(
      matk, 
      cluster_rows = TRUE, 
      cluster_cols = TRUE,
      main = "KEGG -log10(p.adjust) across modules",
      color = colorRampPalette(c("white", "lightblue", "blue"))(50),
      fontsize_row = 8,
      fontsize_col = 12
    )
    dev.off()
  } else {
    message("所有模块的 KEGG 中均无 p.adjust <= 0.05 的条目，未生成 KEGG 热图。")
  }
} else {
  message("未找到任何模块 KEGG 结果文件，跳过 KEGG 合并热图。")
}

# === 生成富集结果汇总报告 ===
message("\n====== 富集分析完成 ======")
message("输出文件如下：")

# 列出生成的GO文件
go_files <- list.files(pattern = paste0(out_prefix, ".*GO.*\\.(csv|png)$"))
if(length(go_files) > 0) {
  message("\nGO富集结果文件：")
  cat(paste(go_files, collapse = "\n"))
}

# 列出生成的KEGG文件
kegg_files <- list.files(pattern = paste0(out_prefix, ".*KEGG.*\\.(csv|png)$"))
if(length(kegg_files) > 0) {
  message("\n\nKEGG富集结果文件：")
  cat(paste(kegg_files, collapse = "\n"))
}

# 列出热图文件
heatmap_files <- list.files(pattern = paste0(out_prefix, ".*heatmap.*\\.png$"))
if(length(heatmap_files) > 0) {
  message("\n\n合并热图文件：")
  cat(paste(heatmap_files, collapse = "\n"))
}

message("\n\n分析总结：")
message(sprintf("1. 分析了 %d 个模块: %s", length(modules_to_analyze), paste(modules_to_analyze, collapse = ", ")))
message(sprintf("2. 输入基因总数: %d", nrow(nodes_df)))
message(sprintf("3. 成功映射到ENTREZID的基因数: %d", length(unique(map_df$SYMBOL))))
message(sprintf("4. 每个模块平均基因数: %.1f", mean(table(nodes_df$module))))

message("\n全部完成。各模块的富集结果 CSV 和图像已保存到当前工作目录。")