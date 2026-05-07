# 修改后的富集脚本：增加 GO_CC 并更改输出路径
# 加载需要的包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)  # 用于字符串处理
library(scales)   # percent_format 用到

nodes_module<-read_xlsx("C:\\毕业论文\\数据处理及绘图\\ppi_build\\modules_output\\第一次运行\\modules_fast_greedy.xlsx",sheet = 2)

# === 基本参数设置 ===
modules_to_analyze <- c(1, 2, 3)

# 输出目录（注意：使用正斜杠，R 能正确识别）
output_dir <- "C:/毕业论文/数据处理及绘图/富集分析"
out_prefix <- file.path(output_dir, "MODULE_enrichment")

# 创建输出目录（如果不存在）
if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

# === 检查数据 ===
if(!exists("nodes_module")) {
  stop("请先在R环境中准备好'nodes_module'数据框，且包含'name'和'module'两列。")
}
if(!all(c("name", "module") %in% colnames(nodes_module))) {
  stop("nodes_module必须包含name和module两列。")
}

# 处理数据
nodes_df <- nodes_module %>%
  dplyr::filter(!is.na(module)) %>%
  mutate(module = as.numeric(module)) %>%
  filter(module %in% modules_to_analyze)

if(nrow(nodes_df) == 0) {
  stop("在nodes_module中未找到模块1、2、3中的基因。")
}

# 所有基因的SYMBOL -> ENTREZID映射
all_genes <- unique(as.character(nodes_df$name))
map_df <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", 
               OrgDb = "org.Hs.eg.db", drop = TRUE)

# === 辅助函数：科学计数法格式化 ===
format_scientific <- function(x) {
  sapply(x, function(val) {
    if (is.na(val) || val == 0) {
      return("0")
    } else if (val < 0.001) {
      return(formatC(val, format = "e", digits = 1))
    } else {
      return(sprintf("%.3f", val))
    }
  })
}

# === 辅助函数：转换GeneRatio为数值 ===
convert_gene_ratio <- function(gene_ratio_str) {
  # 将 "10/100" 转换为数值比例 0.1
  parts <- str_split(gene_ratio_str, "/", simplify = TRUE)
  as.numeric(parts[,1]) / as.numeric(parts[,2])
}

# === 对每个模块进行富集分析并绘图 ===
for(mod in modules_to_analyze) {
  cat("正在分析模块", mod, "...\n")
  
  # 获取当前模块的基因
  genes_mod <- nodes_df %>% 
    filter(module == mod) %>% 
    pull(name) %>% 
    unique()
  
  if(length(genes_mod) < 5) {
    message(sprintf("模块 %s 基因数 < 5，跳过富集。", mod))
    next
  }
  
  # 转换为ENTREZID
  mapped <- map_df %>% filter(SYMBOL %in% genes_mod)
  entrez_ids <- unique(mapped$ENTREZID)
  
  if(length(entrez_ids) == 0) {
    message(sprintf("模块 %s 无 ENTREZ ID，跳过。", mod))
    next
  }
  
  # ========== 1. GO_BP富集分析 ==========
  cat("  进行GO_BP富集分析...\n")
  ego_bp <- tryCatch({
    enrichGO(gene = entrez_ids,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "BP",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.05,
             readable = TRUE)
  }, error = function(e) {
    message(sprintf("模块 %s GO_BP富集失败: %s", mod, e$message))
    NULL
  })
  
  if(!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    # 保存结果
    ego_bp_df <- as.data.frame(ego_bp)
    write_csv(ego_bp_df, file.path(output_dir, paste0(basename(out_prefix), "_module", mod, "_GO_BP.csv")))
    
    # 绘制气泡图
    if(nrow(ego_bp_df) > 0) {
      # 选取前20个最显著的条目
      ego_bp_top <- ego_bp_df %>% 
        arrange(p.adjust) %>% 
        head(min(20, nrow(ego_bp_df)))
      
      # 转换GeneRatio为数值
      ego_bp_top$GeneRatio_numeric <- sapply(ego_bp_top$GeneRatio, convert_gene_ratio)
      
      p_bp <- ggplot(ego_bp_top, aes(x = GeneRatio_numeric, y = reorder(Description, -log10(p.adjust)))) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        scale_color_gradient(low = "blue", high = "red", 
                             name = "-log10(p.adjust)",
                             labels = function(x) format_scientific(10^(-x))) +
        scale_size(range = c(3, 8), name = "Gene Count") +
        scale_x_continuous(labels = percent_format(scale = 1)) +
        labs(
          title = paste0("Module ", mod, " - GO Biological Process"),
          x = "Gene Ratio",
          y = "GO Term"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank()
        )
      
      ggsave(file.path(output_dir, paste0(basename(out_prefix), "_module", mod, "_GO_BP_dotplot.png")),
             plot = p_bp, width = 12, height = 10, dpi = 300)
      
      cat("  GO_BP气泡图已保存\n")
    }
  }
  
  # ========== 2. GO_MF富集分析 ==========
  cat("  进行GO_MF富集分析...\n")
  ego_mf <- tryCatch({
    enrichGO(gene = entrez_ids,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "MF",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.05,
             readable = TRUE)
  }, error = function(e) {
    message(sprintf("模块 %s GO_MF富集失败: %s", mod, e$message))
    NULL
  })
  
  if(!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
    # 保存结果
    ego_mf_df <- as.data.frame(ego_mf)
    write_csv(ego_mf_df, file.path(output_dir, paste0(basename(out_prefix), "_module", mod, "_GO_MF.csv")))
    
    # 绘制气泡图
    if(nrow(ego_mf_df) > 0) {
      # 选取前20个最显著的条目
      ego_mf_top <- ego_mf_df %>% 
        arrange(p.adjust) %>% 
        head(min(20, nrow(ego_mf_df)))
      
      # 转换GeneRatio为数值
      ego_mf_top$GeneRatio_numeric <- sapply(ego_mf_top$GeneRatio, convert_gene_ratio)
      
      p_mf <- ggplot(ego_mf_top, aes(x = GeneRatio_numeric, y = reorder(Description, -log10(p.adjust)))) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        scale_color_gradient(low = "blue", high = "red", 
                             name = "-log10(p.adjust)",
                             labels = function(x) format_scientific(10^(-x))) +
        scale_size(range = c(3, 8), name = "Gene Count") +
        scale_x_continuous(labels = percent_format(scale = 1)) +
        labs(
          title = paste0("Module ", mod, " - GO Molecular Function"),
          x = "Gene Ratio",
          y = "GO Term"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank()
        )
      
      ggsave(file.path(output_dir, paste0(basename(out_prefix), "_module", mod, "_GO_MF_dotplot.png")),
             plot = p_mf, width = 12, height = 10, dpi = 300)
      
      cat("  GO_MF气泡图已保存\n")
    }
  }
  
  # ========== 2b. GO_CC富集分析（新增） ==========
  cat("  进行GO_CC富集分析...\n")
  ego_cc <- tryCatch({
    enrichGO(gene = entrez_ids,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "CC",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.05,
             readable = TRUE)
  }, error = function(e) {
    message(sprintf("模块 %s GO_CC富集失败: %s", mod, e$message))
    NULL
  })
  
  if(!is.null(ego_cc) && nrow(as.data.frame(ego_cc)) > 0) {
    # 保存结果
    ego_cc_df <- as.data.frame(ego_cc)
    write_csv(ego_cc_df, file.path(output_dir, paste0(basename(out_prefix), "_module", mod, "_GO_CC.csv")))
    
    # 绘制气泡图
    if(nrow(ego_cc_df) > 0) {
      # 选取前20个最显著的条目
      ego_cc_top <- ego_cc_df %>% 
        arrange(p.adjust) %>% 
        head(min(20, nrow(ego_cc_df)))
      
      # 转换GeneRatio为数值
      ego_cc_top$GeneRatio_numeric <- sapply(ego_cc_top$GeneRatio, convert_gene_ratio)
      
      p_cc <- ggplot(ego_cc_top, aes(x = GeneRatio_numeric, y = reorder(Description, -log10(p.adjust)))) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        scale_color_gradient(low = "blue", high = "red", 
                             name = "-log10(p.adjust)",
                             labels = function(x) format_scientific(10^(-x))) +
        scale_size(range = c(3, 8), name = "Gene Count") +
        scale_x_continuous(labels = percent_format(scale = 1)) +
        labs(
          title = paste0("Module ", mod, " - GO Cellular Component"),
          x = "Gene Ratio",
          y = "GO Term"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank()
        )
      
      ggsave(file.path(output_dir, paste0(basename(out_prefix), "_module", mod, "_GO_CC_dotplot.png")),
             plot = p_cc, width = 12, height = 10, dpi = 300)
      
      cat("  GO_CC气泡图已保存\n")
    }
  }
  
  # ========== 3. KEGG富集分析 ==========
  cat("  进行KEGG富集分析...\n")
  ekegg <- tryCatch({
    enrichKEGG(gene = entrez_ids,
               organism = "hsa",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)
  }, error = function(e) {
    message(sprintf("模块 %s KEGG富集失败: %s", mod, e$message))
    NULL
  })
  
  if(!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    # 转换为可读格式
    ekegg_readable <- tryCatch(
      setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"),
      error = function(e) ekegg
    )
    
    # 保存结果
    ekegg_df <- as.data.frame(ekegg_readable)
    write_csv(ekegg_df, file.path(output_dir, paste0(basename(out_prefix), "_module", mod, "_KEGG.csv")))
    
    # 绘制气泡图
    if(nrow(ekegg_df) > 0) {
      # 选取前20个最显著的条目
      ekegg_top <- ekegg_df %>% 
        arrange(p.adjust) %>% 
        head(min(20, nrow(ekegg_df)))
      
      # 转换GeneRatio为数值
      ekegg_top$GeneRatio_numeric <- sapply(ekegg_top$GeneRatio, convert_gene_ratio)
      
      p_kegg <- ggplot(ekegg_top, aes(x = GeneRatio_numeric, y = reorder(Description, -log10(p.adjust)))) +
        geom_point(aes(size = Count, color = -log10(p.adjust))) +
        scale_color_gradient(low = "green", high = "purple", 
                             name = "-log10(p.adjust)",
                             labels = function(x) format_scientific(10^(-x))) +
        scale_size(range = c(3, 8), name = "Gene Count") +
        scale_x_continuous(labels = percent_format(scale = 1)) +
        labs(
          title = paste0("Module ", mod, " - KEGG Pathways"),
          x = "Gene Ratio",
          y = "KEGG Pathway"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank()
        )
      
      ggsave(file.path(output_dir, paste0(basename(out_prefix), "_module", mod, "_KEGG_dotplot.png")),
             plot = p_kegg, width = 12, height = 10, dpi = 300)
      
      cat("  KEGG气泡图已保存\n")
    }
  }
  
  cat("模块", mod, "分析完成！\n\n")
}

# === 汇总报告 ===
cat("====== 富集分析完成 ======\n")
cat("分析总结：\n")
cat(sprintf("1. 分析模块数: %d (模块 %s)\n", length(modules_to_analyze), paste(modules_to_analyze, collapse = ", ")))
cat(sprintf("2. 输入基因总数: %d\n", nrow(nodes_df)))
cat(sprintf("3. 成功映射到ENTREZID的SYMBOL数: %d\n", length(unique(map_df$SYMBOL))))
cat(sprintf("4. 每个模块平均基因数: %.1f\n", mean(table(nodes_df$module))))

# 列出生成的文件（在 output_dir 中）
cat("\n生成的文件（位于：", output_dir, "）:\n", sep = "")
csv_files <- list.files(path = output_dir, pattern = paste0(basename(out_prefix), ".*\\.csv$"), full.names = TRUE)
png_files <- list.files(path = output_dir, pattern = paste0(basename(out_prefix), ".*\\.png$"), full.names = TRUE)

if(length(csv_files) > 0) {
  cat("富集结果CSV文件：\n")
  cat(paste("  ", csv_files, collapse = "\n"), "\n")
}

if(length(png_files) > 0) {
  cat("气泡图PNG文件：\n")
  cat(paste("  ", png_files, collapse = "\n"), "\n")
}

cat("\n全部完成！\n")