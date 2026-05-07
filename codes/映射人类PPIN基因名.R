

# 加载包
library(readxl)
library(biomaRt)
library(dplyr)
library(stringr)

# 1. 读取数据
file_path <- "C:/毕业论文/数据/PPIN/PPIN.xlsx"
ppi_data <- read_excel(file_path)

# 查看数据结构
str(ppi_data)
head(ppi_data)

# 2. 提取所有唯一的ENSP ID（去掉前缀"9606."）
# 从第一列提取
protein1_ids <- unique(str_replace(ppi_data$protein1, "^9606\\.", ""))
# 从第二列提取
protein2_ids <- unique(str_replace(ppi_data$protein2, "^9606\\.", ""))

# 合并所有唯一ID
all_ensp_ids <- unique(c(protein1_ids, protein2_ids))
cat("需要转换的蛋白ID数量:", length(all_ensp_ids), "\n")

# 3. 使用biomaRt进行ID转换
# 连接到Ensembl数据库
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# 尝试用多种方法获取基因名
# 方法1：直接从ENSP ID转换
cat("正在从Ensembl数据库获取基因名映射...\n")
tryCatch({
  # 尝试获取ENSP到Gene name的映射
  gene_mapping <- getBM(
    attributes = c("ensembl_peptide_id", "hgnc_symbol", "ensembl_gene_id"),
    filters = "ensembl_peptide_id",
    values = all_ensp_ids,
    mart = ensembl
  )
  
  # 重命名列以便理解
  colnames(gene_mapping) <- c("ensp_id", "gene_name", "ensembl_gene_id")
  
  # 显示部分映射结果
  cat("成功获取映射:", nrow(gene_mapping), "条记录\n")
  head(gene_mapping)
  
}, error = function(e) {
  cat("方法1失败:", e$message, "\n")
  
  # 方法2：如果方法1失败，尝试通过ENSG ID转换
  cat("尝试通过ENSG ID转换...\n")
  
  # 首先获取ENSP到ENSG的映射
  ensp_to_ensg <- getBM(
    attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
    filters = "ensembl_peptide_id",
    values = all_ensp_ids,
    mart = ensembl
  )
  
  # 然后获取ENSG到Gene name的映射
  if(nrow(ensp_to_ensg) > 0) {
    ensg_ids <- unique(ensp_to_ensg$ensembl_gene_id)
    ensg_to_gene <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "ensembl_gene_id",
      values = ensg_ids,
      mart = ensembl
    )
    
    # 合并两个映射
    gene_mapping <- merge(ensp_to_ensg, ensg_to_gene, 
                          by = "ensembl_gene_id", all.x = TRUE)
    colnames(gene_mapping) <- c("ensembl_gene_id", "ensp_id", "gene_name")
  }
})

# 4. 创建映射字典
# 移除空的基因名
gene_mapping_clean <- gene_mapping %>%
  filter(!is.na(gene_name) & gene_name != "")

# 创建从ENSP ID到基因名的映射
ensp_to_gene <- setNames(gene_mapping_clean$gene_name, gene_mapping_clean$ensp_id)

# 统计映射成功率
mapped_ids <- sum(all_ensp_ids %in% gene_mapping_clean$ensp_id)
cat(sprintf("映射成功率: %d/%d (%.1f%%)\n", 
            mapped_ids, length(all_ensp_ids), 
            mapped_ids/length(all_ensp_ids)*100))

# 5. 应用映射到原始数据
# 创建辅助函数来转换ID
convert_id <- function(full_id) {
  # 去掉"9606."前缀
  ensp_id <- str_replace(full_id, "^9606\\.", "")
  
  # 查找对应的基因名
  gene_name <- ensp_to_gene[ensp_id]
  
  # 如果找不到，返回原始ENSP ID
  if(is.na(gene_name) || is.null(gene_name)) {
    return(ensp_id)  # 或者返回NA，这里返回ENSP ID以便追踪
  } else {
    return(gene_name)
  }
}

# 应用转换
ppi_data_converted <- ppi_data %>%
  mutate(
    gene1 = sapply(protein1, convert_id),
    gene2 = sapply(protein2, convert_id)
  ) %>%
  # 重新排列列，将基因名放在前面
  dplyr::select(gene1, gene2, combined_score, protein1, protein2)

# 6. 查看转换结果
cat("\n转换后的数据前几行:\n")
head(ppi_data_converted)

# 统计未成功转换的ID
unmapped1 <- sum(!sapply(ppi_data$protein1, function(x) {
  ensp_id <- str_replace(x, "^9606\\.", "")
  ensp_id %in% gene_mapping_clean$ensp_id
}))

unmapped2 <- sum(!sapply(ppi_data$protein2, function(x) {
  ensp_id <- str_replace(x, "^9606\\.", "")
  ensp_id %in% gene_mapping_clean$ensp_id
}))

cat(sprintf("\n未成功转换的protein1数量: %d/%d\n", unmapped1, nrow(ppi_data)))
cat(sprintf("未成功转换的protein2数量: %d/%d\n", unmapped2, nrow(ppi_data)))

# 7. 保存结果到新文件
output_path <- "C:/毕业论文/数据/PPIN/PPIN_with_genes.xlsx"

# 如果需要保存为Excel
if (!require("writexl")) install.packages("writexl")
library(writexl)
write_xlsx(ppi_data_converted, output_path)

# 或者保存为CSV（更通用）
write.csv(ppi_data_converted, 
          file = "C:/毕业论文/数据/PPIN/PPIN_with_genes.csv",
          row.names = FALSE)

cat("\n转换完成！结果已保存到:\n")
cat("Excel文件:", output_path, "\n")
cat("CSV文件: C:/毕业论文/数据/PPIN/PPIN_with_genes.csv\n")

# 8. 可选：查看一些示例映射
cat("\n=== 示例映射 ===\n")
sample_data <- head(ppi_data_converted, 10)
for(i in 1:nrow(sample_data)) {
  cat(sprintf("%s (%s) -- %s (%s): %d\n",
              sample_data$gene1[i], 
              str_replace(sample_data$protein1[i], "^9606\\.", ""),
              sample_data$gene2[i],
              str_replace(sample_data$protein2[i], "^9606\\.", ""),
              sample_data$combined_score[i]))
}