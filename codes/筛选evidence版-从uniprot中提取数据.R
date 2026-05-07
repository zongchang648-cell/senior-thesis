library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(writexl)

# --------- 配置 URL ----------
url <- "https://www.ebi.ac.uk/proteins/api/variation/P15056?format=json"

# 直接从 URL 读取 JSON 数据
variation_list <- fromJSON(url)

feature <- variation_list$features

# ---------------------------------------------------------------------------
# 提取 clinicalSignificances 第一项的通用函数（你原来的）
extract_first_cell_base <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.atomic(x) && length(x) == 0) return(NA_character_)
  if (is.data.frame(x) && nrow(x) == 0) return(NA_character_)
  
  if (is.data.frame(x)) {
    val <- x[[1]][1]
  } else if (is.list(x)) {
    first_elem <- x[[1]]
    if (is.data.frame(first_elem)) val <- first_elem[[1]][1] else val <- first_elem[1]
  } else {
    val <- x[1]
  }
  
  if (is.factor(val)) val <- as.character(val)
  if (is.null(val) || length(val) == 0 || is.na(val)) return(NA_character_)
  return(trimws(as.character(val)))
}

# 对每一行应用
clinical_type_vec <- vapply(feature$clinicalSignificances, extract_first_cell_base, FUN.VALUE = character(1), USE.NAMES = FALSE)

# 把空字符串设为 NA
clinical_type_vec[clinical_type_vec == ""] <- NA_character_

# 新增列
feature$clinical_type <- clinical_type_vec

# 过滤：不区分大小写精确匹配 "likely pathogenic" 或 "pathogenic"
keep_mask <- grepl("^(likely pathogenic|pathogenic)$", feature$clinical_type, ignore.case = TRUE)
feature_filtered <- feature[which(!is.na(feature$clinical_type) & keep_mask), , drop = FALSE]

# ----------------------------------------------------------------------------------------------
# 提取 association->disease（与你原逻辑类似）
extract_firstcol_concat_base <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.matrix(x) || is.data.frame(x)) {
    if (ncol(x) < 1 || nrow(x) < 1) return(NA_character_)
    col1 <- as.character(x[[1]])
    col1 <- col1[!is.na(col1) & trimws(col1) != ""]
    if (length(col1) == 0) return(NA_character_)
    return(paste(col1, collapse = "/"))
  }
  if (is.list(x)) {
    first_elem <- x[[1]]
    if (is.data.frame(first_elem) || is.matrix(first_elem)) {
      return(extract_firstcol_concat_base(first_elem))
    } else {
      vec <- as.character(unlist(x))
      vec <- vec[!is.na(vec) & trimws(vec) != ""]
      if (length(vec) == 0) return(NA_character_)
      return(paste(vec, collapse = "/"))
    }
  }
  vec <- as.character(x)
  vec <- vec[!is.na(vec) & trimws(vec) != ""]
  if (length(vec) == 0) return(NA_character_)
  paste(vec, collapse = "/")
}

# 对每一行计算 disease 字符串
disease_vec <- vapply(feature_filtered$association, extract_firstcol_concat_base, FUN.VALUE = character(1), USE.NAMES = FALSE)

# 新增列
feature_filtered$disease <- disease_vec

# 过滤：移除 association 为 NULL 或 disease 为 NA/空 的行
keep <- !vapply(feature_filtered$association, is.null, FUN.VALUE = logical(1)) & !is.na(feature_filtered$disease) & feature_filtered$disease != ""
feature_filtered2 <- feature_filtered[keep, , drop = FALSE]

# 选择 consequenceType 为 missense
feature_2 <- feature_filtered2[feature_filtered2$consequenceType == "missense", , drop = FALSE]

# ------------------------------------------------------------------------
# 新增：检查第17列中 data.frame(s) 的 evidences 列是否存在 NULL/空 —— 若存在则过滤掉该行
# 说明：你提到 feature_2 的第17列是若干个 data.frame（或 list-of-data.frames），每个 data.frame 都有名为 evidences 的列（通常是第4列）
check_null_in_evidences <- function(cell) {
  # cell 可能为 NULL、data.frame、或 list-of-data.frames
  if (is.null(cell)) return(TRUE)  # 如果整格是 NULL，则认为有问题（要过滤）
  
  # 规范成 data.frame 列表
  dfs <- list()
  if (is.data.frame(cell)) {
    dfs <- list(cell)
  } else if (is.list(cell) && length(cell) > 0 && any(sapply(cell, function(x) is.data.frame(x) || is.null(x)))) {
    # 可能是 list-of-data.frames（或其中包含NULL）
    dfs <- cell
  } else {
    # 其它类型：尝试当成单个 data.frame（防御性处理）
    if (is.data.frame(cell)) {
      dfs <- list(cell)
    } else {
      # 未能识别结构，视为有问题（安全起见过滤掉）
      return(TRUE)
    }
  }
  
  # 对每个 data.frame 检查 evidences 列
  for (df in dfs) {
    if (is.null(df)) return(TRUE)
    if (!is.data.frame(df)) next
    # 如果没有 evidences 列，则跳过这个 df（不把它当作错误）
    if (!("evidences" %in% names(df))) next
    
    col <- df[["evidences"]]
    # 遍历 evidences 的每个单元格
    for (k in seq_along(col)) {
      v <- col[[k]]
      # 如果完全是 NULL
      if (is.null(v)) return(TRUE)
      # 如果是 list 且长度为0（空列表）也视为问题
      if (is.list(v) && length(v) == 0) return(TRUE)
      # 如果原子向量但长度为0，也视为问题
      if (!is.list(v) && length(v) == 0) return(TRUE)
      # 如果原子向量且全部为 NA，则视为问题
      if (is.atomic(v) && all(is.na(v))) return(TRUE)
      # 如果字符并且全部为空字符串或空白，也视为问题
      if (is.character(v) && all(trimws(v) == "")) return(TRUE)
      # 其他情况认为这一行的 evidences 有内容，继续检查下一行
    }
  }
  # 如果遍历完都没发现问题，返回 FALSE（表示没有 NULL/空）
  FALSE
}

# 确认 feature_2 至少有 17 列，然后对第17列逐行检测并过滤
if (is.data.frame(feature_2) && ncol(feature_2) >= 17) {
  col17 <- feature_2[[17]]  # 取第17列（list-column）
  # 对每个单元应用检查函数（注意：如果 col17 长度不等于 nrow(feature_2)，vapply 仍可用于列表）
  has_null_evidence <- vapply(col17, check_null_in_evidences, logical(1))
  # 过滤掉那些含有 NULL/空 evidence 的行
  feature_2 <- feature_2[!has_null_evidence, , drop = FALSE]
} else {
  stop("找不到第17列：请确认 feature_2 是 data.frame 且至少有 17 列。")
}

# ------------------------------------------------------------------------
# 你原来的 safe_extract 与 locations -> abstract 提取逻辑（保持不变）
safe_extract <- function(df) {
  if (is.null(df) || (is.data.frame(df) && nrow(df) == 0)) return("")
  first_col <- tryCatch(as.character(df[[1]]), error = function(e) character(0))
  if (length(first_col) == 0) return("")
  first_col <- first_col[!is.na(first_col)]
  first_col <- trimws(first_col)
  first_col <- first_col[first_col != ""]
  matches <- first_col[grepl("^p\\.", first_col)]
  if (length(matches) == 0) return("")
  uniq <- unique(matches)
  paste(uniq, collapse = "/")
}

if (is.data.frame(feature_2) && "locations" %in% names(feature_2)) {
  feature_2$abstract <- vapply(feature_2$locations, safe_extract, FUN.VALUE = character(1), USE.NAMES = FALSE)
} else if (is.list(feature_2) && !is.null(feature_2$locations)) {
  # 如果 feature_2 是 list 且有名为 locations 的元素，locations 也可能是一个 list-of-dfs
  feature_2$abstract <- vapply(feature_2$locations, safe_extract, FUN.VALUE = character(1), USE.NAMES = FALSE)
} else {
  stop("找不到 locations 列/元素：请确认 feature_2 是 data.frame/tibble（含 list-column 'locations'）或 list（含 'locations'）")
}

# ------------------------------------------------------------------------
# 选择列并写出结果（保留你原来的列索引）
# 注意：如果列索引越界请根据实际数据调整 c(2,3,9,10,11,22,23,24)
cols_to_keep <- c(2,3,9,10,11,22,23,24)
# 检查列索引是否存在
max_idx <- ncol(feature_2)
if (any(cols_to_keep > max_idx)) {
  warning("某些列索引超出 feature_2 的列数，请确认列索引。将仅选择存在的列。")
  cols_to_keep <- cols_to_keep[cols_to_keep <= max_idx]
}
feature_filtered3 <- feature_2[, cols_to_keep, drop = FALSE]

# 写出 xlsx
output_path <- 'C:\\毕业论文\\数据\\BRAF\\uniprot_processed.xlsx'
write_xlsx(feature_filtered3, output_path)

message("处理完成，结果已写入：", output_path)
