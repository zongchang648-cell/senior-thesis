library(jsonlite)
library(dplyr)
library(purrr)
library(stringr)
library(writexl)

url <- "https://www.ebi.ac.uk/proteins/api/variation/P15056?format=json"

# 直接从URL读取JSON数据
variation_list <- fromJSON(url)

feature<-variation_list$features
#---------------------------------------------------------------------------
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

#----------------------------------------------------------------------------------------------
# 提取函数（与上面逻辑类似）
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

# 对每一行计算
disease_vec <- vapply(feature_filtered$association, extract_firstcol_concat_base, FUN.VALUE = character(1), USE.NAMES = FALSE)

# 新增列
feature_filtered$disease <- disease_vec

# 过滤：移除 association 为 NULL 或 disease 为 NA/空 的行
keep <- !vapply(feature_filtered$association, is.null, FUN.VALUE = logical(1)) & !is.na(feature_filtered$disease) & feature_filtered$disease != ""
feature_filtered2 <- feature_filtered[keep, , drop = FALSE]




feature_2<-feature_filtered2[feature_filtered2$consequenceType=="missense",]


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

feature_filtered3<-feature_2[,c(2,3,9,10,11,22,23,24)]

write_xlsx(feature_filtered3,'C:\\毕业论文\\数据\\BRAF\\uniprot_processed.xlsx')