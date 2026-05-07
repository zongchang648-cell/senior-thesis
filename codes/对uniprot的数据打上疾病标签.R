# 加载必要的包
library(readxl)
library(dplyr)
library(stringr)
library(writexl)

# 读取uniprot_processed.xlsx文件
uniprot_data <- read_excel("C:/毕业论文/数据/PTEN/uniprot_processed.xlsx", sheet = 1)

# 读取disease_type.xlsx文件
disease <- read_excel("C:/毕业论文/数据/PTEN/disease_type.xlsx", sheet = 1)

# 创建type列
uniprot_data$type <- ""

# 遍历uniprot_data的每一行
for (i in 1:nrow(uniprot_data)) {
  disease_string <- uniprot_data$disease[i]
  type_values <- c()
  
  # 检查disease数据框的每一行，看是否在disease_string中
  for (j in 1:nrow(disease)) {
    disease_pattern <- disease$the_data2[j]
    if (grepl(disease_pattern, disease_string, fixed = TRUE)) {
      type_values <- c(type_values, disease$Type[j])
    }
  }
  
  # 用"/"连接所有匹配到的type值
  if (length(type_values) > 0) {
    uniprot_data$type[i] <- paste(unique(type_values), collapse = "/")
  } else {
    uniprot_data$type[i] <- "Not"
  }
}

# 对type列进行最终分类
for (i in 1:nrow(uniprot_data)) {
  type_string <- uniprot_data$type[i]
  type_parts <- unlist(str_split(type_string, "/"))
  
  has_NDD <- "NDD" %in% type_parts
  has_Cancer <- "Cancer" %in% type_parts
  
  if (has_NDD && !has_Cancer) {
    uniprot_data$type[i] <- "NDD"
  } else if (has_Cancer && !has_NDD) {
    uniprot_data$type[i] <- "Cancer"
  } else if (has_NDD && has_Cancer) {
    uniprot_data$type[i] <- "OR"
  } else {
    uniprot_data$type[i] <- "Not"
  }
}
uniprot_data<-as.data.frame(uniprot_data)

write_xlsx(uniprot_data,"C:/毕业论文/数据/PTEN/uniprot.xlsx")