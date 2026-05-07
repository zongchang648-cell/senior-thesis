library(readxl)
library(writexl)
library(dplyr)
mutation<-read_xlsx('C:\\毕业论文\\数据\\PTPN11\\PTPN11_origin.xlsx',sheet = 1)
mutation<-mutation[,1:3]


mutation_agg <- mutation %>%
  group_by(Mutation) %>%
  summarise(
    Disease = {
      vals <- unique(na.omit(Disease))
      paste(vals, collapse = "/")
    },
    Type = {
      vals <- unique(na.omit(Type))
      # 只有 Cancer
      if (length(vals) > 0 && all(vals == "Cancer")) {
        "Cancer"
        # 只有 NDD
      } else if (length(vals) > 0 && all(vals == "NDD")) {
        "NDD"
        # 同时包含 Cancer 和 NDD（不管是否还有其他值） -> "OR"
      } else if ("Cancer" %in% vals && "NDD" %in% vals) {
        "OR"
        # 其它情况：把去重后的值用 "/" 拼接（保持 vals 的顺序）
      } else if (length(vals) == 0) {
        NA_character_
      } else {
        paste(vals, collapse = "/")
      }
    }
  ) %>%
  ungroup()
write_xlsx(mutation_agg,'C:\\毕业论文\\数据\\PTPN11\\PTPN11-new.xlsx')