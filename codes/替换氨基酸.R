
library(readxl)
library(openxlsx)

# ---------- 用户配置 ----------
input_path  <- "C:\\毕业论文\\数据\\SOS1\\SOS1_hgmd.xlsx"     # 输入文件路径，替换为你的文件名/路径
sheet_index <- 2                # 读取的 sheet（可以用名称或索引）
output_path <- "C:\\毕业论文\\数据\\SOS1\\SOS1_hgmd_2.xlsx"    # 输出文件路径

# ---------- 读入数据 ----------
df <- read_xlsx(path = input_path, sheet = sheet_index)

# ---------- 三字母到单字母的映射（包括 Sec/Pyl 可选项） ----------
map3to1 <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C",
  Gln="Q", Glu="E", Gly="G", His="H", Ile="I",
  Leu="L", Lys="K", Met="M", Phe="F", Pro="P",
  Ser="S", Thr="T", Trp="W", Tyr="Y", Val="V"
)

# ---------- 替换函数：对第一列进行全局替换 ----------
# 使用 perl 正则，(?i) 不区分大小写，(?<![A-Za-z])... (?![A-Za-z]) 保证三字母代码不被字母夹在中间时误替换
for (aa3 in names(map3to1)) {
  aa1 <- map3to1[[aa3]]
  pattern <- paste0("(?i)(?<![A-Za-z])", aa3, "(?![A-Za-z])")
  df[[1]] <- gsub(pattern, aa1, df[[1]], perl = TRUE)
}

# ---------- 写出结果 ----------
write.xlsx(df, output_path, overwrite = TRUE)

cat("已完成：已将第一列中三字母氨基酸替换为单字母，并保存为：", output_path, "\n")
