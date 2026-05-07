library(ggplot2)

# 读取数据
file_path <- "C:/毕业论文/数据/core_subnet_topo_features.csv"
df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# 取 Degree 列
deg <- as.numeric(df$Degree)

# 统计每个 Degree 值对应的蛋白数量
deg_freq <- as.data.frame(table(deg), stringsAsFactors = FALSE)
colnames(deg_freq) <- c("Degree", "Count")
deg_freq$Degree <- as.numeric(as.character(deg_freq$Degree))
deg_freq$Count <- as.numeric(deg_freq$Count)

# 绘图
p <- ggplot(deg_freq, aes(x = Degree, y = Count)) +
  geom_point(size = 2, color = "#2C7FB8") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Degree distribution",
    x = "Number of interaction partners",
    y = "Number of proteins"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)


library(ggplot2)

# 读取数据
file_path <- "C:/毕业论文/数据/NCCN_topo_features.csv"
df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# 取 Degree 列
deg <- as.numeric(df$Degree)

# 统计每个 Degree 值对应的蛋白数量
deg_freq <- as.data.frame(table(deg), stringsAsFactors = FALSE)
colnames(deg_freq) <- c("Degree", "Count")
deg_freq$Degree <- as.numeric(as.character(deg_freq$Degree))
deg_freq$Count <- as.numeric(deg_freq$Count)

# 绘图
p <- ggplot(deg_freq, aes(x = Degree, y = Count)) +
  geom_point(size = 2, color = "#2C7FB8") +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Degree distribution",
    x = "Number of interaction partners",
    y = "Number of proteins"
  ) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p)