# compute_prs_from_modules.py
"""
Usage:
  1) 修改下方 INPUT_PATHS 中的文件路径为你的本地路径（已按你给出路径预设）
  2) 运行：python compute_prs_from_modules.py
生成文件（主要）：
  - filtered_edges.csv
  - prs_matrix.npy  (原始 PRS 矩阵，numpy 二进制)
  - prs_matrix_normalized.npy (按节点特异方差归一化后的 PRS)
  - prs_matrix_normalized.csv  (文本，可能很大)
  - prs_plots_normalized.png   (对归一化后的 PRS 作图)
  - responses.npy / responses.csv (基于归一化后的 PRS 计算的响应矩阵)
"""

import os
import numpy as np
import pandas as pd
import logging
import matplotlib.pyplot as plt

# ProDy imports
from prody import GNM  # 如果 import 失败，请确认 prody 已正确安装

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

# -----------------------
# 路径配置（按你的说明预设）
# -----------------------
INPUT_PATHS = {
    'modules_xlsx': r"C:\毕业论文\数据处理及绘图\ppi_build\modules_output\第一次运行\modules_fast_greedy.xlsx",
    'corenet_xlsx': r"C:\毕业论文\数据处理及绘图\ppi_build\top10%_CoreNet.xlsx",
    'perturb_xlsx': r"C:\毕业论文\数据\PCA_compressed_data.xlsx",
}

OUT_DIR = r"C:\毕业论文\数据处理及绘图\ppi_build\python_prs_output"
os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------
# 1) 读取 module gene 列表（sheet2），保留顺序
# -----------------------
logging.info("Reading module gene list from %s", INPUT_PATHS['modules_xlsx'])
modules_df = pd.read_excel(INPUT_PATHS['modules_xlsx'], sheet_name=1)  # sheet2 -> sheet_name=1
# 假设列名为 'name'（你给的样例），否则请相应修改
if 'name' not in modules_df.columns:
    raise ValueError("modules_xlsx sheet2 中找不到 'name' 列，请检查列名。")
module_genes = modules_df['name'].astype(str).tolist()
n_nodes = len(module_genes)
logging.info("Loaded %d module genes (nodes).", n_nodes)

# 建立 gene -> index 映射，保持文件中顺序
node_to_idx = {g: i for i, g in enumerate(module_genes)}
idx_to_node = {i: g for g, i in node_to_idx.items()}

# -----------------------
# 2) 读取 core net 边表并过滤：仅保留两端都在 module_genes 中的边
# -----------------------
logging.info("Reading core network edges from %s", INPUT_PATHS['corenet_xlsx'])
edges_df = pd.read_excel(INPUT_PATHS['corenet_xlsx'], sheet_name=0)  # sheet1
# 期望列名为 'gene1','gene2'
if not {'gene1', 'gene2'}.issubset(edges_df.columns):
    raise ValueError("corenet_xlsx sheet1 中找不到 'gene1' 和 'gene2' 列，请检查列名。")

# 转成字符串并去两端空白
edges_df['gene1'] = edges_df['gene1'].astype(str).str.strip()
edges_df['gene2'] = edges_df['gene2'].astype(str).str.strip()

logging.info("Filtering edges: keep only edges where both genes are in module list...")
mask_both_in = edges_df['gene1'].isin(node_to_idx) & edges_df['gene2'].isin(node_to_idx)
filtered_edges = edges_df[mask_both_in].copy()
logging.info("Filtered edges: %d / %d kept.", filtered_edges.shape[0], edges_df.shape[0])
filtered_edges.to_csv(os.path.join(OUT_DIR, "filtered_edges.csv"), index=False)

# -----------------------
# 3) 构建邻接矩阵 A (n_nodes x n_nodes)，无向、无权重（0/1）
# -----------------------
logging.info("Building adjacency matrix (unweighted, undirected)...")
A = np.zeros((n_nodes, n_nodes), dtype=float)
for _, row in filtered_edges.iterrows():
    g1 = row['gene1']
    g2 = row['gene2']
    i = node_to_idx[g1]
    j = node_to_idx[g2]
    if i == j:
        continue
    A[i, j] = 1.0
    A[j, i] = 1.0

# 可选：检查是否有孤立节点（degree=0）
degrees = A.sum(axis=1).astype(int)
n_isolated = int(np.sum(degrees == 0))
logging.info("Degrees computed. isolated nodes (degree=0): %d", n_isolated)

# -----------------------
# 4) 构建 Kirchhoff (Laplacian) 矩阵：K = diag(degree) - A
#    并把它交给 ProDy 的 GNM
# -----------------------
logging.info("Building Kirchhoff (Laplacian) matrix and creating GNM model...")
K = np.diag(degrees) - A  # NxN

# ProDy GNM：setKirchhoff 接受一个 numpy 矩阵
gnm = GNM("modules_GNM")
gnm.setKirchhoff(K)

# 计算所有非零模态（n_modes=None -> all），默认会忽略零模式（平移模式）
logging.info("Calculating modes (this may take a while for large matrices)...")
try:
    gnm.calcModes(n_modes=None)  # 计算所有模式
except Exception as e:
    logging.warning("gnm.calcModes(n_modes=None) 报错，尝试 calcModes(n_modes = n_nodes - 1). Error: %s", str(e))
    gnm.calcModes(n_modes=max(1, n_nodes - 1))

# -----------------------
# 5) 计算协方差矩阵（GNM 下的 C），把它视作 PRS 矩阵（原始）
# -----------------------
logging.info("Extracting covariance matrix (interpret as PRS matrix)...")
C = gnm.getCovariance()  # NxN
# C 可能是 numpy.ndarray 或 masked array
PRS = np.array(C, dtype=float)  # 强制转成 ndarray

# 保存 原始 PRS
prs_npy = os.path.join(OUT_DIR, "prs_matrix.npy")
prs_csv = os.path.join(OUT_DIR, "prs_matrix.csv")
np.save(prs_npy, PRS)
logging.info("Saved raw PRS matrix to %s", prs_npy)
try:
    pd.DataFrame(PRS, index=module_genes, columns=module_genes).to_csv(prs_csv)
    logging.info("Saved raw PRS matrix to %s", prs_csv)
except Exception as e:
    logging.warning("保存 raw PRS as CSV 失败（可能内存/磁盘限制）：%s", str(e))

# -----------------------
# 5b) 对 PRS 做归一化：每列除以其节点特异方差（diag(PRS)）
#      PRS_norm[:, j] = PRS[:, j] / diag_vars[j]
# -----------------------
logging.info("Normalizing PRS by node-specific variances (diagonal of PRS)...")
diag_vars = np.diag(PRS).astype(float)  # n_nodes 长度的一维数组

# 处理非正或为零的对角元，避免除零
zero_or_nonpos = np.where(diag_vars <= 0)[0]
if zero_or_nonpos.size > 0:
    logging.warning("发现 %d 个对角方差 <= 0（将用小常数替代以避免除零）：%s", zero_or_nonpos.size, zero_or_nonpos[:20].tolist())
    # 用很小的正数替代（保守策略），并记录
    eps = 1e-12
    diag_vars_safe = diag_vars.copy()
    diag_vars_safe[diag_vars_safe <= 0] = eps
else:
    diag_vars_safe = diag_vars

# 按列除以对应的 diag_vars
PRS_norm = PRS / diag_vars_safe[np.newaxis, :]  # 广播：每列除以对应 diag_vars_safe[j]

# 保存归一化后的 PRS
prs_npy_norm = os.path.join(OUT_DIR, "prs_matrix_normalized.npy")
prs_csv_norm = os.path.join(OUT_DIR, "prs_matrix_normalized.csv")
np.save(prs_npy_norm, PRS_norm)
logging.info("Saved normalized PRS matrix to %s", prs_npy_norm)
try:
    pd.DataFrame(PRS_norm, index=module_genes, columns=module_genes).to_csv(prs_csv_norm)
    logging.info("Saved normalized PRS matrix to %s", prs_csv_norm)
except Exception as e:
    logging.warning("保存 normalized PRS as CSV 失败（可能内存/磁盘限制）：%s", str(e))

# -----------------------
# 6) 可视化：使用归一化后的 PRS_norm 绘制热图 + effectiveness/sensitivity 曲线
#    effectiveness = 列求和 (对每个扰动源的总效应)
#    sensitivity   = 行求和 (每个节点对所有扰动的总敏感度)
# -----------------------
logging.info("Plotting normalized PRS heatmap and profiles...")
effectiveness = PRS_norm.sum(axis=0)  # shape (n_nodes,)
sensitivity = PRS_norm.sum(axis=1)   # shape (n_nodes,)

fig = plt.figure(figsize=(14, 10))

# heatmap (PRS_norm matrix)
ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=1, rowspan=2)
im = ax1.imshow(PRS_norm, aspect='auto', interpolation='nearest')
ax1.set_title('PRS matrix (normalized by node-specific variance)')
ax1.set_xlabel('perturbation node index')
ax1.set_ylabel('response node index')
plt.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)

# effectiveness (columns sum)
ax2 = plt.subplot2grid((3, 2), (0, 1))
ax2.plot(effectiveness)
ax2.set_title('Effectiveness (sum over responses for each perturbation source) - normalized')
ax2.set_xlabel('perturbation node index')
ax2.set_ylabel('effectiveness (sum)')

# sensitivity (rows sum)
ax3 = plt.subplot2grid((3, 2), (1, 1))
ax3.plot(sensitivity)
ax3.set_title('Sensitivity (sum over perturbation sources for each response node) - normalized')
ax3.set_xlabel('response node index')
ax3.set_ylabel('sensitivity (sum)')

plt.tight_layout()
plotfile = os.path.join(OUT_DIR, "prs_plots_normalized.png")
plt.savefig(plotfile, dpi=300)
plt.close(fig)
logging.info("Saved normalized PRS plots to %s", plotfile)

# -----------------------
# 7) 读取扰动文件，逐条施加扰动并记录响应（responses matrix）
#    规则：每条记录是独立扰动，force 就是 total_score，扰动只施给单一节点（Gene）
#    注意：使用归一化后的 PRS_norm 来计算响应
# -----------------------
logging.info("Reading perturbation list from %s", INPUT_PATHS['perturb_xlsx'])
perturb_df = pd.read_excel(INPUT_PATHS['perturb_xlsx'], sheet_name=0)
if 'Gene' not in perturb_df.columns or 'total_score' not in perturb_df.columns:
    raise ValueError("perturb_xlsx 中必须包含 'Gene' 和 'total_score' 列，请检查。")

perturb_genes = perturb_df['Gene'].astype(str).str.strip().tolist()
perturb_forces = perturb_df['total_score'].astype(float).tolist()
n_pert = len(perturb_genes)
logging.info("Loaded %d perturbations.", n_pert)

# 为保证形状一致（n_nodes x n_perturbations），我们为找不到映射的 gene 填 0 列，并记录
responses = np.zeros((n_nodes, n_pert), dtype=float)
missing_genes = []
for col, (g, f) in enumerate(zip(perturb_genes, perturb_forces)):
    if g in node_to_idx:
        idx = node_to_idx[g]
        # 单位扰动对应的响应是 PRS_norm[:, idx]， 加权扰动乘以 f
        responses[:, col] = PRS_norm[:, idx] * f
    else:
        missing_genes.append((col, g))
        # responses[:, col] remains zeros

if missing_genes:
    logging.warning("以下扰动的基因在 module 节点中找不到（它们对应的响应列被置为 0），共 %d 个：", len(missing_genes))
    for col, g in missing_genes[:20]:
        logging.warning("  col %d: %s", col, g)
    if len(missing_genes) > 20:
        logging.warning("  ... (另外 %d 个未列出)", len(missing_genes) - 20)

# 保存 responses
resp_npy = os.path.join(OUT_DIR, "responses.npy")
resp_csv = os.path.join(OUT_DIR, "responses.csv")
np.save(resp_npy, responses)
logging.info("Saved responses matrix to %s", resp_npy)
try:
    pd.DataFrame(responses, index=module_genes).to_csv(resp_csv)
    logging.info("Saved responses matrix to %s", resp_csv)
except Exception as e:
    logging.warning("保存 responses CSV 失败（可能太大）：%s", str(e))

logging.info("All done. Outputs are in %s", OUT_DIR)