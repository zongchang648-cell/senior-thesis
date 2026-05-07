import pandas as pd
import mygene
import os

# ========================
# 1. 参数设置
# ========================
input_file = r"C:/毕业论文/数据处理及绘图/ppi_build/PPIN.xlsx"
out_dir = r"C:/毕业论文/数据处理及绘图/ppi_build"

success_file = os.path.join(out_dir, "PPIN_score700_SYMBOL_success.xlsx")
failed_file  = os.path.join(out_dir, "PPIN_score700_SYMBOL_failed.xlsx")

score_cutoff = 700
species = "human"

# ========================
# 2. 读取并筛选数据
# ========================
df = pd.read_excel(input_file)

df_filtered = df[df["combined_score"] >= score_cutoff].copy()

# ========================
# 3. 提取 ENSP ID（去掉 9606. 前缀）
# ========================
df_filtered["protein1_clean"] = df_filtered["protein1"].str.replace(r"^9606\.", "", regex=True)
df_filtered["protein2_clean"] = df_filtered["protein2"].str.replace(r"^9606\.", "", regex=True)

all_ensp = pd.unique(
    df_filtered[["protein1_clean", "protein2_clean"]].values.ravel()
)

# ========================
# 4. ENSP → SYMBOL 映射
# ========================
mg = mygene.MyGeneInfo()

query_res = mg.querymany(
    all_ensp,
    scopes="ensembl.protein",
    fields="symbol",
    species=species,
    as_dataframe=True
)

ensp2symbol = query_res["symbol"].to_dict()

# ========================
# 5. 映射回原表
# ========================
df_filtered["gene1"] = df_filtered["protein1_clean"].map(ensp2symbol)
df_filtered["gene2"] = df_filtered["protein2_clean"].map(ensp2symbol)

# ========================
# 6. 拆分成功 / 失败结果
# ========================
success_df = df_filtered.dropna(subset=["gene1", "gene2"]).copy()
failed_df  = df_filtered[
    df_filtered["gene1"].isna() | df_filtered["gene2"].isna()
].copy()

# ========================
# 7. 整理输出列
# ========================
success_df = success_df[[
    "gene1", "gene2", "combined_score",
    "protein1", "protein2"
]]

failed_df = failed_df[[
    "protein1", "protein2", "combined_score",
    "gene1", "gene2"
]]

# ========================
# 8. 写出结果
# ========================
success_df.to_excel(success_file, index=False)
failed_df.to_excel(failed_file, index=False)

print("处理完成！")
print(f"成功映射结果: {success_file}")
print(f"映射失败结果: {failed_file}")
print(f"成功条数: {len(success_df)}")
print(f"失败条数: {len(failed_df)}")
