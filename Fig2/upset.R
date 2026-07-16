library(dplyr)
library(ComplexHeatmap)

# 读取数据
df <- read.csv("03_merge_all.csv")

# 清洗
df <- df %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(tissue, SYMBOL)

# 按 tissue 提取基因集合
tissue_list <- split(df$SYMBOL, df$tissue)

cat("组织数:", length(tissue_list), "\n")
for (nm in names(tissue_list)) {
  cat(sprintf("  %s: %d genes\n", nm, length(tissue_list[[nm]])))
}

# 构建 UpSet 组合矩阵
outdir <- "F:/屠宰性状-转录组分析/06_基因合并"
m <- make_comb_mat(tissue_list)

pdf(file.path(outdir, "Upset_Tissue_ALL.pdf"), width = 16, height = 8)

UpSet(m,
      comb_order = order(comb_size(m), decreasing = TRUE),
      pt_size          = unit(2, "mm"),          # ⭐ 点变小
      top_annotation   = upset_top_annotation(m, add_numbers = FALSE,    # ⭐ 不标数字
                                              height = unit(5, "cm")),   # ⭐ 柱状图更高
      right_annotation = upset_right_annotation(m, add_numbers = FALSE,  # ⭐ 不标数字
                                                width = unit(4, "cm")),  # ⭐ 柱状图更宽
      row_names_gp     = gpar(fontsize = 10))

dev.off()

cat("Upset PDF 已保存\n")
