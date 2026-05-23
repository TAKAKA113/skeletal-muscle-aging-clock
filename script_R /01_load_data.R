# 01_load_data.R
# 目的：
# 1. Dataフォルダの中身を確認する
# 2. raw counts と metadata を読む
# 3. 行数・列数・列名をざっくり確認する

getwd()
list.files()


counts_path <- "Data/raw_counts.csv.gz"
meta_path   <- "Data/Metadata_Complete.csv"

counts <- read.csv(gzfile(counts_path), check.names = FALSE)
meta <- read.csv(meta_path, check.names = FALSE)



cat("\nmetadata のサイズ（行, 列）:\n")
print(dim(meta))

cat("\nmetadata の列名:\n")
print(colnames(meta))

cat("\nmetadata の最初の5行:\n")
print(head(meta, 5))

