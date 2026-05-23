# 02_prepare_data.R
# 目的：
# 1. counts と metadata を読む
# 2. young / old を作る
# 3. metadata と counts のサンプルをそろえる
# 4. DESeq2用の準備データを保存する



# 1. データを読む

counts_path <- "Data/raw_counts.csv.gz"
meta_path   <- "Data/Metadata_Complete.csv"

counts <- read.csv(gzfile(counts_path), check.names = FALSE)
meta   <- read.csv(meta_path, check.names = FALSE)




# 2. 列名をわかりやすくする

colnames(counts)[1] <- "gene_id"




# 3. young / old を作る

meta_subset <- meta[meta$Age < 50 | meta$Age > 65, ]

meta_subset$age_group <- ifelse(meta_subset$Age < 50, "young", "old")




# 4. counts と metadata で共通サンプルを探す

count_samples <- colnames(counts)[-1]
meta_samples  <- meta_subset$Sample

common_samples <- intersect(count_samples, meta_samples)



# 5. metadata を共通サンプルだけにする

meta_subset <- meta_subset[meta_subset$Sample %in% common_samples, ]




# 6. metadata の順番に counts を並べる

meta_subset <- meta_subset[order(match(meta_subset$Sample, common_samples)), ]

counts_subset <- counts[, c("gene_id", meta_subset$Sample)]


# 7. 確認


cat("metadata のサイズ:\n")
print(dim(meta_subset))

cat("\ncounts のサイズ:\n")
print(dim(counts_subset))

cat("\nage_group の人数:\n")
print(table(meta_subset$age_group))

cat("\nmetadata の最初の5行:\n")
print(head(meta_subset, 5))

cat("\ncounts の最初の列名（10個まで）:\n")
print(colnames(counts_subset)[1:min(10, ncol(counts_subset))])



# 8. 保存


write.csv(meta_subset, "Data/meta_young_old.csv", row.names = FALSE)
write.csv(counts_subset, "Data/counts_young_old.csv", row.names = FALSE)

cat("\n保存完了:\n")
cat("Data/meta_young_old.csv\n")
cat("Data/counts_young_old.csv\n")


#最初に counts と metadata を読みます。
#そのあと、counts の1列目の名前を gene_id に変えて、わかりやすくします。
#次に metadata の Age を見て、50歳未満か65歳超かでサンプルを残し、age_group という新しい列を作ります。
#そのあと、metadata の Sample と counts の列名を比べて、両方に存在するサンプルだけを残します。
#最後に、その揃えたデータを meta_young_old.csv と counts_young_old.csv として Data に保存します。