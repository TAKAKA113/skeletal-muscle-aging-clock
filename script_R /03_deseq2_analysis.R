# 03_deseq2_analysis.R
# 目的：
# 1. young / old 用の counts と metadata を読む
# 2. DESeq2 用の形に整える
# 3. young vs old の DE解析を行う
# 4. 結果を CSV で保存する


BiocManager::install("DESeq2")
library(DESeq2)


# 1. データを読む

counts <- read.csv("Data/counts_young_old.csv", check.names = FALSE)
meta   <- read.csv("Data/meta_young_old.csv", check.names = FALSE)


# 2. metadata を必要な列だけにする
#metaのdfからSampleとAge_gruopの列のみを取り出し、上書き

meta <- meta[, c("Sample", "age_group")]


# 3. counts の1列目を gene_id として使う
#DESeqは遺伝子IDを行名にしたい
#そのためgene_idの行名にした列を作り、ただの列のgene_idを消した
#この操作はRNAseqのデータ整形によくでるらしい

rownames(counts) <- counts$gene_id
counts <- counts[, -1]


# 4. metadata の行名を Sample にする

rownames(meta) <- meta$Sample


# 5. Sample 列はもう不要なので外す
#ここもcountと同じで行名に変更した

meta <- meta[, -1, drop = FALSE]




# 6. counts の列順を metadata に合わせる
#DESeqではcountsの語順＝metaの語順でないといけない
#DESeqあるある

counts <- counts[, rownames(meta)]

all(colnames(counts) == rownames(meta))  #countsとmetaのsample順が完全一致か確認

# 7. age_group を factor にする
# young を基準、old を比較対象にする
#DESeqでは何と比べて変化したかをみるため基準を指定する
#levels = c("X", "Y")で指定。Xが基準となる
#factor()は文字列ではなく、グループ分けされてると指定
meta$age_group <- factor(meta$age_group, levels = c("young", "old"))



# 8. counts を整数にする

counts <- round(as.matrix(counts))


# 9. DESeq2 オブジェクトを作る
#ddsにDESeqようにobjectをまとめる

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ age_group
)



# 10. 低カウント遺伝子を軽く除く


dds <- dds[rowSums(counts(dds)) >= 10, ]


# 11. DESeq2 を実行する

dds <- DESeq(dds)


# 12. old vs young の結果を取る

res <- results(dds, contrast = c("age_group", "old", "young"))


# 13. 結果をデータフレームに変換する

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# gene_id を先頭列に移動
res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]


# 14. padj 順に並べる

res_df <- res_df[order(res_df$padj), ]


# 15. 有意遺伝子も作る

sig_res <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)


# 16. 結果を保存する

write.csv(res_df, "Data/deseq2_results_all.csv", row.names = FALSE)
write.csv(sig_res, "Data/deseq2_results_sig.csv", row.names = FALSE)


# 17. 確認用メッセージ

cat("解析完了\n")

cat("\n結果テーブルのサイズ:\n")
print(dim(res_df))

cat("\n有意遺伝子数:\n")
print(nrow(sig_res))

cat("\n上位5件:\n")
print(head(res_df, 5))

cat("\n保存ファイル:\n")
cat("Data/deseq2_results_all.csv\n")
cat("Data/deseq2_results_sig.csv\n")
