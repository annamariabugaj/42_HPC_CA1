
df_wt <- read.csv('ca1_counts_id_wt.csv', sep = ",")
head(df_wt)
meta_wt <- read.csv('ca1_meta_wt.csv', sep = ",")
head(meta_wt)

countData <- df_wt
head(countData)
metaData <- meta_wt
head(metaData)
metaData

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~age, 
                              tidy = TRUE)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

resultsNames(dds)
head(dds)

res_wt <- results(dds, contrast=c("age","adult","pup"))

ix = which.min(res_wt$padj) # most significant

res_wt <- res_wt[order(res_wt$padj),] # sort

summary(res_wt)
head(res_wt)

par(mar=c(7,4,4,4))
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )

write.csv( as.data.frame(res_wt), file="Res_WT_adult_vs_pup.csv" )


schrinkedres_wt <- lfcShrink(dds, "age_adult_vs_pup")
write.csv( as.data.frame(schrinkedres_wt), file="lfcSchrinked_RES_MOD2_age_in_WT_id.csv" )
