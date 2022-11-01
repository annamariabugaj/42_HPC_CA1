if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(knitr)
df <-read.csv('ca1_counts_id.csv', header = TRUE, sep = ",")
head(df)


meta <- read.csv('ca1_meta.csv', header = TRUE, sep = ',')
head(meta)

countData <- df
head(countData)
metaData <- meta
head(metaData)


########## age effect specific for genotype GROUP method ############

dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData, 
                              design = ~age, tidy = TRUE)

dds$genotype
dds$genotype = relevel( as.factor(dds$genotype), "wt")
dds$genotype

dds$age
dds$age = relevel( as.factor(dds$age), "adult")
dds$age
dds <- estimateSizeFactors(dds)

dds$group <- factor(paste0(dds$genotype, dds$age))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)
results(dds, contrast=c("group", "wtpup", "mutpup"))
results(dds, contrast=c("group", "wtadult", "mutadult"))
results(dds, contrast=c("group", "wtpup", "mutadult"))

data <- plotCounts(dds, 
           gene = "Tuba1a", 
           intgroup = c("genotype", "age"),
           returnData = TRUE)

ggplot(data, aes(x=interaction(genotype, age), y=count, color=genotype, shape=age)) + 
  geom_point(size=3) + 
  scale_y_log10()

################# GOOD multifactorial design #################

#Here we have two genotypes, wild-type (wt), and mutant (mut). 
#Two conditions, genotype (genotype) and age (age)
#We are interested in gene expression of both wild-type and mutant across age.
#We are also interested in the differences in gene expression across age between genotypes, 
#which is captured by the interaction term in linear models.


dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~genotype+age+genotype:age, 
                              tidy = TRUE)
dds
dds$genotype
dds$genotype = relevel( dds$genotype, "wt")
dds$genotype

dds$age
dds$age = relevel( dds$age, "adult")
dds$age

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

resultsNames(dds)
head(dds)
res <- results(dds)
summary(res)


#Below, we are going to use the combination of the different results 
#(“genotype_mut_vs_wt”, “age_pup_vs_adult”, “genotypemut.agepup” ) 
#to derive biologically meaningful comparisons.

########### The effect of age in WT - the main effect ###################
#This is for WT, pup compared with adult. 
#Note that WT is not mentioned, because it is the reference level. 
#In other words, this is the difference between WT pups(11:14) and adults(5:7)

res1 = results(dds, contrast=c("age","pup","adult"))
ix = which.min(res1$padj) # most significant
res1 <- res1[order(res1$padj),] # sort
res1 <- res1[order(res1$log2FoldChange)]

kable(res1[1:5,-(3:4)])
head(res1)
d<- cbind(rownames(res1), data.frame(res1, row.names=NULL))
head(d)

summary(res1)


par(mar=c(7,4,4,4))
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
write.csv( as.data.frame(res1), file="Results_age_in_wt_id.csv" )

EnhancedVolcano(d,
                lab = d$`rownames(res1)`,
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0, 400),
                pCutoff = 0.05,
                title = 'WT pup vs adult',
                subtitle = 'Padj = 0.05 log2FC = 1',
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                legendLabels = c('NS', bquote(Log[2]~FC),
                                 'padj', bquote(padj~and~log[2]~FC)),
                FCcutoff = 0.3,
                legendLabSize = 10,
                legendIconSize = 3.0,
                pointSize = 3.5,
                labSize = 2.5,
                col = c("grey30", "grey30", "royalblue", "red2"))

#################The effect of age in MUTANT ###############################
#This is, by definition, the main effect plus the interaction term 
#(the extra age effect in genotype Mutant compared to genotype WT).

res2 <- results(dds, list( c("age_pup_vs_adult","genotypemut.agepup") ))
ix = which.min(res2$padj) # most significant
res2 <- res2[order(res2$padj),] # sort
kable(res2[1:5,-(3:4)])

#This measures the effect of age in mutant. 
#In other words, samples 1:4 compared with samples 8:10.

#Here we show the most significant gene, which is upregulated expressed in samples 1:4.


par(mar=c(7,4,4,4))
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
write.csv( as.data.frame(res2), file="Results_age_in_mut_id.csv" )
summary(res2)

##############difference between mutant and wild-type independent of age########

res3 = results(dds, contrast=c("genotype","mut","wt"))
ix = which.min(res3$padj) # most significant
res3 <- res3[order(res3$padj),] # sort
kable(res3[1:5,-(3:4)])

par(mar=c(7,4,4,4))
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
write.csv( as.data.frame(res3), file="Results_adult_MUT_vs_WT_id.csv" )
summary(res3)


############ With age, what is the difference between mutant and wild-type ######
res4 = results(dds, list( c("genotype_mut_vs_wt","genotypemut.agepup") ))
ix = which.min(res4$padj) # most significant
res4 <- res4[order(res4$padj),] # sort
kable(res4[1:5,-(3:4)])

summary(res4)

#This gives us the difference between genotype MUT and WT, in age pup.
#this is the sampless 1:4 compared with samples 11:14.

#Here we show the most significant gene.
par(mar=c(7,4,4,4))
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )
write.csv( as.data.frame(res4), file="Results_pup_MUT_vs_WT_id.csv" )

########The different response in genotypes (interaction term)##########
#Is the effect of age different across genotypes? 
#This is the interaction term.

res5 = results(dds, name="genotypemut.agepup")
ix = which.min(res5$padj) # most significant
res5 <- res5[order(res5$padj),] # sort
kable(res5[1:5,-(3:4)])
summary(res5)
barplot(assay(dds)[ix,],las=2, main=rownames(dds)[ ix  ]  )

#In wild-type, this gene is downregulated in pup, while for the mutant, it is upregulated.
write.csv( as.data.frame(res5), file="Results_interaction_id.csv" )


######################### VISUALIZATION ###################

ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("genotype","age")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#################### SAMPLE DISTANCE ####################

library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$age, vsd$genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("genotype", "age"))

pcaData <- plotPCA(vsd, intgroup=c("age", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=age, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

