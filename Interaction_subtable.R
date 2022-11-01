wt_df <- read.csv("Age_in_WT.csv", sep = ",")
head(wt_df)
class(wt_df)
dim(wt_df)

int <- read.table('Interaction_genelist.txt')
summary(int)
class(int)
dim(int)


subtable_WT <- wt_df[wt_df$Symbol %in% int$V1, ]
head(subtable_WT)
dim(subtable_WT)

write.csv(subtable_WT, 'age_WT_Interaction.csv')


mut_df <- read.csv("Age_in_MUT.csv", sep = ",")
head(mut_df)
class(mut_df)
dim(mut_df)

int <- read.table('Interaction_genelist.txt')
summary(int)
class(int)
dim(int)


subtable_MUT <- mut_df[mut_df$Symbol %in% int$V1, ]
head(subtable_MUT)
dim(subtable_MUT)

write.csv(subtable_MUT, 'age_MUT_Interaction.csv')


######################### SYNAPTIC CC INTERACTION LFC 0.6 #########################

# synaptic genes (CC, lfc > 0.6, synapse, neuron projection, somatodendritic compartment)
# subtable in Interaction genes to see "fodl change" in interaction



pldf <- read.csv("bigtable.csv")
head(pldf)
dim(pldf)

synaptic <- read.table('synaptic_genes_int06.txt')
summary(synaptic)
class(synaptic)
dim(synaptic)

subtable_synaptic <- pldf[pldf$Symbol %in% synaptic$V1, ]
head(subtable_synaptic)
dim(subtable_synaptic)

write.csv(subtable_synaptic, 'synaptic_genes_CC_Interaction_06.csv')


######### SYNAPTIC genes in Adult MUT vs Adult Wt ###########################

addf <- read.csv("Adult_mut_wt_23032021.csv")
head(addf)


subtable_synaptic1 <- addf[addf$Symbol %in% synaptic$V1, ]
head(subtable_synaptic1)
dim(subtable_synaptic1)

write.csv(subtable_synaptic1, 'synaptic_genes_CC_Adult_MUT_vs_WT.csv')

###################### SYNAPTIC genes in Pup MUT vs Pup WT ######################

pddf <- read.csv("Pups_mut_wt_23032021.csv")
head(pddf)

subtable_synaptic2 <- pddf[pddf$Symbol %in% synaptic$V1, ]
head(subtable_synaptic2)
dim(subtable_synaptic2)


write.csv(subtable_synaptic2,'synaptic_genes_CC_Pups_MUT_vs_WT.csv')


####################### SYNAPTIC genes ADULT only ########################

mutdf <- read.csv('0.05_synaptic_genes_CC_Adult_MUT_vs_WT.csv')
head(mutdf)
dim(mutdf)

t1 <- read.table('synaptic_adult_only.txt')
summary(t1)
class(t1)
dim(t1)

subtable_mutdf <- mutdf[mutdf$Symbol %in% t1$V1, ]
head(subtable_mutdf)
dim(subtable_mutdf)

write.csv(subtable_mutdf,'synaptic_genes_CC_ADULT_ONLY.csv')


####################### SYNAPTIC genes PUP ONLY ########################

mutdf2 <- read.csv('synaptic_genes_CC_Pups_MUT_vs_WT.csv')
head(mutdf2)
dim(mutdf2)

t2 <- read.table('synaptic_pup_only.txt')
summary(t2)
class(t2)
dim(t2)

subtable_mutdf2 <- mutdf2[mutdf2$Symbol %in% t2$V1, ]
head(subtable_mutdf2)
dim(subtable_mutdf2)

write.csv(subtable_mutdf2,'synaptic_genes_CC_PUP_ONLY.csv')


####################### SYNAPTIC genes ADULT AND PUP  ########################

mutdf3 <- read.csv('Adult_mut_wt_23032021.csv')
head(mutdf3)
dim(mutdf3)

t3 <- read.table('synaptic_adult-and-pup.txt')
summary(t3)
class(t3)
dim(t3)

subtable_mutdf3 <- mutdf3[mutdf3$Symbol %in% t3$V1, ]
head(subtable_mutdf3)
dim(subtable_mutdf3)

write.csv(subtable_mutdf3,'synaptic_genes_CC_ADULT_AND_PUP_adult.csv')
