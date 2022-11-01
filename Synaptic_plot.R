############### SYNAPTIC PLOT INTERACTION ###############################


syndf <- read.csv("synaptic_genes_CC_Interaction_06.csv")
head(syndf)
dim(syndf)
syndf[order(syndf$padj_int),]


library(ggplot2)
library(ggpubr)

theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

# Initiate a ggplot
b <- ggplot(syndf, aes(x = log2FC_WT, y = log2FC_MUT))

# Basic scatter plot
b + geom_point()

# Change color, shape and size
b + geom_point(color = "#00AFBB", size = 1, shape = 23)

b <- b + geom_point(aes(color = padj_int), size = 4, alpha = 0.5) +
  scale_color_gradientn(colors = c("red2", "white")) +
  theme(legend.position = "right")
b

library(ggrepel)

# Add text to the plot
row.names(syndf) <- (syndf$Symbol)

.labs <- rownames(syndf)
.labs
b + geom_text_repel(aes(label = .labs), size = 2, max.overlaps = 15)
  
b + geom_text_repel(aes(label = .labs), size = 3)


