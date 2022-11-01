
pldf <- read.csv("bigtable.csv")
head(pldf)
dim(pldf)
pldf[order(pldf$padj_int),]


library(ggplot2)
library(ggpubr)

theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

# Initiate a ggplot
b <- ggplot(pldf, aes(x = log2FoldChange_age_in_WT, y = log2FoldChange_age_in_MUT))

# Basic scatter plot
b + geom_point()

# Change color, shape and size
b + geom_point(color = "#00AFBB", size = 1, shape = 23)

b + geom_point(aes(color = padj_int), size = 1, alpha = 0.5) +
  scale_color_gradientn(colors = c("red","grey"),  values = c(0.00:0.05, 0.05:1.00), breaks = c(0.00, 0.05, 0.25, 0.50, 0.75,1.00)) +
  theme(legend.position = "right")

b + geom_point(aes(color = padj_int), size = 1, alpha = 1) +
  scale_color_gradientn(colors = c("red","red", "grey"),  values = c(0.00:0.05, 0.05:1.00)) +
  theme(legend.position = "right")

b + geom_point(aes(color = padj_int), size = 1, alpha = 1) +
  scale_color_gradientn(colors = c("red", "grey"), values = c(0.00:0.05, 0.05:1.00), breaks = c(0.00, 0.05, 0.25, 0.50, 0.75,1.00)) +
  theme(legend.position = "right")





#library(ggrepel)

# Add text to the plot
#.labs <- df$Symbol
#b + geom_point(aes(color = cyl)) +
#  geom_text_repel(aes(label = .labs,  color = cyl), size = 3)+
#  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
