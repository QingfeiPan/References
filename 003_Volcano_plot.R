## This script is made for plotting volcano figure in DE/DA analysis

## Configure Environment
library(dplyr)
library(ggplot2)
library(ggrepel)

## Read the DE/DA output, including GeneSymbol/GeneSet, G1_Expr/Activity, G2_Expr/Activity, Foldchang, Pvalue and FDR
dir <- "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/01_Liver_cancer/03_cluster_genes/02_celltypes/01_Tumor_Cells/lowQDTCs_excluded/11_BulkEarly_subcluster/BulkEarly_DE"
input <- read.table(paste0(dir, "/01_BulkEarly_DE_Gene_Left_VS_Right.txt"), header = T, sep = "\t")

## Formatting
input <- input[,c(1,4,5)] # Pick the GeneSymbol/GeneSet, Foldchange and Pvalue
row.names(input) <- input$GeneSymbol
input[,3] <- ifelse(input[,3]==0, 1e-300, input[,3]) # Avoid the error by log-transform of '0'
input <- mutate(input, Sig=ifelse(input[,3]<0.01, "Yes", "No")) # Add a new column to save significance

## Get the subgrounps
all <- input
up <- subset(input, input[,3]<=0.01 & input[,2]>=log2(1.5))
down <- subset(input, input[,3]<=0.01 & input[,2]<=log2(2/3))
marked <- filter(input, input[,3]<=0.01 & (input[,2]>=1.5 | input[,2]<=-0.7)) # Select the genes to mark names on

## Visulization
p <- ggplot(all, aes(x = all[,2], y = -log10(all[,3])))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white", color = "black", size = 1)) +
    geom_point(size = 0.8, colour = "grey") + labs(x = "Log2 (Fold Change)", y = "-Log10 (P-value)") +
    xlim(-3,3) + ylim(0,310) + # To zoon in specific area
    theme(legend.position = "none", legend.text = element_text(size=12), legend.title = element_text(size=14,face="bold"),
          axis.title = element_text(size = 15, face = "bold", colour = "black"), axis.text = element_text(size = 12, face = "bold", colour = "black")) +
    geom_point(up, mapping = aes(x = up[,2], y = -log10(up[,3])), col = "red", size = 0.8) +
    geom_point(down, mapping = aes(x = down[,2], y = -log10(down[,3])), col = "blue", size = 0.8) +
    geom_hline(aes(yintercept=2), linetype = "dashed", colour = "black") +
    geom_vline(aes(xintercept=log2(3/2)), linetype = "dashed", colour = "black") + 
    geom_vline(aes(xintercept=log2(2/3)), linetype = "dashed", colour = "black") +
    geom_text_repel(data=marked, aes(label=marked$GeneSymbol, x = marked[,2], y = -log10(marked[,3])), size = 3) +
    annotate("text", x=2, y=15, label=paste0("Up-regulated: ", dim(up)[1]), colour="red", fontface="bold", size=5) +
    annotate("text", x=-2, y=15, label=paste0("Down-regulated: ", dim(down)[1]), colour="blue", fontface="bold", size=5)
ggsave(filename = paste0(dir, "/01_BulkEarly_DE_Gene_Left_VS_Right.pdf"), p, width = 10, height = 10, units = "in")
