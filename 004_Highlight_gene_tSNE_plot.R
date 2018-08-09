require("ggplot2")
require("easyGgplot2")

## Read the two input files: expressioin matrix and cluster file/masterTable with tSNE info
expression <- read.table("/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/02_scMINER/Combined_Total/3_1_Combined_Total_mm10_scMINER_total.txt", header = T, sep = "\t")
colnames(expression)[1] <- "Cell_ID"
mastertable <- read.table("/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/05_cluster_definition/Combined_Total/00_highlight_master_table_Combined_Total_mm10.txt", header = T, sep = "\t")
input <- merge(mastertable, expression, by = "Cell_ID", all = T)
rm(expression); rm(mastertable)

## Pick the gene
target_gene <- "Ptprc" ## Be carefult with the genes: 1) with "-" embedded; 2) begains with numbers
if (!(target_gene %in% colnames(input))) {print(paste0(target_gene, " is not found in the expression matrix."));}

## Calculate the number and percentage of expressed/nonexpressed cells
expressed_cells <- sum(input[, target_gene] > 0); nonexpressed_cells <- sum(input[, target_gene] == 0);
if (expressed_cells <= 0) {print(paste0(target_gene, " is not found in the expression matrix."));}

per.expressed_cells <- paste((round(((expressed_cells/(expressed_cells+nonexpressed_cells))*100),2)),"%");
per.nonexpressed_cells <- paste((round(((nonexpressed_cells/(expressed_cells+nonexpressed_cells))*100),2)),"%");
text_label <- paste0("Expressed Cells: ", expressed_cells, " (", per.expressed_cells, ")\n", "None-expr Cells: ", nonexpressed_cells, " (", per.nonexpressed_cells, ")")

## Highlight the expressed cells
input_expressed <- subset(input, input[, target_gene] > 0);
p1 <- ggplot(input, aes(x = tSNE_X, y = tSNE_Y, colour = input[, target_gene]))
p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_point(size = 0.5, col = "grey") + ggtitle(target_gene, subtitle = text_label) + labs(x = "tSNE_X", y = "tSNE_Y") +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"),
          plot.subtitle = element_text(size = 12, face = "bold", color = "black"),
          axis.title = element_text(size = 15, face = "bold", color = "black"),
          axis.text = element_text(size = 12, face = "bold", color = "black")) +
    geom_point(input_expressed, mapping = aes(x = input_expressed$tSNE_X, y = input_expressed$tSNE_Y), col = "blue", size = 0.5)
##output_file_1 <- paste0(outdir, "/", "01_highlight_byGene_", sample, "_", species, "_", genes[i], "_merged",".pdf")
##ggsave(output_file_1, plot = p1, scale = 1, width = 10, height = 10, units = "in", dpi = 300)

## Violin plot with dots
p2 <- ggplot(input, aes(factor(Cluster), input[, target_gene])) + geom_violin(trim = FALSE, aes(fill = factor(Cluster)))
p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_jitter(height = 0, width = 0.2, size = 0.001) + labs(title = target_gene, x = "Clusters", y = "Log (100K-scaled UMI)") +
    ylim(0.1, NA) +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"),
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))
##output_file_2 <- paste0(outdir, "/", "01_highlight_byGene_", sample, "_", species, "_", genes[i], "_violin_dot",".pdf")
##ggsave(output_file_2, plot = p2, scale = 1, width = 10, height = 10, units = "in", dpi = 300)

## Violin plot without dots
p3 <- ggplot(input, aes(factor(Cluster), input[, target_gene])) + geom_violin(trim = FALSE, aes(fill = factor(Cluster)))
p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title = target_gene, x = "Clusters", y = "Log (100K-scaled UMI)") +
    ylim(0.1, NA) +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"),
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text = element_text(size = 12, face = "bold", hjust = 0.5, color = "black")) + geom_boxplot(width=0.2)
##output_file_3 <- paste0(outdir, "/", "01_highlight_byGene_", sample, "_", species, "_", genes[i], "_violin",".pdf")
##ggsave(output_file_3, plot = p3, scale = 1, width = 10, height = 10, units = "in", dpi = 300)
