input <- read.table("/Volumes/project_space/yu3grp/scRNASeq/yu3grp/metastasis/01_Liver_cancer/03_cluster_genes/02_celltypes/02_Microenvironment_Cells/RBCs_excluded/05_fibroblast/00_fibroblast_expressionMatrix_clean.txt",
                    header = T, sep = "\t", row.names = 1)
input_t <- t(input)
df <- data.frame(row.names = row.names(input_t), Hif1a=input_t$Hif1a)
for (i in 1:nrow(df)) {
    df[i,2] <- strsplit(row.names(df)[i], "_")
}
names(df)[2] <- "Sample"
head(df)

p <- ggplot(df, aes(factor(Sample), df[, "Hif1a"])) + geom_violin(trim = FALSE, aes(fill = factor(Sample)))
p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_jitter(height = 0, width = 0.2, size = 0.001) + labs(title = "Hif1a", x = "Sample", y = "Log (100K-scaled UMI)") +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"), 
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))+
    ylim(0.1,6) + geom_boxplot(width=0.2)
    stat_summary(fun.y=median, geom="point", size=2, color="red")
