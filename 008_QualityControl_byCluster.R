dir = "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/05_cluster_definition/Combined_Total"
masterTable <- read.table(paste0(dir, "/00_highlight_master_table_Combined_Total_mm10.txt"), header = T, sep = "\t", row.names = 1)

p <- ggplot(masterTable, aes(factor(Cluster), nGene)) + geom_violin(trim = FALSE, aes(fill = factor(Cluster)))
p1 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_jitter(height = 0, width = 0.2, size = 0.001) +
    labs(x = "Cluster", y = "Number of Genes") + geom_boxplot(width=0.2) +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"), 
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))
ggsave(filename = paste0(dir, "/02_highlight_byQuality_Gene_count_violin_dot.pdf"), p1, width = 10, height = 10, units = "in")

p2 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        labs(x = "Cluster", y = "Number of Genes") + geom_boxplot(width=0.2) +
        theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"), 
              axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
              axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
              axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))
ggsave(filename = paste0(dir, "/02_highlight_byQuality_Gene_count_violin.pdf"), p2, width = 10, height = 10, units = "in")       

p <- ggplot(masterTable, aes(factor(Cluster), nUMI)) + geom_violin(trim = FALSE, aes(fill = factor(Cluster)))
p1 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_jitter(height = 0, width = 0.2, size = 0.001) +
    labs(x = "Cluster", y = "Number of UMIs") + geom_boxplot(width=0.2) +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"), 
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))
ggsave(filename = paste0(dir, "/01_highlight_byQuality_UMI_count_violin_dot.pdf"), p1, width = 10, height = 10, units = "in")

p2 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x = "Cluster", y = "Number of UMIs") + geom_boxplot(width=0.2) +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"), 
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))
ggsave(filename = paste0(dir, "/01_highlight_byQuality_UMI_count_violin.pdf"), p2, width = 10, height = 10, units = "in")

p <- ggplot(masterTable, aes(factor(Cluster), pMito)) + geom_violin(trim = FALSE, aes(fill = factor(Cluster)))
p1 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_jitter(height = 0, width = 0.2, size = 0.001) +
    labs(x = "Cluster", y = "Percentage of Mitochondrial Genes") + geom_boxplot(width=0.2) +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"), 
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))
ggsave(filename = paste0(dir, "/03_highlight_byQuality_pMito_violin_dot.pdf"), p1, width = 10, height = 10, units = "in")

p2 <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(x = "Cluster", y = "Percentage of Mitochondrial Genes") + geom_boxplot(width=0.2) +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"), 
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))
ggsave(filename = paste0(dir, "/03_highlight_byQuality_pMito_violin.pdf"), p2, width = 10, height = 10, units = "in")
