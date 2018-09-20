dir1 = "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/05_cluster_definition/Combined_Total"
masterTable1 <- read.table(paste0(dir1, "/00_highlight_master_table_Combined_Total_mm10.txt"), header = T, sep = "\t", row.names = 1)

dir2 = "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/05_cluster_definition/Combined_Total/Subclusters"
masterTable2 <- read.table(paste0(dir2, "/00_highlight_master_table.txt"), header = T, sep = "\t", row.names = 1)
group <- subset(masterTable, Cluster == 18)
group1 <- group[,1:2]
group2 <- merge(group1, masterTable1, by = "row.names", all.x = T)

ggplot(masterTable1, aes(x = tSNE_X, y = tSNE_Y)) +
    geom_point(color="grey", size=0.5) +
    xlab("MICA-1") + ylab("MICA-2") +
    theme(plot.title = element_text(size=25, margin=margin(t=20, b=10)),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14,face="bold")) +
    geom_point(group2, mapping = aes(x = group2$tSNE_X, y = group2$tSNE_Y), col = "blue", size = 0.5)
