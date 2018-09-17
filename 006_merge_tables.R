output_merge <- read.table("/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/04_cluster_pathways/01_cluster_vs_rest/Combined_Total/Combined_Total_mm10_cluster_pathways_Cluster_1.txt", header = T, sep = "\t", row.names = 1)

for (i in 2:16) {
    infile = paste0("/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/04_cluster_pathways/01_cluster_vs_rest/Combined_Total/Combined_Total_mm10_cluster_pathways_Cluster_", i, ".txt")
    input <- read.table(infile, header = T, sep = "\t", row.names = 1)
    output_merge <- merge(output_merge, input, by = "row.names", all.x = TRUE)
    row.names(output_merge) <- output_merge$Row.names; output_merge <- output_merge[-1]
}
write.table(output_merge, file = "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/04_cluster_pathways/01_cluster_vs_rest/Combined_Total/Combined_Total_mm10_cluster_pathways_Cluster_Merged.txt", quote = F, row.names = T, col.names = T, sep = "\t")
