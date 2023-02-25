library(infercnv)

obj <- readRDS("/research_jude/rgs01_jude/groups/yu3grp/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_ALL/by_Seurat/01_seuratObj/01_seuratObj_combinedALL.s3_clustered.rds")
dir <- "/research_jude/rgs01_jude/groups/yu3grp/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_ALL/by_Seurat/03_inferCNV"

# 1. prepare raw count matrix
rawCounts <- as.matrix(obj@assays$RNA@counts)
output <- data.frame(GeneSymbol = row.names(rawCounts), rawCounts)
write.table(output, file = paste0(dir, "/00_inferCNV.input_rawCounts.txt"), col.names = T, row.names = F, sep = "\t", quote = F)

# 2. prepare the annotation file
fd <- data.frame(cellID = gsub("-", ".", row.names(obj@meta.data)), clusterID = as.numeric(obj@meta.data$seurat_clusters)); head(fd)
table(fd$cellID %in% colnames(output))
write.table(fd, file = paste0(dir, "/00_inferCNV.input_cellAnnotation.txt"), col.names = F, row.names = F, sep = "\t", quote = F)

# 3. prepare gene coordinate file
system("cp /research_jude/rgs01_jude/groups/yu3grp/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_TYR/by_Seurat/03_inferCNV/00_inferCNV.input_geneAnnotation.txt /research_jude/rgs01_jude/groups/yu3grp/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_ALL/by_Seurat/03_inferCNV/00_inferCNV.input_geneAnnotation.txt")

# 4. create inferCNV object
dir <- "/Volumes/groups/yu3grp/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_ALL/by_Seurat/03_inferCNV"
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = paste0(dir, "/00_inferCNV.input_rawCounts.txt"),
                                    annotations_file = paste0(dir, "/00_inferCNV.input_cellAnnotation.txt"),
                                    delim="\t",
                                    max_cells_per_group = 1000,
                                    gene_order_file = paste0(dir, "/00_inferCNV.input_geneAnnotation.txt"),
                                    ref_group_names=c("15","16")) 
# 5. run inferCNV
options(scipen = 100) # better to add this when the cell count is large
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=dir, 
                             num_threads = 1,
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE)

## 6. plot single gene CNV: SMARCB1 or SMARCA4
cellAnno <- read.table("/Volumes/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_ALL/by_Seurat/03_inferCNV/00_inferCNV.input_cellAnnotation.txt", header = F, sep = "\t", quote = "", stringsAsFactors = F)
infercnv <- read.table("/Volumes/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_ALL/by_Seurat/03_inferCNV/infercnv.observations.txt", header = T, sep = "", quote = "", stringsAsFactors = F)

row.names(infercnv) <- gsub('\\"(.+)\\"', "\\1", row.names(infercnv))
colnames(infercnv) <- gsub("X\\.(.+)\\.", "\\1", colnames(infercnv))

table(colnames(infercnv) %in% cellAnno$V1)
table(cellAnno$V1 %in% colnames(infercnv))

master <- merge(cellAnno, t(infercnv), by.x = "V1", by.y = "row.names", all.y = T);
if (!"SMARCB1" %in% colnames(master)) { master$SMARCB1 <- 1 }
if (!"SMARCA4" %in% colnames(master)) { master$SMARCA4 <- 1 }
p1 <- ggplot(master, aes(factor(V2), log2(SMARCB1))) + geom_violin(trim = FALSE, aes(fill = factor(V2)), scale = "width") +
  geom_jitter(height = 0, width = 0.4, size = 0.001) + labs(x = "Clusters",y = "SMARCB1 CNV Score") + #ylim(0.75,1.25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 10, face = "bold", hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", hjust = 1, color = "black"))
p1
ggsave("/Volumes/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_ALL/by_Seurat/03_inferCNV.SMARCB1.pdf", p1, width = 20, height = 10, units = "in", useDingbats = F)

p2 <- ggplot(master, aes(factor(V2), log2(SMARCA4))) + geom_violin(trim = FALSE, aes(fill = factor(V2)), scale = "width") +
  geom_jitter(height = 0, width = 0.4, size = 0.001) + labs(x = "Clusters",y = "SMARCA4 CNV Score") + #ylim(0,3500) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.title = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"),
        axis.text.x = element_text(size = 10, face = "bold", hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 10, face = "bold", hjust = 1, color = "black"))
p2
ggsave("/Volumes/projects/brainTumor_JY/yu3grp/ATRT/qpan/Combined_ALL/by_Seurat/03_inferCNV.SMARCA4.pdf", p2, width = 20, height = 10, units = "in", useDingbats = F)
