library("ggplot2")
library("dplyr")
library("monocle")
library("stringr")
library("pracma")

## Input 
cluster <- read.table("/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/02_scMINER/Combined_Total/scMINER_Combined_Total/scMINER_Combined_Total_MDS_16/scMINER_MICA_out/Combined_Total.ggplot.txt", header = T, sep = "\t")
expression <- read.table("/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/02_scMINER/Combined_Total/3_1_Combined_Total_mm10_scMINER_total.txt", header = T, sep = "\t", row.names = 1)
outdir <- "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/03_cluster_genes/01_cluster_vs_rest/Combined_Total/small_islands"

## Prepare the expression matrix, normalized and non-log transformed
d <- t(expression); d <- data.frame(d); d <- exp(d) - 1

## Plot the cells to get the x- and y- coordinate boundaries
ggplot(cluster, aes(x = X, y = Y)) + geom_point(size=0.5) + xlim(-150, 150) + ylim(-150, 150) + geom_vline(xintercept = c(-52, -38)) + geom_hline(yintercept = c(35, 0))
ggplot(cluster[cluster$label==5,], aes(x = X, y = Y)) + geom_point(size=0.5) + vline(x=-35)

## Define the small clusters with x- and y- coordinate combinations
subcluster_7 <- filter(cluster, (X > -15 & X < 30) & (Y > -10 & Y < 40));
subcluster_center <- filter(cluster, (X > -38 & X < -20) & (Y > -5 & Y < 15))
subcluster_1_5_main <- filter(cluster, (X > -40 & X < 25) & (Y > -100 & Y < -10))
subcluster_1_5_side <- filter(cluster, (X > -75 & X < -40) & (Y > -80 & Y < -40))
subcluster_2_4_main <- filter(cluster, ((X > -100 & X < -52) & (Y > -25 & Y < 55) | (X > -52 & X < -38) & (Y > 0 & Y < 35)))
subcluster_2_4_side <- filter(cluster, (X > -52 & X < -30) & (Y > 35 & Y < 55))

## build the phenotype file
subcluster_7_pd <- data.frame(CellID = subcluster_7$GeneSymbol, Group = "subcluster_7")
subcluster_center_pd <- data.frame(CellID = subcluster_center$GeneSymbol, Group = "subcluster_center")
subcluster_1_5_main_pd <- data.frame(CellID = subcluster_1_5_main$GeneSymbol, Group = "subcluster_1_5_main")
subcluster_1_5_side_pd <- data.frame(CellID = subcluster_1_5_side$GeneSymbol, Group = "subcluster_1_5_side")
subcluster_2_4_main_pd <- data.frame(CellID = subcluster_2_4_main$GeneSymbol, Group = "subcluster_2_4_main")
subcluster_2_4_side_pd <- data.frame(CellID = subcluster_2_4_side$GeneSymbol, Group = "subcluster_2_4_side")
pd <- rbind(subcluster_7_pd, subcluster_center_pd, subcluster_1_5_main_pd, subcluster_1_5_side_pd, subcluster_2_4_main_pd, subcluster_2_4_side_pd)

## Check the gene coverage of each subcluster, most islands were caused by the batch effect.
subcluster_7_d <- select(d, subcluster_7_pd$CellID); subcluster_7_geneCoverage <- colSums(subcluster_7_d != 0)
subcluster_center_d <- select(d, subcluster_center_pd$CellID); subcluster_center_geneCoverage <- colSums(subcluster_center_d != 0)
subcluster_1_5_main_d <- select(d, subcluster_1_5_main_pd$CellID); subcluster_1_5_main_geneCoverage <- colSums(subcluster_1_5_main_d != 0)
subcluster_1_5_side_d <- select(d, subcluster_1_5_side_pd$CellID); subcluster_1_5_side_geneCoverage <- colSums(subcluster_1_5_side_d != 0)
subcluster_2_4_main_d <- select(d, subcluster_2_4_main_pd$CellID); subcluster_2_4_main_geneCoverage <- colSums(subcluster_2_4_main_d != 0)
subcluster_2_4_side_d <- select(d, subcluster_2_4_side_pd$CellID); subcluster_2_4_side_geneCoverage <- colSums(subcluster_2_4_side_d != 0)

Sample_ID <- c(rep("subcluster_7", length(subcluster_7_geneCoverage)),
               rep("subcluster_center", length(subcluster_center_geneCoverage)),
               rep("subcluster_1_5_main", length(subcluster_1_5_main_geneCoverage)),
               rep("subcluster_1_5_side", length(subcluster_1_5_side_geneCoverage)),
               rep("subcluster_2_4_main", length(subcluster_2_4_main_geneCoverage)),
               rep("subcluster_2_4_side", length(subcluster_2_4_side_geneCoverage)))
geneCoverage <- c(subcluster_7_geneCoverage, subcluster_center_geneCoverage, subcluster_1_5_main_geneCoverage, subcluster_1_5_side_geneCoverage, subcluster_2_4_main_geneCoverage, subcluster_2_4_side_geneCoverage)
boxplot_df <- data.frame(Sample_ID = Sample_ID, geneCoverage = geneCoverage)
p <- ggplot(boxplot_df, aes(x = Sample_ID, y = geneCoverage, fill = Sample_ID)) + geom_boxplot()
p <- ggplot(boxplot_df, aes(x = Sample_ID, y = geneCoverage, fill = Sample_ID)) + geom_violin() + geom_jitter(height = 0, width = 0.2, size = 0.001)
p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    labs(title = "Gene Coverage", x = "Sample", y = "Number of Identified Genes") +
    theme(legend.position = "none", plot.title = element_text(size = 20, face = "bold", hjust = 0.5, color = "black"), 
          axis.title = element_text(size = 15, face = "bold", hjust = 0.5, color = "black"),
          axis.text.x = element_text(size = 12, face = "bold", hjust = 1, color = "black", angle = 45),
          axis.text.y = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"))

## Prepare the compare file for DE analysis

compare <- matrix(c("subcluster_7", "subcluster_1_5_main",
                    "subcluster_center", "subcluster_1_5_main",
                    "subcluster_center", "subcluster_2_4_main",
                    "subcluster_1_5_side", "subcluster_1_5_main",
                    "subcluster_2_4_side", "subcluster_2_4_main"
                    ), ncol = 2, byrow = T)
colnames(compare) <- c("Forground", "Background")
compare <- data.frame(compare)

output_merge <- data.frame(GeneSymbol = row.names(d), row.names = row.names(d))
for (i in 1:nrow(compare)) {
    cat("Comparison", i, ':', as.character(compare[i, 1]), 'VS', as.character(compare[i, 2]), '\n')
    
    ## Expression table of CellDataSet Object
    if (str_count(as.character(compare[i,1]), ",") > 0) {c_fg <- unlist(strsplit(as.character(compare[i,1]), ","))} else {c_fg <- as.character(compare[i,1])}
    if (str_count(as.character(compare[i,2]), ",") > 0) {c_bg <- unlist(strsplit(as.character(compare[i,2]), ","))} else {c_bg <- as.character(compare[i,2])};
    d_sel <- select(d, pd$CellID[pd$Group%in%c(c_fg, c_bg)]) ## Select data with specific conditions
    
    ## Phenotype table of CellDataSet Object
    pd_sel <- pd[pd$Group%in%c(c_fg, c_bg),]
    group_fg <- which(is.element(pd_sel$Group, c_fg)); group_bg <- which(is.element(pd_sel$Group, c_bg));
    group <- c(); group[group_fg] <- 1; group[group_bg] <- 0; group <- factor(group); names(group) <- colnames(d_sel)
    pd_sel <- data.frame(Group=group); rownames(pd_sel) <- colnames(d_sel)
    
    ## FeatureInfo table of CellDataSet Object
    fd_sel <- data.frame(gene_short_name = row.names(d_sel), row.names = row.names(d_sel)) # To avoid the warnings, since the 'gene_short_name' column is required
    
    ## Create the CellDataSet Object
    pd_sel_obj <- new("AnnotatedDataFrame", data = pd_sel); fd_sel_obj <- new("AnnotatedDataFrame", data = fd_sel);
    DE_Obj <- newCellDataSet(as.matrix(d_sel), phenoData=pd_sel_obj, featureData = fd_sel_obj, expressionFamily=negbinomial.size())
    DE_Obj <- estimateSizeFactors(DE_Obj) # Calculate the size factors, which are scaling factors use as 'offsets' by the statistical model to make the different samples comparable.
    DE_Obj <- estimateDispersions(DE_Obj) # Calculate the dispersion. Some outlies may be found.
    
    ## DE analysis - calculate the Pval and Qval
    DE_res <- differentialGeneTest(DE_Obj, fullModelFormulaStr="~Group")
    DE_res <- DE_res[rownames(d_sel),];
    
    ## DE analysis - calculate the gene expression values with lsqnonneg()
    d_sel_log <- log(d_sel + 1)
    Ngenes <- dim(d_sel_log)[1]; Cells <- colnames(d_sel_log); Cell_levels <- levels(group[Cells]);
    Cell_groups <- match(group[Cells], Cell_levels);
    k <- max(Cell_groups); I <- 1:length(Cell_groups); W <- as.matrix(sparseMatrix(i=I, j=Cell_groups, x=rep.int(1,length(Cell_groups))));
    Exp_matrix <- matrix(0, k, Ngenes);
    for (h in 1:Ngenes){
        exp_gene <- as.vector(t(d_sel_log[h,])); Res_fit <- lsqnonneg(W, exp_gene); Exp_matrix[,h] <- Res_fit$x; # The Expression values are log-transformed.
    }
    colnames(Exp_matrix) <- rownames(d_sel_log); rownames(Exp_matrix) <- Cell_levels; H_0vs1 <- t(Exp_matrix); H_1vs0 <- H_0vs1[, c(2, 1)]
    
    ## DE analysis - Calculate foldchange
    H_1vs0 <- cbind(H_1vs0, FC = 0);
    for(j in 1:nrow(H_1vs0)) {
        H_1vs0[j,3] <- exp(H_1vs0[j,1] - H_1vs0[j,2])
    }
    
    ## Combine info together
    monocle_res <- merge(H_1vs0, DE_res[, c("pval", "qval")], by = "row.names", all = T); row.names(monocle_res) <- monocle_res$Row.names;
    monocle_res[,2] <- monocle_res[,2]*log(2.718281828)/log(2); monocle_res[,3] <- monocle_res[,3]*log(2.718281828)/log(2); ## LOG2 Transformation of Expression Data
    
    ## Add names of output matrix
    monocle_res <- data.frame(monocle_res);
    names(monocle_res) <- c("GeneSymbol", paste0("Log2Expr_", compare[i,1]), paste0("Log2Expr_", compare[i,2]), paste0("FoldChange_", compare[i,1], "vs", compare[i,2]), paste0("Pval_", compare[i,1], "vs", compare[i,2]), paste0("FDR_", compare[i,1], "vs", compare[i,2]))
    monocle_res <- monocle_res[order(monocle_res[,4], decreasing = T),];
    
    ## Write the output
    output_file = paste0(outdir, "/01_small_islands_", compare[i, 1], "_VS_", compare[i, 2], ".txt")
    write.table(monocle_res, file = output_file, quote = F, row.names = F, col.names = T, sep = "\t")
    
    ## Merge the result into matrix
    output_merge <- merge(output_merge, monocle_res, by='GeneSymbol', all.x=TRUE)
}
output_merge_file = paste0(outdir, "/01_small_islands__cellType_gene_Merged.txt")
output_merge = output_merge[order(output_merge[,4], decreasing = T),]
write.table(output_merge, file = output_merge_file, quote = F, row.names = F, col.names = T, sep = "\t")
