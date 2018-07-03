
## Input files
input_expression <- "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/01_Liver_cancer/03_cluster_genes/02_celltypes/01_Tumor_Cells/00_tumorCells_expressionMatrix.txt"
input_phenotype <- "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/01_Liver_cancer/03_cluster_genes/02_celltypes/01_Tumor_Cells/01_tumorCells_Loci_phenoData.txt"
input_comparison <- "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/01_Liver_cancer/03_cluster_genes/02_celltypes/01_Tumor_Cells/01_tumorCells_Loci_comps.txt"

## Output directory
output_dir <- "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/01_Liver_cancer/03_cluster_genes/02_celltypes/01_Tumor_Cells/01_tumorCells_Loci/new"

## Load the required packages and functions
library(monocle)
library(pracma)
library(stringr)

## Input 1: ExpressioinData
# Format: 1) Genes by Samples; 2) Counts for Monocle, and Log2Counts for NetBID.
d <- read.table(input_expression, header = T, sep = "\t", stringsAsFactors = F, row.names = 1) # Table of expression (Nomalized/Scaled & Log-transformed)
d <- data.frame(d);
d <- exp(d) - 1 # De-log-transformed to get the Nomalized/Scaled raw counts

## Input 2: PhenotypeData
# Format: 
pd <- read.table(input_phenotype, header = T, sep = "\t", stringsAsFactors = F, row.names = 1) # Table of sample/cell phenotypes
pd <- data.frame(pd); pd$CellID <- row.names(pd)
table(row.names(pd)==colnames(d)) # Check if all cellls involved in pd could be found in d
table(pd$Group) # The number of cells in each cluster

# To compare one cluster VS other specific cluster(s) by pre-prepared compare file
compare <- read.table(input_comparison, header = T, sep = "\t"); 
compare <- data.frame(compare)

# To compare one cluster VS other specific cluster(s) by pre-prepared compare file
## compare <- read.table(input_compare, header = T, sep = "\t"); compare <- data.frame(compare)

output_merge <- data.frame(GeneSymbol = row.names(d), row.names = row.names(d))
for (i in 1:nrow(compare)) {
  cat("Comparison", i, ':', as.character(compare[i, 1]), 'VS', as.character(compare[i, 2]), '\n')

  ## Expression table of CellDataSet Object
  if (str_count(as.character(compare[i,1]), ",") > 0) {c_fg <- unlist(strsplit(as.character(compare[i,1]), ","))} else {c_fg <- as.character(compare[i,1])}
  if (str_count(as.character(compare[i,2]), ",") > 0) {c_bg <- unlist(strsplit(as.character(compare[i,2]), ","))} else {c_bg <- as.character(compare[i,2])};
  d_sel<-d[,pd$Group%in%c(c_fg, c_bg)] ## Select data with specific conditions

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
  output_file = paste0(output_dir, "/01_tumorCells_Loci_cellType_gene_", compare[i, 1], "_VS_", compare[i, 2], ".txt")
  write.table(monocle_res, file = output_file, quote = F, row.names = F, col.names = T, sep = "\t")

  ## Merge the result into matrix
  output_merge <- merge(output_merge, monocle_res, by='GeneSymbol', all.x=TRUE)
}

output_merge_file = paste0(output_dir, "/01_tumorCells_Loci_cellType_gene__Merged.txt")
output_merge = output_merge[order(output_merge[,4], decreasing = T),]
write.table(output_merge, file = output_merge_file, quote = F, row.names = F, col.names = T, sep = "\t")
