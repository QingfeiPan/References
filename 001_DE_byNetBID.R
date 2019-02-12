## Load the required packages and functions
require(affy)
require(arm)
require(MCMCglmm)
require(plyr)
require(dplyr)
require(XLConnect)
require(stringr)

source("/Users/qpan/Projects/Training/Activity/include_R/de.inc.R")
source("/Users/qpan/Projects/Training/Activity/include_R/metaanalysis.inc.R")

load("/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/01_mouseModel/02_microArray.29mouseModel.8538_385.activity_scaled.gene.eset")
eset <- microArray.29mouseModel.8538_385.activity_scaled.gene.eset; rm(microArray.29mouseModel.8538_385.activity_scaled.gene.eset)

outdir <- "/Volumes/yu3grp/solidTumor_JY/yu3grp/BRCA/HDAC6/Qingfei/04_Datasets/01_mouseModel/04_DA"
mouseModels <- as.vector(sort(unique(pData(eset)$Class.Name)))

compare <- data.frame()
for (i in 1:length(mouseModels)) {
    compare[i,1] <- mouseModels[i]
    compare[i,2] <- str_c(mouseModels[-i], collapse = "|")
}
names(compare) <- c("Foreground","Background")

## DE analysis for multi-comparison of different conditions
output_merge <- fData(eset)
for (i in 1:nrow(compare)) {
    cat("Comparison", i, ':', compare[i,1], 'VS', compare[i,2], '\n')
    if (str_count(compare[i,1], "\\|") > 0) {c_fg <- unlist(strsplit(compare[i,1], "\\|"))} else {c_fg <- compare[i,1]}
    if (str_count(compare[i,2], "\\|") > 0) {c_bg <- unlist(strsplit(compare[i,2], "\\|"))} else {c_bg <- compare[i,2]}
    eset.sel<-eset[,pData(eset)$Class.Name%in%c(c_fg, c_bg)] ## Select data with specific conditions
    
    ## comp <- factor(as.numeric(gsub(paste('^',comps[i,2],'$',sep=''),0,gsub(paste('^',comps[i,1],'$',sep=''),1,(pData(eset.sel)$Conditions))))) # Grouping. Make sure 1 for foreground, and 0 for background.
    cluster_labs <- pData(eset.sel)$Class.Name
    group_fg <- which(is.element(cluster_labs, c_fg)); group_bg <- which(is.element(cluster_labs, c_bg));
    group <- c(); group[group_fg] <- 1; group[group_bg] <- 0; group <- factor(group)
    
    ## DE-analysis
    df <- data.frame(GeneSet=as.character(fData(eset.sel)$GeneSet),exprs(eset.sel),stringsAsFactors=FALSE) ## Add a new column to store geneSymbols
    de.cur <- ddply(df, .(GeneSet), 'combRowEvid.2grps', comp=group, family=gaussian, method='Bayesian', n.iter=5000, nitt=25000, burnin=5000, thin=1,
                    pooling=c('full'), logTransformed=TRUE, restand=FALSE, average.method=c('geometric')) # The data in 'd' should be normalized and LOG2-transformed.
    head(de.cur) # FC here is not log-transformed. The sign of it only means Up or Down. The expressioin values have been logE transformed.
    output <- de.cur[, c(1,8)]
    FDR <- p.adjust(de.cur$pval.full,'BH')
    de.stat <- data.frame(row.names = de.cur$GeneSet, pval = de.cur$pval.full, fdr = FDR)
    
    ## Calculate expression and FDR
    df.tmp <- df[,-1]
    AveExp_fg <- as.vector(apply(as.data.frame(df.tmp[, group == 1]), 1, mean, na.rm = T))
    AveExp_bg <- as.vector(apply(as.data.frame(df.tmp[, group == 0]), 1, mean, na.rm = T))
    AveExp_FC <- AveExp_fg - AveExp_bg
    de.exp <- data.frame(row.names = df$GeneSet, AveExp_fg = AveExp_fg, AveExp_bg = AveExp_bg, AveExp_FC = AveExp_FC)
    
    output_f <- merge(de.exp, de.stat, by = "row.names", all.x = T)
    names(output_f) <- c("GeneSet", paste0("log2Exp_", compare[i,1]), paste0("log2Exp_", compare[i,2]), paste0("Log2FC_", compare[i,1], "vs", compare[i,2]), paste0("Pval_", compare[i,1], "vs", compare[i,2]), paste0("FDR_", compare[i,1], "vs", compare[i,2]))
    output_f <- output_f[order(output_f[,4], decreasing = T),]
    
    fileName <- gsub(pattern = ", ", replacement = "_", compare[i,1])
    fileName <- gsub(pattern = "\\/", replacement = "_", fileName)
    output_file = paste0(outdir, "/byClass/DA_NetBID_classVSrest_", fileName, ".txt")
    write.table(output_f, file = output_file, quote = F, row.names = F, col.names = T, sep = "\t")
    
    output_merge <- merge(output_merge, output_f, by='GeneSet', all.x=TRUE)
}

output_merge_file = paste0(outdir, "/byClass/DA_NetBID_classVSrest_Merged.txt")
output_merge = output_merge[order(output_merge[,4], decreasing = T),]
write.table(output_merge, file = output_merge_file, quote = F, row.names = F, col.names = T, sep = "\t")
