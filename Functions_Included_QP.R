
#################### Matrix Normalizatioin ####################
# 1) d: the matrix for normalization, the header and row.names are assigned.
# 2) total: the library size specified.
# 3) pseudoCount: additive smoothing. '0' only means not observed, but not an impossible event.
normalize.scale <- function(d, total=NULL, pseudoCount=1) {
    if (!is.data.frame(d)) d <- data.frame(d) # Must be data frame
    if (!all(d>0)) d <- d+pseudoCount # 'pseudoCount' will be used if there are '0' in the matrix.
    s <- apply(d, 2, sum) # sum of columns
    m <- ifelse(is.null(total), as.integer(mean(s)), as.integer(total)) # library size: mean or total
    options(digits = 2+nchar(m)) # control the number of significant digits to print
    fac <- m/s
    for (i in 1:length(s)) {
        d[,i] <- round(d[,i]*fac[i], 0) # No decimal places.
    }
    if (!all(d>0)) d <- d+pseudoCount
    d
}


#################### Quality Control Report ####################
# Find "arrayQualityMetrics" package for details.
# 1) eset: ExpressionSet created by Biobase package.
# 2) outputdir: Output directory. This directory will be created unless it exists.
# 3) groupsize: the number of conditions will be shown in one QC report
# 4) intgroup: name of the column with sample grouping information.
# 5) do.logtransform: do logrithm transformation or not.
# 6) do.all: ???
# 7) grouprep: if there are replicates in each group. Use "TRUE" always.
# 8) force: only works for the situation in which the outputdir already exists. If TRUE, the outputdir will be overwritten, otherwise an error is thrown.
getIndex <- function(query, vector, first=FALSE) {
    index <- which(vector==query)
    if(first) index[1]
    else index
}

QC.expresso <- function(eset, outputdir="output", groupsize=99, intgroup='Group', do.logtransform=TRUE, do.all=TRUE, grouprep=TRUE, force=TRUE, ...) {
    group <- unique(pData(eset)[, intgroup]) # Extract the grouping types
    if (length(group) > groupsize) { # if the number of conditions is more than 'groupsize', subfolders will be made.
        step <- groupsize-1
        i <- 1
        while (i>0 && i<=(length(group))) {
            idx <- unlist(sapply(group[i:(i+step)], getIndex, pData(eset)[, intgroup]))
            if (length(idx)>1) {
                if (i <= length(group)-step) {outdir = file.path(outputdir, paste('group_', i, '_', i+step, sep=""))}
                else {outdir = file.path(outputdir, paste('group_', i, '_', length(group), sep=""))}
                arrayQualityMetrics(expressionset=eset[,idx], outdir=outdir, force=force, do.logtransform=do.logtransform, intgroup=intgroup, grouprep=grouprep, ...)
            }
            i <- i+step+1
        }
    }
    
    if (length(group)<=groupsize || (length(group)>groupsize && do.all)) { # Once 'do.all' is true, all conditions will be shown in one QC report.
        arrayQualityMetrics(expressionset=eset, outdir=file.path(outputdir), force=force, do.logtransform=do.logtransform, intgroup=intgroup, grouprep=grouprep, ...)
    }
}
