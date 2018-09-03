library(utils)
library(ggplot2)


input <- "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/02_scMINER/Combined_Total/scMINER_Combined_Total/scMINER_Combined_Total_MDS_16/scMINER_MICA_out/Combined_Total.ggplot.txt"
output <- "/Volumes/yu3grp/scRNASeq/yu3grp/metastasis/02_Breast_cancer/02_scMINER/Combined_Total/scMINER_Combined_Total/scMINER_MICA_figures/Combined_Total_clust_k16_mod.rplot.pdf"

cc = t(read.table(input, sep="\t"))
ncols = as.numeric(ncol(cc))
cells = t(cc[1, 2:ncols])

x_tsne <- as.numeric(t(cc[2, 2:ncols]))
y_tsne <- as.numeric(t(cc[3, 2:ncols]))
c <- as.numeric(t(cc[4, 2:ncols])) + 1
uc <- unique(c)
n = as.numeric(length(uc))
df <- data.frame(x_tsne, y_tsne, c)
x_pose = rep(c(0), n)
y_pose = rep(c(0), n)
label = c(1:n)
l_size = rep(c(0), n)
for(i in 1:n){
    x_pose[i] <- sum(x_tsne[c==i]) / length(c[c==i]);
    y_pose[i] <- sum(y_tsne[c==i]) / length(c[c==i]);
    l_size[i] <- paste(c(i, "(", length(x_tsne[c==i]), ")"), collapse = '')
    label[i] <- i;
}
df1 <- data.frame(x_pose, y_pose, label)

## 1 Specific colors automatically
p <- ggplot() +
    geom_point(data=df, aes(x=x_tsne, y=y_tsne, color=as.factor(c)), size=0.3) +
    xlab("MICA-1") + ylab("MICA-2") +
    guides(col=guide_legend(override.aes = list(size=10), title=paste(c("Clusters\n(", toString(ncols-1), ")"), collapse = '')[1]), ncol=ceiling(n/10)) +
    geom_text(data=df1, aes(x=x_pose, y=y_pose, label=NA), size=9) +
    scale_color_discrete(labels=l_size) +
    ggtitle(paste0("Combined_Total", " (", toString(n), ")", collapse = '')[1]) +
    theme(plot.title = element_text(size=25, margin=margin(t=20, b=10)),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.text = element_text(size=12),
          legend.title = element_text(size=14,face="bold"))
ggsave(output, plot = p, width = 12, height = 12)

## Specific colors manually
p <- ggplot() +
    geom_point(data=df[df$c %in% c(1, 3), ], aes(x=x_tsne, y=y_tsne), color="darkorchid1", size=0.5) +
    geom_point(data=df[df$c %in% c(2), ], aes(x=x_tsne, y=y_tsne), color="magenta", size=0.5) +
    geom_point(data=df[df$c %in% c(4), ], aes(x=x_tsne, y=y_tsne), color="indianred1", size=0.5) +
    geom_point(data=df[df$c %in% c(5), ], aes(x=x_tsne, y=y_tsne), color="red", size=0.5) +
    geom_point(data=df[df$c %in% c(6), ], aes(x=x_tsne, y=y_tsne), color="chartreuse3", size=0.5) +
    geom_point(data=df[df$c %in% c(7), ], aes(x=x_tsne, y=y_tsne), color="royalblue", size=0.5) +
    geom_point(data=df[df$c %in% c(8), ], aes(x=x_tsne, y=y_tsne), color="mediumseagreen", size=0.5) +
    geom_point(data=df[df$y_tsne>120, ], aes(x=x_tsne, y=y_tsne), color="turquoise3", size=0.5) +
    geom_point(data=df[(df$x_tsne>0 & df$y_tsne < -130), ], aes(x=x_tsne, y=y_tsne), color="orange", size=0.5) +
    xlab("MICA-1") + ylab("MICA-2") +
    ggtitle(paste0("Combined_Total", " (", toString(n), ")", collapse = '')[1]) +
    theme(plot.title = element_text(size=25, margin=margin(t=20, b=10)),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.position = "NA")
ggsave(output, plot = p, width = 12, height = 12)
