#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
#install.packages("jaccard")

library(cowplot)
library(Seurat)
library(gplots)
library(jaccard)
library(RColorBrewer)
display.brewer.all()

out_objects_dir <- "./results/R.out/data/Robjects"
out_plot_dir <- "./results/R.out/plots"
setwd("~/work/scRNA")

sce.2I <- readRDS(file=file.path(out_objects_dir, "PBMC.sce.rds"))
seurat.2I <- as.Seurat(sce.2I, counts="logcounts")

sce.2I.imp <- readRDS(file=file.path(out_objects_dir, 
                                     "PBMC.imputed.sce.rds"))
seurat.2I.imp <- as.Seurat(sce.2I.imp, counts="logcounts")

#get cluster information from seurat objects
df.all <- as.data.frame(x=as.matrix(x= seurat.2I$seurat_clusters))
df.all.imp <- as.data.frame(x=as.matrix(x= seurat.2I.imp$seurat_clusters))
names(df.all)[1] <- "Before_Cluster"
names(df.all.imp)[1] <- "After_Cluster"
#merged before and after imputation data tables
df.all["After_Cluster"] <- df.all.imp["After_Cluster"]
#convert cols to numerical values for easy comparison
df.all$Before_Cluster <- as.numeric(as.character(df.all$Before_Cluster))
df.all$After_Cluster <- as.numeric(as.character(df.all$After_Cluster))

#create list of vectors of all barcodes contained by each cluster, before and after
b.v.list <- list()
a.v.list <- list()
for (i in 0:max(df.all$After_Cluster))
{
  a <- df.all[which(df.all$After_Cluster == i),]
  a.v <- rownames(a)
  a.v.list[[i+1]] <- a.v
}
for (i in 0:max(df.all$Before_Cluster))
{
  b <- df.all[ which(df.all$Before_Cluster == i),]
  b.v <- rownames(b)
  b.v.list[[i+1]] <- b.v
}

#make empty matrix
nclust.after <- seq(0,max(df.all$After_Cluster))
nclust.before <- seq(0,max(df.all$Before_Cluster))
jaccard.m <- matrix(nrow=(max(df.all$Before_Cluster)+1),
                    ncol=(max(df.all$After_Cluster)+1),
                    dimnames = list(nclust.before, nclust.after))

#compute jaccard index, store in matrix
for (i in 0:max(df.all$Before_Cluster))
{
  for (j in 0:max(df.all$After_Cluster))
  {
    inter <- intersect(b.v.list[[i+1]],a.v.list[[j+1]])
    u <- union(b.v.list[[i+1]],a.v.list[[j+1]])
    u.l <- length(u)
    inter.l <- length(inter)
    jaccard.val <- inter.l / u.l
    jaccard.m[i+1,j+1] <- jaccard.val
  }
}

#output heatmap for jaccard index
pdf("PBMC.jaccard.heatmap.pdf")
map <- heatmap.2(jaccard.m, Colv = NA, Rowv = NA, xlab = "After Imputation",
                 ylab = "Before Imputation", key = TRUE, cexRow = 0.7, cexCol = 1.2, 
                 tracecol = NA, col = brewer.pal(n=8, name = "Reds"), density.info = 'none',
                 key.xlab = "Jaccard Index",
                 main = "PBMC Cluster Identity by Barcode")
dev.off()

#compute absolute value of cell's cluster identity before and after imputation
nclust.after <- seq(0,max(df.all$After_Cluster))
nclust.before <- seq(0,max(df.all$Before_Cluster))
m <- matrix(0L,nrow=(max(df.all$Before_Cluster)+1),
                    ncol=(max(df.all$After_Cluster)+1),
                    dimnames = list(nclust.before, nclust.after))

for (row in 1:nrow(df.all))
{
  m[(df.all$Before_Cluster[row]+1),(df.all$After_Cluster[row]+1)] = (
    m[(df.all$Before_Cluster[row]+1),(df.all$After_Cluster[row]+1)] + 1)
}


identical(as.numeric(sum(m)), as.numeric(nrow(df.all)))

#create heatmap for absolute value with heatmap function
color = rev(heat.colors(256))
pdf("2Iheatmap.pdf")
map <- heatmap(m, Colv = NA, Rowv = NA, xlab= "After Imputation",
               ylab = "Before Imputation", col = color)
dev.off()

#create heatmap for absolute value with heatmap.2 function
color = rev(heat.colors(256))
pdf("2Iheatmap.2.pdf")
map <- heatmap.2(m, Colv = NA, Rowv = NA, xlab = "After Imputation",
                 ylab = "Before Imputation", key = FALSE, cellnote = m,
                 notecol = "black", cexRow = 0.7, cexCol = 1.2, 
                 tracecol = NA, keysize = 0.5, col = color)
dev.off()