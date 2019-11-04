library("cowplot")
library("Seurat")
library("scran")
library("rsvd")
library(ggplot2)
library(gridExtra)

setwd("~/work/scRNA")

out_objects_dir <- "./results/R.out/data/Robjects"
out_plot_dir <- "./results/R.out/ALRA"
if(!dir.exists(out_plot_dir)) dir.create(out_plot_dir)

set.seed(1234)

dataList <- readRDS(file=file.path(out_objects_dir, "ExpressionList_script_4.rds"))

dataList <- readRDS(file=file.path(out_objects_dir, "ExpressionList_QC.rds"))

m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
fD$keep[is.na(fD$keep)] <- FALSE
rm(dataList)

m <- m[fD$keep, pD$PassAll]   ## 11707 genes X 15847 cells
pD <- pD[pD$PassAll, ]
rownames(pD)  <- pD[, 1]
fD <- fD[fD$keep, ]

## pig gene names
pig_genes <- read.delim("All.pig.gene.plus.HGNC.names.txt", header = TRUE, as.is = TRUE)
pig_genes_HGNC <- pig_genes$Gene.name
names(pig_genes_HGNC) <- pig_genes$Gene.stable.ID
rownames(m) <- pig_genes_HGNC[rownames(m)]

## !!change sample_name, num_samples here to choose samples!!
## options: 2I, IPP, No-PP, PBMC
sample_name <- "PBMC"
num_samples <- 2


i <- which(pD$SampleID == paste0("intestine_trial_6_Pig-1-", sample_name))
j <- which(pD$SampleID == paste0("intestine_trial_6_Pig-2-", sample_name))
index2i <- c(i,j)
pD <- pD[index2i,] #use indices to keep only the rows that correspond to this sample
m <- m[,index2i] #use indices from pD to keep only the barcode cols that correspond
identical((length(i)+length(j)), ncol(m), nrow(pD)) #check that size is correct

samp <- CreateSeuratObject(counts = m, meta.data = pD)
samp.list <- SplitObject(object = samp, split.by = "SampleID")

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
samp.list <- lapply(samp.list, function(.x){
    temp <- NormalizeData(object = .x)
    temp <- FindVariableFeatures(object = temp, do.plot = FALSE)
    temp
})

# ALRA is undeterministic
# apply ALRA to impute non-technical zeroes
#  A: k=18,  7.49% nonzero to 47.66% nonzero
#  B: k=18,  6.97% nonzero to 46.35% nonzero
#  C: k=20,  7.10% nonzero to 43.39% nonzero
#  CT2-1NOV: k=23, 9.31% nonzero to 44.40% nonzero
#  CT2-30OCT: k=24, 9.48% nonzero to 46.32% nonzero

# Choose k.

choose_k <- function (A_norm,K=100, pval_thresh=1E-10, noise_start=80,q=2) {
    #  Heuristic for choosing rank k for the low rank approximation based on
    #  statistics of the spacings between consecutive singular values. Finds
    #  the smallest singular value \sigma_i such that $\sigma_i - \sigma_{i-1}
    #  is significantly different than spacings in the tail of the singular values.
    #
    #
    # Args:
    #   A_norm: The log-transformed expression matrix of cells (rows) vs. genes (columns)
    #   K: Number of singular values to compute. Must be less than the smallest dimension of the matrix.
    #   pval_thresh : The threshold for ``significance''
    #   noise_start : Index for which all smaller singular values are considered noise
    #   q : Number of additional power iterations
    #
    # Returns:
    #   A list with three items
    #       1) Chosen k
    #       2) P values of each possible k
    #       3) Singular values of the matrix A_norm

    if (K > min(dim(A_norm))) {
        stop("For an m by n matrix, K must be smaller than the min(m,n).\n")
    }
    if (noise_start > K-5) {
        stop("There need to be at least 5 singular values considered noise.\n")
    }
    noise_svals <- noise_start:K
    rsvd_out <- rsvd(A_norm,K,q=q)
    diffs <- diff(rsvd_out$d)
    pvals <- pnorm(diffs,mean(diffs[noise_svals-1]),sd(diffs[noise_svals-1]))
    k <- max(which( pvals  <pval_thresh))
    return (list( k=k, pvals =pvals, d=rsvd_out$d))
}

k_choice <- lapply(samp.list, function(.x) {
    choose_k(as.matrix(.x@assays$RNA@data))
})
# For the results in the paper, automatically chosen k worked quite well, but in
# some cases you might want to take a closer look, as we do here. The k is
# chosen based on the spacings between the singular values, as it can be quite
# hard to identify the ``beginning of noise'' from just looking at the spectrum
# itself. Uncomment the code below to plot them


for (i in 1:num_samples)
{
    pdf(file.path(out_plot_dir, paste0("k-value.chose.samp", i, ".", sample_name, ".pdf")), width = 25 , height = 8)
    df <- data.frame(x=1:100,y=k_choice[[i]]$d)
    g1<-ggplot(df,aes(x=x,y=y)) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice[[i]]$k)   + theme( axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_i') + ggtitle('Singular values')
    df <- data.frame(x=2:100,y=diff(k_choice[[i]]$d))[3:99,]
    g2<-ggplot(df,aes(x=x,y=y)) + geom_point(size=1)  + geom_line(size=0.5)+ geom_vline(xintercept=k_choice[[i]]$k+1)   + theme(axis.title.x=element_blank() ) + scale_x_continuous(breaks=seq(10,100,10)) + ylab('s_{i} - s_{i-1}') + ggtitle('Singular value spacings')
    grid.arrange(g1,g2,nrow=1)

    dev.off()
}

#visualize k_choice, if it looks reasonable, use these automatically generated k's
ks <- c()
for (i in 1:num_samples)
{
	ks <- c(ks, k_choice[[i]]$k)
}

## otherwise, set your own k's
## ks <- c(35,34)

## for No-PP samples
## ks <- c(30,33)

samp.list <- mapply(function(.x, .y){
    temp <- RunALRA(object=.x, k = .y, q = 10, assay = NULL,
            slot = "data", setDefaultAssay = TRUE, genes.use = NULL,
            K = NULL, p.val.th = 1e-10, noise.start = NULL, q.k = 2,
            k.only = FALSE)
    temp
}, samp.list, ks)

save(file=file.path(out_objects_dir, paste0(sample_name, ".ALRA.imputed.RData")))


###  Apply Seurat to imputed data

load(file=file.path(out_objects_dir, paste0(sample_name, ".ALRA.imputed.RData")))
 
# ###### reconstruct Seurat object using imputed expression
# 
counts <- do.call(cbind, lapply(samp.list, function(.x){
    as.matrix(.x@assays$alra@data)
}))
samp <- CreateSeuratObject(counts = counts, meta.data = pD)
samp.list <- SplitObject(object = samp, split.by = "SampleID")

samp.list <- lapply(samp.list, function(.x){
    temp <- FindVariableFeatures(object = .x, do.plot = FALSE)
    temp
})

samp_int <- FindIntegrationAnchors(object.list = samp.list, scale = TRUE, dims = 1:30)

## get integrated expression values for all genes
samp.integrated <- IntegrateData(anchorset = samp_int,
                                 features.to.integrate = rownames(m), dims = 1:30)


## integrated analysis
DefaultAssay(object = samp.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
regression <- TRUE

if (regression){
    samp.integrated <- ScaleData(object = samp.integrated,
                                 vars.to.regress = c("UmiSums", "prcntMito"), verbose = FALSE)
} else {
    ## no regression out
    samp.integrated <- ScaleData(object = samp.integrated, verbose = FALSE)
}

# # pdf(file.path(out_plot_dir, paste0(sample_name, " Highly variable gene plot.pdf")), width =15, height =5)
# # samp.integrated <- FindVariableFeatures(object = samp.integrated, do.plot = TRUE)
# # dev.off()
# 
samp.integrated <- RunPCA(object = samp.integrated, features = samp.integrated$integrated@var.features, npcs = 50, verbose = FALSE)

# ## plot variance
# sd <- samp.integrated@reductions$pca@stdev
# var <- sd^2/(sum(sd^2))*100
# 
# pdf(file.path(out_plot_dir, "vairance for PCA-ALRA.pdf"), height = 10, width = 10)
# plot(x=1:50, y=var, pch = 16, type= "b", ylab= "Variance (%)", xlab = "Principle component")
# dev.off()
# 
# ## UMAP: This depends on python package umap-learn
# 
samp.integrated <- RunUMAP(object = samp.integrated,  seed.use = 423,
                           reduction = "pca", dims = 1:22)
## TSNE
samp.integrated <- RunTSNE(object = samp.integrated, reduction = "pca", k.seed = 2,
                           dims = 1:22)
samp.integrated <- FindNeighbors(object = samp.integrated, reduction = "pca",
                                 dims = 1:22 )
samp.integrated <- FindClusters(object = samp.integrated, reduction = "pca",
                                dims = 1:22, save.SNN = TRUE)

								
## convert to Single Cell Experiment object and save rds
s.sce <-as.SingleCellExperiment(samp.integrated)
saveRDS(s.sce,file=file.path(out_objects_dir, paste0(sample_name,".imputed.sce.rds")))

## save Seurat object, too
saveRDS(samp.integrated,file=file.path(out_objects_dir, paste0(sample_name, ".imputed.seurat.rds")))


pdf(file.path(out_plot_dir, paste0(sample_name, " 18 PCA-Tsne and Umap plot of cell clusters-ALRA.pdf")), height = 12, width =15,)
p1 <- DimPlot(object = samp.integrated, reduction = "tsne", group.by = "SampleID", pt.size =0.5) + xlim(-60,60) + ylim(-60,60)
p2 <- DimPlot(object = samp.integrated, reduction = "tsne", do.return = TRUE, label = TRUE,  pt.size = 0.5) + xlim(-50,60) + ylim(-60,60)
p3 <- DimPlot(object = samp.integrated, reduction = "umap", group.by = "SampleID", pt.size =0.5) + xlim(-15,15) + ylim(-20,20)
p4 <- DimPlot(object = samp.integrated, reduction = "umap", do.return = TRUE, label = TRUE,  pt.size = 0.5) + xlim(-15,15) + ylim(-20,20)

plot_grid(p1, p2, p3, p4, nrow =2)
dev.off()


save.image(file=file.path(out_objects_dir, paste0(sample_name, ".Seurat.integrated.all.features-ALRA.RData")))

## create PDFs for each marker, tSNE and UMAP
overlay <- function(file_name, width, height, features, reduction, xlim, ylim)
{
    pdf(file.path(out_plot_dir,file_name), height = height, width = width)
    print(FeaturePlot(object = samp.integrated, features = features,
                      dims = c(1, 2), cells = NULL,
                      cols = c("lightgrey", "red"), pt.size = 1, min.cutoff = "q9",
                      max.cutoff = NA, reduction = reduction, split.by = NULL,
                      shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
                      order = NULL, label = FALSE, label.size = 4, ncol = NULL,
                      combine = TRUE, coord.fixed = FALSE, aspect) + xlim(-15,15) + ylim(-20,20))
    dev.off()
}
					 
markers1 <- c("CD3E", "CD4","CD5", "CD8A", "CD8B", "TRDC",  
             "GZMB", "IFNG", "CD79A", "CD79B", "CD19") 

markers2 <- c("CD69", "MS4A1", "FCER1G", "MS4A2", "JCHAIN",  
             "ITGAM","FCGR1A", "CD14", "SERPING1", "MX1", "EPCAM")
			 
markers3 <- c("IL1RAP", "IFNGR1", "CST3", "TLR4",
             "NCR1", "KLRB1",  "GNLY", "LYZ", "MCM2",
             "MCM3", "TOP2A", "CCNB1", "PCNA")
## print 1st
overlay(file_name = paste0("1.overlay of markers on tSNE-ALRA.", sample_name, ".pdf"),
        height =25, width = 25, reduction = "tsne",
        features = markers1, xlim = "-60,60", ylim = "-60,60")
overlay(file_name = paste0("1.overlay of markers on UMAP-ALRA.", sample_name, ".pdf"),
        height =25, width = 25, reduction = "umap",
        features = markers1, xlim = "-15,15", ylim = "20,20")
## print 2nd		
overlay(file_name = paste0("2.overlay of markers on tSNE-ALRA.", sample_name, ".pdf"),
        height =25, width = 25, reduction = "tsne",
        features = markers2, xlim = "-60,60", ylim = "-60,60")
overlay(file_name = paste0("2.overlay of markers on UMAP-ALRA.", sample_name, ".pdf"),
        height =25, width = 25, reduction = "umap",
        features = markers2, xlim = "-15,15", ylim = "-20,20")
## print 3rd	
overlay(file_name = paste0("3.overlay of markers on tSNE-ALRA.", sample_name, ".pdf"),
        height =25, width = 25, reduction = "tsne",
        features = markers3, xlim = "-60,60", ylim = "-60,60")
overlay(file_name = paste0("3.overlay of markers on UMAP-ALRA.", sample_name, ".pdf"),
        height =25, width = 25, reduction = "umap",
        features = markers3, xlim = "-15,15", ylim = "-20,20")

overlay.p2 <- FeaturePlot(object = samp.integrated, features = markers1,
                          dims = c(1, 2), cells = NULL,
                          cols = c("grey", "red"), 
                          pt.size = 0.3, min.cutoff = "q9",
                          max.cutoff = NA, reduction = "umap", split.by = NULL,
                          shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
                          order = NULL, label = FALSE, label.size = 4, ncol = 6,
                          combine = FALSE, coord.fixed = FALSE)
						  
p2 <- lapply(overlay.p2, function(.x){
    .x + xlim(-15,15) + ylim(-20,20)
})

pdf(file.path(out_plot_dir,paste0(sample_name, ".UMAP markers1.pdf")), width = 22, height = 25)
grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], p2[[5]], p2[[6]],
             p2[[7]],  p2[[8]], p2[[9]], p2[[10]],p2[[11]])
dev.off()
		
		
### plot all naostring detected genes
nanostring_genes <- read.delim("nanostring target ensembl gene ID.txt", header = TRUE, as.is = TRUE)
nanostring_genes <- nanostring_genes$Gene.name[nanostring_genes$Gene.stable.ID != ""]
nanostring_genes <- nanostring_genes[nanostring_genes %in% rownames(m)]

for (i in seq(1, length(nanostring_genes), 24))
{
    if( i <= 160)
    {
        overlay(file_name= paste0("Genes.expr.overlay.on.clusters-ALRA", i, ".pdf"),
                width = 25, height =25, reduction = "tsne",
               features = nanostring_genes[i:(i + 23)])
        next
    }else{
        overlay(file_name= paste0("Genes.expr.overlay.on.clusters-ALRA", i, ".pdf"),
                width = 25, height = 20, reduction = "tsne",
                features = nanostring_genes[i:length(nanostring_genes)])
    }
}

save.image(file=file.path(out_objects_dir, "ALRA.imputed.RData"))
##
all_markers <- FindAllMarkers(object = samp.integrated, test.use = "wilcox")
write.table(all_markers, file.path(out_plot_dir, "cluster-specific.markers.genes.across.all.clusters-ALRA.txt"), sep ="\t", quote = FALSE, row.names =TRUE)


#all_markers <- read.delim(file.path(out_plot_dir, "cluster-specific.markers.genes.across.all.clusters-ALRA.txt"), as.is = TRUE)
markers.use  <- subset(all_markers, avg_logFC >= 1 & p_val_adj <= 0.05)$gene

markers.use <- markers.use[markers.use %in% rownames(samp.integrated@assays$integrated@scale.data) ]

pdf(file.path(out_plot_dir,"clusterwise.markers.heatmap-ALRA.pdf"), height= 25, width = 15)
DoHeatmap(object = samp.integrated, features = markers.use, cells = NULL, 
          group.by = "ident", size = 2.5,
          group.bar = TRUE, disp.min = -3, disp.max = 3,
          slot = "scale.data", assay = "integrated", label = TRUE,
          hjust = 0, angle = 90, combine = TRUE)

dev.off()