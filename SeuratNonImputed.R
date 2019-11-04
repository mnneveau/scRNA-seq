library(cowplot)
library(Seurat)
setwd("~/work/scRNA")
out_objects_dir <- "./results/R.out/data/Robjects"
out_plot_dir <- "./results/R.out"

dataList <- readRDS(file=file.path(out_objects_dir, "ExpressionList_QC.rds"))
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
fD$keep[is.na(fD$keep)] <- FALSE
rm(dataList)

# Gene and cell filtering

m <- m[fD$keep, pD$PassAll]   ## 11707 genes X 15847 cells
pD <- pD[pD$PassAll, ]
rownames(pD)  <- pD[, 1]
fD <- fD[fD$keep, ]

#begins unique sample code
#"-IPP" can also be "-No-PP", "PBMC", "2I"
i <- which(pD$SampleID == "intestine_trial_6_Pig-1-IPP")
j <- which(pD$SampleID == "intestine_trial_6_Pig-2-IPP")
indexIPP <- c(i,j)
pDIPP <- pD[indexIPP,] #use indices to keep only the rows that correspond to this sample
mIPP <- m[,indexIPP] #use indices from pD to keep only the barcode cols that correspond
identical((length(i)+length(j)), ncol(mIPP), nrow(pDIPP)) #check that size is correct

sipp <- CreateSeuratObject(counts = mIPP, meta.data = pDIPP)
sipp.list <- SplitObject(object = sipp, split.by = "SampleID")

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here
sipp.list <- lapply(sipp.list, function(.x){
    temp <- NormalizeData(object = .x)
    temp <- FindVariableFeatures(object = temp, do.plot = FALSE)
    temp
})

### Integration of 3 pancreatic islet cell datasets
sipp_int <- FindIntegrationAnchors(object.list = sipp.list, dims = 1:30)
sipp.integrated <- IntegrateData(anchorset = sipp_int, dims = 1:30)


DefaultAssay(object = sipp.integrated) <- "integrated"

sipp.integrated <- ScaleData(object = sipp.integrated, vars.to.regress = c("UmiSums", "prcntMito"), verbose = FALSE)

sipp.integrated <- RunPCA(object = sipp.integrated, features = sipp.integrated$integrated@var.features, npcs = 50, verbose = FALSE)

## plot variance
sd <- sipp.integrated@reductions$pca@stdev
var <- sd^2/(sum(sd^2))*100

pdf(file.path(out_plot_dir, "variance for PCA-IPP.pdf"), height = 10, width = 10)
plot(x=1:50, y=var, pch = 16, type= "b", ylab= "Variance (%)", xlab = "Principle component")
dev.off()

## UMAP: This depends on python package umap-learn

sipp.integrated <- RunUMAP(object = sipp.integrated, reduction = "pca", dims = 1:18)
						   
## TSNE
sipp.integrated <- RunTSNE(object = sipp.integrated, reduction = "pca", 
                           dims = 1:18)
sipp.integrated <- FindNeighbors(object = sipp.integrated, reduction = "pca", 
              dims = 1:18 )
sipp.integrated <- FindClusters(object = sipp.integrated, reduction = "pca", 
                                dims = 1:18, save.SNN = TRUE)
								
pdf(file.path(out_plot_dir, "18 PCA-Tsne and Umap plot of cell clusters-IPP.pdf"), width =15, height =12)
p1 <- DimPlot(object = sipp.integrated, reduction = "tsne", group.by = "SampleID", pt.size =0.5)
p2 <- DimPlot(object = sipp.integrated, reduction = "tsne", do.return = TRUE, label = TRUE,  pt.size = 0.5)
p3 <- DimPlot(object = sipp.integrated, reduction = "umap", group.by = "SampleID", pt.size =0.5)
p4 <- DimPlot(object = sipp.integrated, reduction = "umap", do.return = TRUE, label = TRUE,  pt.size = 0.5)

plot_grid(p1, p2, p3, p4, nrow =2)
dev.off()

##takes long time, ~25 min for 2I 
all_markers <- FindAllMarkers(object = sipp.integrated, test.use = "wilcox")

markers.use=subset(all_markers, avg_logFC >= 1)$gene

pdf("IPP.Markers.plot.pdf", height= 20, width =15)
DoHeatmap(object = sipp.integrated, features = markers.use, cells = NULL, 
          group.by = "ident", size =1.5,
          group.bar = TRUE, disp.min = -2.5, disp.max = NULL,
          slot = "scale.data", assay = NULL, label = TRUE,
          hjust = 0, angle = 90, combine = TRUE)

dev.off()

pdf(file.path(out_plot_dir,"IPP overlay of markers.pdf"), height =15, width = 25)
FeaturePlot(object = sipp.integrated, features = c("ENSSSCG00000008228", ## GNLY
                         "ENSSSCG00000013115", ## CD5
                         "ENSSSCG00000017283", ## CD79B
                         "ENSSSCG00000021812", ## MS4A1
                         "ENSSSCG00000000653", ## CD69
                         "ENSSSCG00000000492", ## LYZ
                         "ENSSSCG00000040140", ## CD3E
                         "ENSSSCG00000008217", ## CD8A
                         "ENSSSCG00000033684", ## CD79A
                         "ENSSSCG00000022512", ## TRDC 
                         "ENSSSCG00000034506", ## MS4A2
                         "ENSSSCG00000035379", ## JCHAIN
                         "ENSSSCG00000022675", ## NCR1
                         "ENSSSCG00000006357", ## FCER1G
                         "ENSSSCG00000037360", ## CST3
                         "ENSSSCG00000007978"), dims = c(1, 2), cells = NULL,
            cols = c("lightgrey", "red"), pt.size = 1, min.cutoff = "q9",
            max.cutoff = NA, reduction = "tsne", split.by = NULL,
            shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
            order = NULL, label = TRUE, label.size = 4, ncol = NULL,
            combine = TRUE, coord.fixed = FALSE)
dev.off()

##Warning: Could not find ENSSSCG00000017283 CD79B in the default search locations, found in RNA assay instead
##Warning: Could not find ENSSSCG00000033684 CD79A in the default search locations, found in RNA assay instead


pdf(file.path(out_plot_dir,"IPP overlay of markers on umap.pdf"), height =15, width = 25)
FeaturePlot(object = sipp.integrated, features = c("ENSSSCG00000008228", ## GNLY
                         "ENSSSCG00000013115", ## CD5
                         "ENSSSCG00000017283", ## CD79B
                         "ENSSSCG00000021812", ## MS4A1
                         "ENSSSCG00000000653", ## CD69
                         "ENSSSCG00000000492", ## LYZ
                         "ENSSSCG00000040140", ## CD3E
                         "ENSSSCG00000008217", ## CD8A
                         "ENSSSCG00000033684", ## CD79A
                         "ENSSSCG00000022512", ## TRDC 
                         "ENSSSCG00000034506", ## MS4A2
                         "ENSSSCG00000035379", ## JCHAIN
                         "ENSSSCG00000022675", ## NCR1
                         "ENSSSCG00000006357", ## FCER1G
                         "ENSSSCG00000013115", ## CD5
                         "ENSSSCG00000037360", ## CST3
                         "ENSSSCG00000007978"), dims = c(1, 2), cells = NULL,
            cols = c("lightgrey", "red"), pt.size = 1, min.cutoff = "q9",
            max.cutoff = NA, reduction = "umap", split.by = NULL,
            shape.by = NULL, blend = FALSE, blend.threshold = 0.5,
            order = NULL, label = TRUE, label.size = 4, ncol = NULL,
            combine = TRUE, coord.fixed = FALSE)
dev.off()

save.image(file = "Seurat.RData")


