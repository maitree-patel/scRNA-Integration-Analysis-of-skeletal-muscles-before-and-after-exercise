library(Seurat)

#data: scRNA-seq data analysis to compare the skeletal muscle transcriptome of an individual pre and post exercise
pre.exercise <- Read10X("/Users/maitreepatel/Desktop/R/scRNA analysis/GSE214544_GSM6611297 (2)")
pre.seurat <- CreateSeuratObject(counts = pre.exercise,
                                 project = "muscle_pre")

post.exercise <- Read10X("/Users/maitreepatel/Desktop/R/scRNA analysis/GSE214544_GSM6611298 (3)")
post.seurat <- CreateSeuratObject(counts = post.exercise,
                                  project = "muscle_post")

#creating a list
muscle.list <- list(pre.seurat,
                    post.seurat)

#merging the data
merged.muscle <- merge(x = pre.seurat,
                       y = post.seurat,
                       add.cell.ids = c("pre", "post"))

#getting the metadata
library(stringr)

sample <- names(merged.muscle@active.ident)
head(sample)

sample.detect <- ifelse(str_detect(sample, "pre"),
                        "pre",
                        "post")

#addding a column to the metadata with the conditions
merged.muscle@meta.data$sample <- sample.detect

#manipulating identity class?
Idents(object = merged.muscle) <- "sample"

#normalizing dataset
#dividing the datasets as pre or post
muscle.list <- SplitObject(merged.muscle,
                           split.by = "sample")

#performing pre-processing on both the elements in list using a loop
for (i in 1:length(muscle.list)) {
  muscle.list[[i]] <- NormalizeData(muscle.list[[i]], 
                                    verbose = FALSE)
  
  muscle.list[[i]] <- subset(muscle.list[[i]], 
                             downsample = 1000)

  muscle.list[[i]] <- FindVariableFeatures(muscle.list[[i]], 
                                           selection.method = "vst",
                                           nfeatures = 2000, verbose = FALSE)
}
#use lapply!!
#especially with more than two datasets

#selecting features that are variable across the datasets
features.var <- SelectIntegrationFeatures(object.list = muscle.list)

#running PCA on highly variable features
muscle.list <- lapply(X = muscle.list, FUN = function(x) {
  x <- ScaleData(x, 
                 features = features.var,
                 verbose = FALSE)
  
  x <- RunPCA(x, 
              
              features = features.var,
              verbose = FALSE)
})

#finding anchor genes and integrating data
anchors <- FindIntegrationAnchors(object.list = muscle.list)

merged.muscle <- IntegrateData(anchorset = anchors)

#data manipulation and cleaning
#standard workflow
merged.muscle <- ScaleData(merged.muscle,
                           verbose = FALSE)

merged.muscle <- FindVariableFeatures(merged.muscle,
                                      selection.method = "vst",
                                      nfeatures = 2000,
                                      verbose = FALSE)
#dimensionality reduction
merged.muscle <- RunPCA(merged.muscle,
                        npcs = 30,
                        verbose = FALSE)

merged.muscle <- RunUMAP(merged.muscle,
                         reduction = "pca",
                         dims = 1:20)

#clustering
merged.muscle <- FindNeighbors(merged.muscle,
                               reduction = "pca",
                               dims = 1:20)

merged.muscle <- FindClusters(merged.muscle,
                              resolution = 0.5)

DimPlot(merged.muscle,
        reduction = "umap")

Idents(object = merged.muscle) <- "sample"

DimPlot(merged.muscle,
        reduction = "umap")

#finding markers
Idents(object = merged.muscle) <- "sample"

#finding differentially expressed genes
sammple.markers <- FindMarkers(merged.muscle, 
                               ident.1 = "pre", 
                               ident.2 = "post")
head(sammple.markers)

library(ComplexHeatmap)

heatmapdf <- sammple.markers[1:25,]

row_ha = rowAnnotation("pre" = anno_barplot(heatmapdf$pct.1),
                       "post"= anno_barplot(heatmapdf$pct.2),
                       width = unit(10, "cm"))

ht0 <- Heatmap(heatmapdf$avg_log2FC,
               name = "Log2FC",
               cluster_rows = TRUE, 
               row_labels = rownames(heatmapdf), 
               right_annotation = row_ha,
               width = unit(1, "cm"))
ht0

#finding differentially expressed genes between cluster 0 and the rest
Idents(object = merged.muscle) <- "seurat_clusters"

cluster0.markers <- FindMarkers(merged.muscle, 
                                ident.1 = "0", 
                                ident.2 = NULL, 
                                only.pos = TRUE)
#heatmap cluster 0 versus others
head(cluster0.markers)

heatmapdf <- cluster0.markers[1:25,]

row_ha = rowAnnotation("Cluster 0" = anno_barplot(heatmapdf$pct.1),
                       "Others"= anno_barplot(heatmapdf$pct.2),
                       width = unit(10, "cm"))

ht1 <- Heatmap(heatmapdf$avg_log2FC,
               name = "Log2FC",
               cluster_rows = TRUE, 
               row_labels = rownames(heatmapdf), 
               right_annotation = row_ha,
               width = unit(1, "cm"))

ht1


DefaultAssay(merged.muscle) <- "RNA"

gc()

merged.muscle[['groups']] <- sample(x = c('pre', 'post'), 
                                    size = ncol(x = merged.muscle), 
                                    replace = TRUE)
merged.muscle$groups[1] <- "pre"
merged.muscle$groups[2] <- "post"

FindConservedMarkers(merged.muscle, 
                     ident.1 = 0, 
                     ident.2 = 1, 
                     grouping.var = "groups")

str(merged.muscle)
