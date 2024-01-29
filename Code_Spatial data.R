library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)

setwd("/Users/mcui/Documents/R-studio/Visium/PD1_project/PD1KO")
### slide PDHet
PDHet_1_dir <- "/Users/mcui/Documents/R-studio/Visium/PD1_project/PD1KO/PD1Het_1"
list.files(PDHet_1_dir)
PDHet_1 <- Load10X_Spatial(data.dir= PDHet_1_dir, filename = "filtered_feature_bc_matrix.h5")
## data processing
plot1 <- VlnPlot(PDHet_1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(PDHet_1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

### use SCTransform to normalize data
PDHet_1<- SCTransform(PDHet_1, assay = "Spatial", verbose = FALSE)
### gene visualization
SpatialFeaturePlot(PDHet_1, features = c("Dmd", "Col1a1","Cdh5","Npr3"))
SpatialFeaturePlot(PDHet_1, features = c("Msln","Lyz2","Pdgfrb","Dach1","Cd274","Pdcd1"))
### Dimensionality reduction, clustering, and visualization
PDHet_1 <- RunPCA(PDHet_1, assay = "SCT", verbose = FALSE)
PDHet_1 <- FindNeighbors(PDHet_1, reduction = "pca", dims = 1:30)
PDHet_1 <- FindClusters(PDHet_1, verbose = FALSE,resolution = 2)
PDHet_1 <- RunUMAP(PDHet_1, reduction = "pca", dims = 1:30)
### plot
p1 <- DimPlot(PDHet_1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(PDHet_1, label = TRUE, label.size = 3)
p1 + p2


### slide PDKO
PD1KO_dir <- "/Users/mcui/Documents/R-studio/Visium/PD1_project/PD1KO/PD1KO"
list.files(PD1KO_dir)
PD1KO <- Load10X_Spatial(data.dir= PD1KO_dir, filename = "filtered_feature_bc_matrix.h5")
## data processing
plot1 <- VlnPlot(PD1KO, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(PD1KO, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

### use SCTransform to normalize data
PD1KO<- SCTransform(PD1KO, assay = "Spatial", verbose = FALSE)
### gene visualization
SpatialFeaturePlot(PD1KO, features = c("Dmd", "Col1a1","Cdh5","Npr3"))
SpatialFeaturePlot(PD1KO, features = c("Msln","Lyz2","Pdgfrb","Dach1","Cd274","Pdcd1"))
### Dimensionality reduction, clustering, and visualization
PD1KO <- RunPCA(PD1KO, assay = "SCT", verbose = FALSE)
PD1KO <- FindNeighbors(PD1KO, reduction = "pca", dims = 1:30)
PD1KO <- FindClusters(PD1KO, verbose = FALSE,resolution = 2)
PD1KO <- RunUMAP(PD1KO, reduction = "pca", dims = 1:30)
### plot
p1 <- DimPlot(PD1KO, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(PD1KO, label = TRUE, label.size = 3)
p1 + p2
SpatialFeaturePlot(PD1KO, features = c("Col1a1"),crop = F)



PDHet[["orig.ident"]] <- "PDHet"
PDKO[["orig.ident"]] <- "PDKO"
DefaultAssay(PDHet_1) <- "SCT"
DefaultAssay(PDHet_2) <- "SCT"
DefaultAssay(PD1KO) <- "SCT"

### merge samples
PD <- merge(PDHet, PD1KO)

DefaultAssay(PD) <- "SCT"
PD<- SCTransform(PD, assay = "SCT", verbose = FALSE)


### plot fibrosis marker Postn

SpatialFeaturePlot(PD, features = c("Postn"),crop = F)

## cell type deconvolution see code for Cell2Location
