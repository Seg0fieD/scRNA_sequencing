# Library _____________________________________________________________________________
suppressPackageStartupMessages({
  library(dplyr)
  library(spatstat.core)
  library(Seurat)
  library(patchwork)
  library(DoubletFinder)
  library(SingleR)
  library(enrichR)
  library(Matrix)
  library(ggplot2)
  library(gridExtra)
  library(cowplot)
  library(SingleCellExperiment)
  library(SeuratWrappers)
  library(tidyverse)
  library(monocle3)
  library(celldex)
  library(EnhancedVolcano)
})
# ____________________________________________________________________________________

# set working dir 
setwd("~/scRNA_sequencing/")

# file path to save rds files
file_path <- "saved_rds/"

# Read files
bmmc1 <- readRDS("data/GSM4138872_scRNA_BMMC_D1T1.rds")
bmmc2 <- readRDS("data/GSM4138873_scRNA_BMMC_D1T2.rds")
cd1 <- readRDS("data/GSM4138874_scRNA_CD34_D2T1.rds")
cd2 <- readRDS("data/GSM4138875_scRNA_CD34_D3T1.rds")

# Check dimensions
dim(bmmc1)
dim(bmmc2)
dim(cd1)
dim(cd2)

# Create Seurat object
bmmc1s <- CreateSeuratObject(counts = bmmc1, project = "GSM4138872_scRNA_BMMC_D1T1")
bmmc2s <- CreateSeuratObject(counts = bmmc2, project = "GSM4138873_scRNA_BMMC_D1T2")
cd1s <- CreateSeuratObject(counts = cd1, project = "GSM4138874_scRNA_CD34_D2T1")
cd2s <- CreateSeuratObject(counts = cd2, project = "GSM4138875_scRNA_CD34_D3T1")

# check features and column
bmmc1s[[]]
bmmc2s[[]]
cd1s[[]]
cd2s[[]]

# Adding Metadeta
view(bmmc1s@meta.data)
view(bmmc2s@meta.data)
view(cd1s@meta.data)
view(cd2s@meta.data)

bmmc1s <- AddMetaData(bmmc1s, c("D1", "T1", "F"), col.name = c("Donor", "Replicate", "Sex"))
bmmc2s <- AddMetaData(bmmc2s, c("D1", "T2", "F"), col.name = c("Donor", "Replicate", "Sex"))
cd1s <- AddMetaData(cd1s, c("D2", "T1", "M"), col.name = c("Donor", "Replicate", "Sex"))
cd2s <- AddMetaData(cd2s, c("D3", "T1", "F"), col.name = c("Donor", "Replicate", "Sex"))



#dataframe
bmmc1_df <- data.frame(bmmc1s@meta.data)
bmmc2_df <- data.frame(bmmc2s@meta.data)
cd1_df <- data.frame(cd1s@meta.data)
cd2_df <- data.frame(cd2s@meta.data)



# Calculate Mitochondrial Percentage and add to metadata

bmmc1s[["percent.mt"]] <- PercentageFeatureSet(bmmc1s, pattern = "^MT")
bmmc2s[["percent.mt"]] <- PercentageFeatureSet(bmmc2s, pattern = "^MT")
cd1s[["percent.mt"]] <- PercentageFeatureSet(cd1s, pattern = "^MT")
cd2s[["percent.mt"]] <- PercentageFeatureSet(cd2s, pattern = "^MT")

view(bmmc1s@meta.data)
view(bmmc2s@meta.data)
view(cd1s@meta.data)
view(cd2s@meta.data)

# pre-filteration Visulization
## violin plot
png(file = "results/img/pre_filter_violiinPLT.png", width = 800, height = 700)
grid.arrange(
  VlnPlot(bmmc1s, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3),
  VlnPlot(bmmc2s, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3),
  VlnPlot(cd1s, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3),
  VlnPlot(cd2s, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3), 
  nrow = 2, ncol = 2)
dev.off()


# list of Seurat objects 
seurat_obj <- list (
  bmmc1s = bmmc1s,
  bmmc2s = bmmc2s,
  cd1s   = cd1s,
  cd2s   = cd2s
)

# Filtration 
bmmc1s_fl <- subset(bmmc1s, subset = nFeature_RNA  > 200 & nFeature_RNA < 2500  & nCount_RNA  < 5500 &  percent.mt < 0.5)
bmmc2s_fl <- subset(bmmc2s, subset = nFeature_RNA  > 200 & nFeature_RNA < 2500  & nCount_RNA  < 5500 &  percent.mt < 0.5)
cd1s_fl <- subset(cd1s, subset = nFeature_RNA  > 200 & nFeature_RNA < 2500  & nCount_RNA  < 5500 &  percent.mt < 0.5)
cd2s_fl <- subset(cd2s, subset = nFeature_RNA  > 200 & nFeature_RNA < 2500  & nCount_RNA  < 5500 &  percent.mt < 0.5)


# post filteration visualization
png(file = "results/img/post_filtr_violin.png", width = 800, height = 700)
grid.arrange(VlnPlot(bmmc1s_fl, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3),
             VlnPlot(bmmc2s_fl, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3),
             VlnPlot(cd1s_fl, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3),
             VlnPlot(cd2s_fl, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3), 
             nrow = 2, ncol =2)
dev.off()



# Normalization 
bmmc1s_fl <- NormalizeData(bmmc1s_fl)
bmmc2s_fl <- NormalizeData(bmmc2s_fl)
cd1s_fl <- NormalizeData(cd1s_fl)
cd2s_fl <- NormalizeData(cd2s_fl)


# Feature Selection 
## Variation Identification
bmmc1s_fl <- FindVariableFeatures(bmmc1s_fl, selection.method = "vst", nfeatures = 2000)
bmmc2s_fl <- FindVariableFeatures(bmmc2s_fl, selection.method = "vst", nfeatures = 2000)
cd1s_fl<- FindVariableFeatures(cd1s_fl, selection.method = "vst", nfeatures = 2000)
cd2s_fl<- FindVariableFeatures(cd2s_fl, selection.method = "vst", nfeatures = 2000)

## check variable features
head(VariableFeatures((bmmc1s_fl), 10))
head(VariableFeatures((bmmc2s_fl), 10))
head(VariableFeatures((cd1s_fl), 10))
head(VariableFeatures((cd2s_fl), 10))

### plot variable features 
png(file = "results/img/Variable_feature_Plot.png", height = 1000, width = 1500)
grid.arrange(LabelPoints(plot = VariableFeaturePlot(bmmc1s_fl),
                         points = head(VariableFeatures(bmmc1s_fl), 15),
                         repel = TRUE),
             LabelPoints(plot = VariableFeaturePlot(bmmc2s_fl),
                         points = head(VariableFeatures(bmmc2s_fl), 15),
                         repel = TRUE),
             LabelPoints(plot = VariableFeaturePlot(cd1s_fl),
                         points = head(VariableFeatures(cd1s_fl), 15),
                         repel = TRUE),
             LabelPoints(plot = VariableFeaturePlot(cd2s_fl),
                         points = head(VariableFeatures(cd2s_fl), 15),
                         repel = TRUE),
              ncol = 2, nrow = 2
             )
dev.off()


# Scaling --
bmmc1s_fl <- ScaleData(bmmc1s_fl, features =  rownames(bmmc1s_fl))
bmmc2s_fl <- ScaleData(bmmc2s_fl, features =  rownames(bmmc2s_fl))
cd1s_fl <- ScaleData(cd1s_fl, features =  rownames(cd1s_fl))
cd2s_fl <- ScaleData(cd2s_fl, features =  rownames(cd2s_fl))



# Dimensionality reduction ---
## PCA--
bmmc1s_fl <- RunPCA(bmmc1s_fl, features = VariableFeatures(object = bmmc1s_fl))
bmmc2s_fl <- RunPCA(bmmc2s_fl, features = VariableFeatures(object = bmmc2s_fl))
cd1s_fl <- RunPCA(cd1s_fl, features = VariableFeatures(object = cd1s_fl))
cd2s_fl <- RunPCA(cd2s_fl, features = VariableFeatures(object =  cd2s_fl))

## PCA results --
print(bmmc1s_fl[["pca"]], dims = 1:5, nfeatures = 10)
print(bmmc2s_fl[["pca"]], dims = 1:5, nfeatures = 10)
print(cd1s_fl[["pca"]], dims = 1:5, nfeatures = 10)
print(cd2s_fl[["pca"]], dims = 1:5, nfeatures = 10)


## PCA PLot 
png(file = "results/img/PCA_dim_plot.png", height = 1000, width = 1200)
grid.arrange(DimPlot(bmmc1s_fl, reduction = "pca", cols = "red"),
            DimPlot(bmmc2s_fl, reduction = "pca", cols = "blue"),
            DimPlot(cd1s_fl, reduction = "pca", cols = "green"),
            DimPlot(cd2s_fl, reduction = "pca", cols = "yellow"), 
            nrow = 2, ncol=2)
dev.off()


## Dimension determination - ELBOW Plot
png(file = "results/img/Elbow_plot.png", height = 900, width = 1200)
grid.arrange(ElbowPlot(bmmc1s_fl),
             ElbowPlot(bmmc2s_fl),
             ElbowPlot(cd1s_fl),
             ElbowPlot(cd2s_fl),
             nrow = 2, ncol = 2)
dev.off()


# Clustering --
## Find Neighbours 
bmmc1s_fl <- FindNeighbors(bmmc1s_fl,  dims = 1:15)
bmmc2s_fl <- FindNeighbors(bmmc2s_fl, dims = 1:15)
cd1s_fl <- FindNeighbors(cd1s_fl,  dims = 1:15)
cd2s_fl <- FindNeighbors(cd2s_fl, dims = 1:15)

  
## Clustering based on resolution
bmmc1s_fl <- FindClusters(bmmc1s_fl, resolution = 0.1)
bmmc2s_fl <- FindClusters(bmmc2s_fl, resolution = 0.1)
cd1s_fl <- FindClusters(cd1s_fl, resolution = 0.1)
cd2s_fl <- FindClusters(cd2s_fl, resolution = 0.1)


view(bmmc1s_fl@meta.data)
view(bmmc2s_fl@meta.data)
view(cd1s_fl@meta.data)
view(cd2s_fl@meta.data)

head(bmmc1s_fl@meta.data)

## Plot for cluster 
blk <- theme_minimal() +
  theme(plot.background = element_rect(fill = "black"),
        panel.background = element_rect(fill = "black"),
        legend.background = element_rect(fill = "black"),
        text = element_text(color = "white"),
        axis.text = element_text(color = "white"),
        axis.line = element_line(color = "white"),  # Set axis lines to white
        panel.grid = element_blank())


png(file = "results/img/Cluster.png", height = 900, width = 1200)
grid.arrange(DimPlot(bmmc1s_fl, group.by = "RNA_snn_res.0.1", label = TRUE) + blk,
             DimPlot(bmmc2s_fl, group.by = "RNA_snn_res.0.1", label = TRUE) + blk,
             DimPlot(cd1s_fl, group.by = "RNA_snn_res.0.1", label = TRUE) + blk,
             DimPlot(cd2s_fl, group.by = "RNA_snn_res.0.1", label = TRUE) + blk,
             nrow = 2, ncol = 2)
dev.off()



## Set identity of cluster based on selected resolution
Idents(bmmc1s_fl)
Idents(bmmc1s_fl) <- "RNA_snn_res.0.1" 
Idents(bmmc1s_fl)

Idents(bmmc2s_fl)
Idents(bmmc2s_fl) <- "RNA_snn_res.0.1"
Idents(bmmc2s_fl)

Idents(cd1s_fl)
Idents(cd1s_fl) <- "RNA_snn_res.0.1"
Idents(cd1s_fl)

Idents(cd2s_fl)
Idents(cd2s_fl) <- "RNA_snn_res.0.1"
Idents(cd2s_fl)





## Non_linear Dimensionality Reduction
bmmc1s_fl <- RunUMAP(bmmc1s_fl, dims = 1:15)
bmmc2s_fl <- RunUMAP(bmmc2s_fl, dims = 1:15)
cd1s_fl<- RunUMAP(cd1s_fl, dims = 1:15)
cd2s_fl<- RunUMAP(cd2s_fl, dims = 1:15)


## Individual Cluster visualization 
png(file =  "results/img/individual_cluster.png", width = 1200, height = 800)
grid.arrange(DimPlot(bmmc1s_fl, reduction = "umap", label = TRUE) + blk, 
             DimPlot(bmmc2s_fl, reduction = "umap") + blk,
             DimPlot(cd1s_fl, reduction = "umap") + blk, 
             DimPlot(cd2s_fl, reduction = "umap") + blk, 
             nrow=2, ncol=2)
dev.off()



# Doublet Removal
## pK identification (no ground truth)
bmmc1_dr <- find.pK(summarizeSweep(paramSweep(bmmc1s_fl, PCs = 1:15, sct = FALSE), GT = FALSE))
bmmc2_dr <- find.pK(summarizeSweep(paramSweep(bmmc2s_fl, PCs = 1:15, sct = FALSE), GT = FALSE))
cd1_dr <- find.pK(summarizeSweep(paramSweep(cd1s_fl, PCs = 1:15, sct = FALSE), GT = FALSE))
cd2_dr <- find.pK(summarizeSweep(paramSweep(cd2s_fl, PCs = 1:15, sct = FALSE), GT = FALSE))

## pK indentification ggplot
png(filename = "results/img/pk_plot.png", height = 1000, width = 1800)
grid.arrange(ggplot(bmmc1_dr, aes(pK, BCmetric, group = 1)) + geom_point()+ geom_line(),
             ggplot(bmmc2_dr, aes(pK, BCmetric, group = 1)) + geom_point()+ geom_line(),
             ggplot(cd1_dr, aes(pK, BCmetric, group = 1)) + geom_point()+ geom_line(),
             ggplot(cd2_dr, aes(pK, BCmetric, group = 1)) + geom_point()+ geom_line(), 
             nrow = 2, ncol = 2)
dev.off()


## filter....
bmmc1_pk <- bmmc1_dr %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
bmmc1_pk <- as.numeric(as.character(bmmc1_pk[1]))

bmmc2_pk <- bmmc2_dr %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
bmmc2_pk <- as.numeric(as.character(bmmc2_pk[1]))

cd1_pk <- cd1_dr %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
cd1_pk <- as.numeric(as.character(cd1_pk[1]))

cd2_pk <- cd2_dr %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
cd2_pk <- as.numeric(as.character(cd2_pk[1]))



# Homotypic Doublet Proportion Estimation (_hmpr)
bmmc1_hmpr <- modelHomotypic(bmmc1s_fl@meta.data$seurat_clusters)
bmmc1_nExp <- round(0.075 * nrow(bmmc1s_fl@meta.data))
bmmc1_nExp_adj <- round(bmmc1_nExp * (1 - bmmc1_hmpr))

bmmc2_hmpr <- modelHomotypic(bmmc2s_fl@meta.data$seurat_clusters)
bmmc2_nExp <- round(0.075 * nrow(bmmc2s_fl@meta.data))
bmmc2_nExp_adj <- round(bmmc2_nExp * (1 - bmmc2_hmpr))

cd1_hmpr <- modelHomotypic(cd1s_fl@meta.data$seurat_clusters)
cd1_nExp <-round(0.075 * nrow(cd1s_fl@meta.data))
cd1_nExp_adj <- round(cd1_nExp * (1 - cd1_hmpr))

cd2_hmpr <- modelHomotypic(cd2s_fl@meta.data$seurat_clusters)
cd2_nExp <-round(0.075 * nrow(cd2s_fl@meta.data))
cd2_nExp_adj <- round(cd2_nExp * (1 - cd2_hmpr))


# Doublet finder 
## nExp
bmmc1_dbf <- doubletFinder(bmmc1s_fl, PCs = 1:15 , pN = 0.25, 
                           pK = bmmc1_pk, nExp = bmmc1_nExp, 
                           reuse.pANN = FALSE, sct = FALSE)

bmmc2_dbf <- doubletFinder(bmmc2s_fl, PCs = 1:15 , pN = 0.25, 
                           pK = bmmc2_pk, nExp = bmmc2_nExp, 
                           reuse.pANN = FALSE, sct = FALSE)

cd1_dbf <- doubletFinder(cd1s_fl, PCs = 1:15 , pN = 0.25, 
                           pK = cd1_pk, nExp = cd1_nExp, 
                           reuse.pANN = FALSE, sct = FALSE)

cd2_dbf <- doubletFinder(cd2s_fl, PCs = 1:15 , pN = 0.25, 
                         pK = cd2_pk, nExp = cd2_nExp, 
                         reuse.pANN = FALSE, sct = FALSE)


## nExp _adj
bmmc1_dbf <- doubletFinder(bmmc1_dbf, PCs = 1:15 , pN = 0.25, 
                           pK = bmmc1_pk, nExp = bmmc1_nExp_adj, 
                           reuse.pANN = FALSE, sct = FALSE)

bmmc2_dbf <- doubletFinder(bmmc2_dbf, PCs = 1:15 , pN = 0.25, 
                           pK = bmmc2_pk, nExp = bmmc2_nExp_adj, 
                           reuse.pANN = FALSE, sct = FALSE)

cd1_dbf <- doubletFinder(cd1_dbf, PCs = 1:15 , pN = 0.25, 
                         pK = cd1_pk, nExp = cd1_nExp_adj, 
                         reuse.pANN = FALSE, sct = FALSE)

cd2_dbf <- doubletFinder(cd2_dbf, PCs = 1:15 , pN = 0.25, 
                         pK = cd2_pk, nExp = cd2_nExp_adj, 
                         reuse.pANN = FALSE, sct = FALSE)

# dataframe update 
view(bmmc1_dbf@meta.data)
view(bmmc2_dbf@meta.data)
view(cd1_dbf@meta.data)
view(cd2_dbf@meta.data)

bmmc1_df <- data.frame(bmmc1_dbf@meta.data)
bmmc2_df <- data.frame(bmmc2_dbf@meta.data) 
cd1_df <- data.frame(cd1_dbf@meta.data)
cd2_df <- data.frame(cd2_dbf@meta.data)


## doublets visualization
png(file =  "results/img/Doublet _visualization_bmmc1.png", width = 1200, height = 850)
grid.arrange(DimPlot(bmmc1_dbf, reduction = 'umap', group.by = "pANN_0.25_2_404"),
             DimPlot(bmmc1_dbf, reduction = 'umap', group.by = "pANN_0.25_2_302", cols = c("red")),
             DimPlot(bmmc1_dbf, reduction = 'umap', group.by = "DF.classifications_0.25_2_404", 
                     cols = c("black", "#28E2E5")),
             DimPlot(bmmc1_dbf, reduction = 'umap', group.by = "DF.classifications_0.25_2_302", 
                     cols = c("black", "#2297E6")))
dev.off()


png(file =  "results/img/Doublet _visualization_bmmc2.png", width = 1200, height = 850)
grid.arrange(DimPlot(bmmc2_dbf, reduction = 'umap', group.by = "pANN_0.25_1_407"),
             DimPlot(bmmc2_dbf, reduction = 'umap', group.by = "pANN_0.25_1_299",
                     cols = c("red")),
             DimPlot(bmmc2_dbf, reduction = 'umap', group.by = "DF.classifications_0.25_1_407",
                     cols = c("black", "#28E2E5")),
             DimPlot(bmmc2_dbf, reduction = 'umap', group.by = "DF.classifications_0.25_1_299",
                     cols = c("black", "#61D04F")))
dev.off()

png(file =  "results/img/Doublet _visualization_cd2.png", width = 1200, height = 850)
grid.arrange(DimPlot(cd1_dbf, reduction = 'umap', group.by = "pANN_0.25_4_99",cols = c("blue")),
             DimPlot(cd1_dbf, reduction = 'umap', group.by = "pANN_0.25_4_68",cols = c("red")),
             DimPlot(cd1_dbf, reduction = 'umap', group.by = "DF.classifications_0.25_4_99"),
             DimPlot(cd1_dbf, reduction = 'umap', group.by = "DF.classifications_0.25_4_68", 
                     cols = c("#CD0BBC", "#61D04F")))
dev.off()


png(file =  "results/img/Doublet _visualization_cd3.png", width = 1200, height = 850)
grid.arrange(DimPlot(cd2_dbf, reduction = 'umap', group.by = "pANN_0.25_2_347",cols = c("orange")),
             DimPlot(cd2_dbf, reduction = 'umap', group.by = "pANN_0.25_2_216",cols = c("green")),
             DimPlot(cd2_dbf, reduction = 'umap', group.by = "DF.classifications_0.25_2_347"),
             DimPlot(cd2_dbf, reduction = 'umap', group.by = "DF.classifications_0.25_2_216",
                     cols = c("black", "#F5C710")))
dev.off()



### determine singlet and doublet    -------
table(bmmc1_dbf@meta.data$DF.classifications_0.25_2_404)
table(bmmc2_dbf@meta.data$DF.classifications_0.25_1_407)
table(cd1_dbf@meta.data$DF.classifications_0.25_4_99)
table(cd2_dbf@meta.data$DF.classifications_0.25_2_347)


# table(bmmc1_dbf@meta.data$DF.classifications_0.25_2_404)
# Doublet Singlet 
# 404    4978 
# 
# table(bmmc2_dbf@meta.data$DF.classifications_0.25_1_407)
# Doublet Singlet 
# 407    5018 
# 
# table(cd1_dbf@meta.data$DF.classifications_0.25_4_99)
# Doublet Singlet 
# 99    1225 
# 
# table(cd2_dbf@meta.data$DF.classifications_0.25_2_347)
# Doublet Singlet 
# 347    4275

# table(bmmc1_dbf@meta.data$DF.classifications_0.25_2_302)
# Doublet Singlet 
# 302      5080 

# table(bmmc2_dbf@meta.data$DF.classifications_0.25_1_299)
# Doublet Singlet 
# 299      5126 

# table(cd1_dbf@meta.data$DF.classifications_0.25_4_68) 
# Doublet Singlet 
# 68       1256 

# table(cd2_dbf@meta.data$DF.classifications_0.25_2_216)
# Doublet Singlet 
# 216      4406 


# save seurat objects in rds format
saveRDS(bmmc1_dbf, file = paste0(file_path, "bmmc1.rds"))
saveRDS(bmmc2_dbf, file = paste0(file_path, "bmmc2.rds"))
saveRDS(cd1_dbf, file = paste0(file_path, "cd1.rds"))
saveRDS(cd2_dbf, file = paste0(file_path, "cd2.rds"))

# ------------------------- Merging(WITH-N0-BATCH-CORRECTION) ----------------------------------.
merg_no_batch <- merge(x = bmmc1_dbf, y = list(bmmc2_dbf, cd1_dbf, cd2_dbf), merge.data = FALSE)
view(merg_no_batch@meta.data)


## find variable features 
merg_no_batch <- FindVariableFeatures(merg_no_batch, selection.method = "vst", nfeatures = 2000)
head(VariableFeatures(merg_no_batch),15)

## variable feature plot 
png(file =  "results/img/variable_feature_NoBatchMerg.png",width = 1200, height = 800)
LabelPoints(plot = VariableFeaturePlot(merg_no_batch), 
            points = head(VariableFeatures(merg_no_batch), 20), 
            repel = T)
dev.off()

## Scale and normalize 
merg_no_batch <- NormalizeData(merg_no_batch)
merg_no_batch <- ScaleData(merg_no_batch, features = rownames(merg_no_batch))

# dimensionality reduction
merg_no_batch <- RunPCA(merg_no_batch)

print(merg_no_batch[['pca']], dims= 1:5, nfeatures = 15)

## plot
png(file = "results/img/PCA_NoBatchMerge.png", width = 800, height = 600)
DimPlot(merg_no_batch, reduction = "pca") + blk
dev.off()

## Elbow plot 
png(file = "results/img/Elbow_NoBatchMerge.png", width = 800, height = 600)
ElbowPlot(merg_no_batch)
dev.off()

## Find Neighbour 
merg_no_batch <- FindNeighbors(merg_no_batch, dims = 1:20)

## find clusters 
merg_no_batch <- FindClusters(merg_no_batch)

## UMAP reduction
merg_no_batch <- RunUMAP(merg_no_batch, dims = 1:20)

## UMAP Plot
png(file = "results/img/UMAP_NoBatchMerg.png", width = 1000, height = 700)
DimPlot(object = merg_no_batch, reduction = "umap") + DarkTheme() + 
  labs(title ="\t\t\t\t   Cluster of Merged Seurat objects with no Batch Correction redcution: UMAP")
dev.off()


# Saved the seurat objects 

saveRDS(merg_no_batch, file = paste0(file_path, "Merged_seurat_object_with_NO_batch_correction.rds"))




# --------------------------------- Merging Seurat Objects (WITH-BATCH-CORRECTION) ---------------------------------------.

##--load--data--
bmmc1_dbf <- readRDS(file.path(file_path,"bmmc1_postprocessed_doubletFinder.rds"))
bmmc2_dbf <- readRDS(file.path(file_path,"bmmc2_postprocessed_doubletFinder.rds"))
cd1_dbf   <- readRDS(file.path(file_path,"cd1_postprocessed_doubletFinder.rds"))
cd2_dbf   <- readRDS(file.path(file_path,"cd2_postprocessed_doubletFinder.rds"))

#view(bmmc1_dbf@meta.data)

#add-meta_data/ creating a coloum "sample" in the metadata
bmmc1_mr <-AddMetaData(bmmc1_dbf, c("BMMC_D1T1"),col.name = c("sample"))
bmmc2_mr <- AddMetaData(bmmc2_dbf, c("BMMC_D1T2"),col.name = c("sample")) 
cd1_mr <- AddMetaData(cd1_dbf, c("CD34_D2T1"),col.name = c("sample"))
cd2_mr <- AddMetaData(cd2_dbf, c("CD34_D3T1"), col.name = c("sample"))

#view(bmmc2_mr@meta.data)

merg_wBatch <- merge(x = bmmc1_mr, y = list(bmmc2_mr, cd1_mr, cd2_mr), merge.data = FALSE)
view(merg_wBatch@meta.data)


## Split 
merg_wBatch_split <- SplitObject(merg_wBatch, split.by = "sample")


## SCTransform
merg_wBatch_split <- lapply(X = merg_wBatch_split, FUN = SCTransform)


# feature selection that are continuously variable in the data
merg_wBatch_featur <- SelectIntegrationFeatures(object.list =  merg_wBatch_split)


# Identification of anchor 
merg_wBatch_split <- PrepSCTIntegration(object.list = merg_wBatch_split, 
                                        anchor.features = merg_wBatch_featur)

merg_wBatch_anchor <- FindIntegrationAnchors(object.list = merg_wBatch_split, 
                                             normalization.method = "SCT",
                                            anchor.features = merg_wBatch_featur)
# Integrated Data Assay creation 
merg_wBatch_intgrt <- IntegrateData(anchorset = merg_wBatch_anchor, normalization.method = "SCT")


# Computing PCA for Integrated Assay 
merg_wBatch_intgrt <- RunPCA(merg_wBatch_intgrt)

# Elbow Plot 
png(file = "results/img/Elbow_Merg_with_batch_correction.png", width = 1000, height = 700)
ElbowPlot(merg_wBatch_intgrt)
dev.off()

# UMAP reduction 
merg_wBatch_intgrt <- RunUMAP(merg_wBatch_intgrt, reduction = "pca", dims = 1:20) ####

view(merg_wBatch_intgrt@meta.data)


## plots 
plt1 <- DimPlot(merg_wBatch_intgrt, reduction = "umap", 
                group.by = "sample") + DarkTheme() + 
  labs(title ="Cluster grouped by sample with batch effect, redution: UMAP")

plt2 <- DimPlot(merg_wBatch_intgrt, reduction = "umap", group.by = "seurat_clusters", 
                label = T, repel = T) + DarkTheme() + 
  labs(title ="Clustering based on the Seurat clusters, reduction: UMAP")

png(file = "results/img/batch_corrected_clusters_comparison.png", width = 2000, height = 900)
grid.arrange(plt1, plt2, ncol = 2)
dev.off()


## Clustering
merg_wBatch_intgrt <- FindNeighbors(merg_wBatch_intgrt, dims = 1:20)
merg_wBatch_intgrt <- FindClusters(merg_wBatch_intgrt, resolution = 0.1)

### plot 
png(file = "results/img/cluster_integrated_assay_with_batch_correction.png", width = 1400, height = 800)
DimPlot(merg_wBatch_intgrt, label =  TRUE) +DarkTheme() + 
  labs(title = "\t \t \t \t Cluster of Integrated Assay with Batch Correction")
dev.off()

# save batch corrected seurat object
saveRDS(merg_wBatch_intgrt, file = paste0(file_path, "Integrated_merged_samples_with_Batch_correction.rds"))

# -------------------------- Automatic-Cell-Type-Annotation ---------------------------------------------.
## load data
merg_wBatch_intgrt <- readRDS(file.path(file_path,"Integrated_merged_samples_with_Batch_correction.rds"))

## Automatic Annotation
hpcad <- celldex::HumanPrimaryCellAtlasData()
hpcad 

scE <- as.SingleCellExperiment(DietSeurat(merg_wBatch_intgrt))
scE


hpcad_main <- SingleR(test = scE, assay.type.test = 1, ref = hpcad, labels = hpcad$label.main)
hpcad_fine <- SingleR(test = scE, assay.type.test = 1, ref = hpcad, labels = hpcad$label.fine)

table(hpcad_main$pruned.labels)
table(hpcad_fine$pruned.labels)

merg_wBatch_intgrt@meta.data$hpcad_main <- hpcad_main$pruned.labels
merg_wBatch_intgrt@meta.data$hpcad_fine <- hpcad_fine$pruned.labels

view(merg_wBatch_intgrt@meta.data)

## set identity 
auto_annot <- merg_wBatch_intgrt
auto_annot <- SetIdent(auto_annot, value = "hpcad_main")



### plot img
png(file = "results/img/automatic_Cell_Type_Annotation.png", 
    width = 1400, height = 800)
DimPlot(auto_annot, reduction = "umap", group.by = "hpcad_main", 
        label = TRUE) + DarkTheme() + ggtitle("\t\t Automatic Cell-Type Annotation \n")
dev.off()

## save automatic annotations 
saveRDS(auto_annot, file = file.path(file_path, "automatic_annotated_obj.rds"))


# -------------------------- Manual-Cell-Type-Annotation ---------------------------------------------.
## Determine differential expressed gene from manual annotations------.
DE_markers <- FindAllMarkers(auto_annot, only.pos = TRUE, 
                             min.pct = 0.1, logfc.threshold = 0.25)

# saveRDS(DE_markers, file = file.path(file_path, "differential_markers.rds"))
DE_markers <- readRDS(file.path(file_path,"differential_markers.rds"))

## create a list of marker genes from the provided table
cell_markers <- list(
  "HSC" = c("CD34", "CD38", "Sca1", "Kit"),
  "LMPP" = c("CD38", "CD52", "CSF3R", "Kit", "CD34", "Flk2"),
  "CLP" = "IL7R",
  "GMP" = "ELANE",
  "CMP" = c("IL3", "GM-CSF", "M-CSF"),
  "B_cells" = c("CD19", "CD20", "CD38"),
  "Pre_B" = c("CD19", "CD34"),
  "Plasma" = c("SDC1", "IGHA1", "IGLC1", "MZB1", "JCHAIN"),
  "T_cells" = "CD3D",
  "CD8_T" = c("CD3D", "CD3E", "CD8A", "CD8B"),
  "CD4_T" = c("CD3D", "CD3E", "CD4"),
  "NK" = c("FCGR3A", "NCAM1", "NKG7", "KLRB1"),
  "Erythrocytes" = c("GATA1", "HBB", "HBA1", "HBA2"),
  "pDC" = c("IRF8", "IRF4", "IRF7"),
  "cDC" = c("CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA"),
  "CD14_Mono" = c("CD14", "CCL3", "CCL4", "IL1B"),
  "CD16_Mono" = c("FCGR3A", "CD68", "S100A12"),
  "Basophils" = "GATA2"
)


## Determination of average expression of marker genes in each cluster to assign labels 
expression_avg <- AverageExpression(auto_annot, features = unique(unlist(cell_markers)))
expression_avg$integrated


## assigning cell type to cluster based on maker expression pattern
new_clusterID <- rep("Unknown", 34)  # Initialize with 34 elements
names(new_clusterID) <- levels(auto_annot)

# cluster mapping based on the labels
new_clusterID[which(names(new_clusterID) == "T_cells")] <- "T_cells"
new_clusterID[which(names(new_clusterID) == "Monocyte")] <- "CD14_Mono"
new_clusterID[which(names(new_clusterID) == "B_cell")] <- "B_cells"
new_clusterID[which(names(new_clusterID) == "NK_cell")] <- "NK"
new_clusterID[which(names(new_clusterID) == "Pro-B_cell_CD34+")] <- "Pre_B"
new_clusterID[which(names(new_clusterID) == "Macrophage")] <- "Macrophages"
new_clusterID[which(names(new_clusterID) == "Myelocyte")] <- "Myelocyte"
new_clusterID[which(names(new_clusterID) == "GMP")] <- "GMP"
new_clusterID[which(names(new_clusterID) == "CMP")] <- "CMP"
new_clusterID[which(names(new_clusterID) == "DC")] <- "cDC"
new_clusterID[which(names(new_clusterID) == "HSC_CD34+")] <- "HSC"
new_clusterID[which(names(new_clusterID) == "Erythroblast")] <- "Erythrocytes"

auto_annot <- RenameIdents(auto_annot, new_clusterID)

view(auto_annot@meta.data)


## plots
png(file = "results/img/Manual_cell_type_annotations3.png", 
    width = 1400, height = 800)
DimPlot(auto_annot, reduction = "umap",label = TRUE) + DarkTheme() + ggtitle("\t\t Manual Cell-Type Annotation \n")
dev.off()

# Visualize marker gene expression
mrkr_genes <- c("CD34", "KLRB1", "IL7R", "NKG7")

# violin plot
png(file = "results/img/Violin_plot_Marker_gene_expression.png", 
    width = 2200, height = 1000)
VlnPlot(auto_annot, features = mrkr_genes, ncol = 4)
dev.off()

# UMAP feature plot
png(file = "results/img/UMAP_Marker_gene_expression.png", 
    width = 2200, height = 1000)
FeaturePlot(auto_annot, features = mrkr_genes, ncol = 4)
dev.off()

# Cell-Type Proportion per sample
cell_prop <- prop.table(table(Idents(auto_annot), auto_annot$sample))

# dataframe conversion
cell_prop_df <- as.data.frame(cell_prop)
colnames(cell_prop_df) <- c("Cell_Type", "Sample", "Proportion")

# plot
png(file = "results/img/cell_type_per_sample.png", 
    width = 1400, height = 800)
ggplot(cell_prop_df, aes(x = Sample, y = Proportion, fill = Cell_Type)) +
         geom_bar(stat =  "identity") + 
         theme_dark() +
         theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "bold"),  
               axis.text.y = element_text(size = 14, face = "bold"),  
               axis.title = element_text(size = 16, face = "bold"),   
               plot.title = element_text(size = 20, face = "bold"),   
               legend.text = element_text(size = 14, face = "bold"), 
               legend.title = element_text(size = 14, face = "bold")
               ) + ggtitle("\t \t \t\t\t\t\t\t\t\t Cell Type Proportions per Sample")
dev.off()

saveRDS(auto_annot, file = file.path(file_path, "Manual_annotated_obj.rds"))

# --------------------------------- Differential Expression Analysis ------------------------.
## volcano plot function def 
volcano_plt <- function(DE_results, plt_title){
  EnhancedVolcano(DE_results,lab = rownames(DE_results),
                  x = "avg_log2FC", y = "p_val_adj",
                  title = plt_title, 
                  pCutoff = 0.05,
                  FCcutoff = 0.58,             #LOG2(1.5)
                  pointSize = 3.0,
                  labSize = 4.0
  )
}

# Set identity 
Idents(auto_annot) <- auto_annot@meta.data$hpcad_main
unique(Idents(auto_annot))

## B_cell vs T_cell comparison
B_v_T <- FindMarkers(object = auto_annot, 
                     ident.1 = "B_cell",
                     ident.2 = "T_cells",
                     min.pct = 0.25,
                     test.use = "wilcox")


## T_cell vs Monocytes comparisons 
T_v_Mono <- FindMarkers(object = auto_annot, 
                        ident.1 = "T_cells",
                        ident.2 = "Monocyte",
                        min.pct = 0.25,
                        test.use = "wilcox")

## Volcano plots 
png(file = "results/img/Volcano_plots_B_Cells vs T_cells.png", 
    width = 1400, height = 800)
volcano_plt(B_v_T, "B_Cells vs T_cells")
dev.off()

png(file = "results/img/Volcano_plots_T_cells vs Monocytes.png", 
    width = 1400, height = 800)
volcano_plt(T_v_Mono, "T_Cells vs Monocytes")
dev.off()


## TOP 5 Differential Expressed genes plot 
### Filterating to ensure valid selection is performed
filtr_validGene <- function(DE_results, p_val_cutoff = 0.05){
  DE_results %>% 
    tibble::rownames_to_column("gene") %>% 
    filter(!is.na(p_val_adj) & p_val_adj < p_val_cutoff) %>% 
    arrange(p_val_adj) %>% 
    head(5)
}


## Using a value smaller than p_val_adj = 0 to avoid issues during plotting.
zero_handling <- function(DE_results){
  DE_results$p_val_adj[DE_results$p_val_adj == 0] <- 1e-10
  return(DE_results)
}


## Get top 5 gene (filtered and no zeros)
top5_B_v_T <- filtr_validGene(B_v_T)
top5_T_v_Mono <- filtr_validGene(T_v_Mono)

top5_B_v_T <- zero_handling(B_v_T)
top5_T_v_Mono <- zero_handling(T_v_Mono)

top5_B_v_T
top5_T_v_Mono

### combine data
BTM_plt <- do.call(rbind, list(
  data.frame(top5_B_v_T, comparison = "B_Cell vs T_Cells"),
  data.frame(top5_T_v_Mono, comparison = "T_Cells vs Monocyte")
))

BTM_plt
unique(auto_annot@meta.data$sample)

### dot plot
top5_dot_plt <- ggplot(BTM_plt, aes(x = comparison, y = gene, size = -log10(p_val_adj), color = avg_log2FC)) + 
  geom_point() + scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + 
  theme_dark() + labs(x = "Comparison", y = "Genes", size = "-log10(p_val_adj)", color = "log2 Fold Change") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggsave("results/img/TOP5_dot_plot.png", top5_dot_plt)





# --------------------------------------------- Pathway Analysis ------------------------------------------.
## Set ident
unique(auto_annot@meta.data$orig.ident)
Idents(auto_annot) <- auto_annot@meta.data$orig.ident

bmmc_v_cd34 <- FindMarkers(object = auto_annot,
                           ident.1 = "BMMC",
                           ident.2 = "CD34",
                           min.pct = 0.1,
                           test.use = "wilcox")
## --- TOP 5 DEGs
top5_overall <- bmmc_v_cd34 %>% tibble::rownames_to_column("gene") %>% 
                arrange(p_val_adj) %>% head(5) %>%
                select(gene, avg_log2FC, p_val_adj)

## Plot 
plt_ovrall <- ggplot(top5_overall, 
                     aes(x = reorder(gene, avg_log2FC), 
                         y = avg_log2FC)) +
  geom_bar(stat = "identity", aes(fill = -log10(p_val_adj))) +
  coord_flip() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(title = "Top 5 DEGs: Overall BMMC vs CD34",
       x = "Gene", 
       y = "Log2 Fold Change",
       fill = "-log10(adj.P)") +
  theme(axis.text.y = element_text(face = "italic"))

ggsave("results/img/DEG_overall_TOP_5.png", plt_ovrall, width = 12, height = 10,
       dpi = 300)

## For Monocytes only
### Set ident
Idents(auto_annot) <- auto_annot@meta.data$hpcad_main

monotypes <- WhichCells(auto_annot, idents = "Monocyte")

monotype_subset <- subset(auto_annot, cells = monotypes)

Idents(monotype_subset) <- "orig.ident"

## --- DEG analysis on Monocytes
Mono_bmmc_v_cd34 <- FindMarkers(object = monotype_subset,
                                ident.1 = "BMMC",
                                ident.2 = "CD34",
                                min.pct = 0.1,
                                test.use = "wilcox")

## Top 5 DEGS for Monocytes 
top5_mono <- Mono_bmmc_v_cd34 %>% tibble::rownames_to_column("gene") %>% 
              arrange(p_val_adj) %>%
              head(5) %>% select(gene, avg_log2FC, p_val_adj)

## Plot 
plt_mono <- ggplot(top5_mono, 
                   aes(x = reorder(gene, avg_log2FC), 
                       y = avg_log2FC)) +
  geom_bar(stat = "identity", aes(fill = -log10(p_val_adj))) +
  coord_flip() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(title = "Top 5 DEGs: Monocytes BMMC vs CD34",
       x = "Gene", 
       y = "Log2 Fold Change",
       fill = "-log10(adj.P)") +
  theme(axis.text.y = element_text(face = "italic"))

ggsave("results/img/DEG_MONOCYTES_TOP_5.png", plt_ovrall, width = 12, height = 10,
       dpi = 300)

## Pathway way analysis on group
Idents(auto_annot) <- "orig.ident"

enrich_plt <- DEenrichRPlot(object = auto_annot,
                ident.1 = "BMMC",
                ident.2 = "CD34",
                enrich.database = "GO_Biological_2001",
                max.genes = 1000,
                figure.size = c(10,8))

png(file = "results/img/pathway_analysis_enrichment_plot.png", 
    width = 3000, height = 2400, res = 300)
print(enrich_plt)
dev.off()





