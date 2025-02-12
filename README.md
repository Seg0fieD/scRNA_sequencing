# Single-Cel RNA Sequencing Analysis

This project focused on single-cell RNA sequencing (scRNA-seq) data analysis, specifically from human bone marrow cells and CD34+ enriched bone marrow cells. The analysis was performed using the Seurat package in R.

## Key Tasks:

1. **Loading the Data**: Expression matrices from two datasets (BMMC and CD34) were loaded and a Seurat object was constructed.
   
2. **Metadata Addition**: Relevant metadata was added to each sample, including the number of cells and genes in the expression matrices.

3. **Preprocessing**: Data preprocessing included filtering, doublet removal using DoubletFinder, normalization, and feature selection.

4. **Batch Correction**: Two different methods were applied: merging without batch correction and integrating the data using Seurat’s batch correction method. Comparisons were made to evaluate the necessity of batch correction.

5. **Dimensionality Reduction**: Principal Component Analysis (PCA) and Uniform Manifold Approximation and Projection (UMAP) were used to reduce dimensions and visualize the data in 2D.

6. **Clustering**: Clustering was performed on the reduced data, resulting in 7-15 clusters, which were visualized in 2D.

7. **Cell Type Annotation**: Both automatic and manual cell type annotations were performed. Differential expression analysis was conducted to identify cell-type-specific markers.

8. **Differential Expression Analysis**: Differential expression between cell types (B cells vs T cells, T cells vs Monocytes) was performed, and top differentially expressed genes were plotted.

9. **Pathway Analysis**: Gene ontology (GO) pathway analysis was performed for differential expression between BMMC and CD34 datasets, focusing on the top pathways and their biological implications.

10. **Trajectory Analysis**: A subset of cells was selected for trajectory analysis using Monocle 3 to study the progression of cell states.

11. **Cell-Cell Communication**: Cell-cell communication was analyzed using CellChat for both BMMC and CD34 samples, focusing on common signaling pathways.

## R Script:
The R script used for analysis is named **scRNA_seq_script**. 

## Results:
The resulting plots and images can be found in the `results/img` directory of this repository.
