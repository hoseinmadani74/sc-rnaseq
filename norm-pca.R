library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

#we have all samples singlet after qc and demux, merged together

singlets_1$sample_id <- "Sample_1"
singlets_2$sample_id <- "Sample_2"
singlets_3$sample_id <- "Sample_3"
singlets_4$sample_id <- "Sample_4"
singlets_5$sample_id <- "Sample_5"

allsample <- merge(singlets_1, y = list(singlets_2, singlets_3, singlets_4, singlets_5), add.cell.ids = c("S1", "S2", "S3", "S4", "S5"))

# Normalize the counts
seurat_phase <- NormalizeData(allsample)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)


# Identify the ? most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:5]

# Plot the average expression and variance of these genes

#ensure

seurat_phase <- FindVariableFeatures(seurat_phase)  # ensure it's run
seurat_phase@assays$RNA@meta.features <- seurat_phase@assays$RNA@meta.features %>%
  dplyr::filter(variance > 0)

#safer? no


p <- VariableFeaturePlot(seurat_phase)

top_genes <- head(VariableFeatures(seurat_phase), 10)  # or your own gene list

# Only label genes that exist and won’t throw errors
top_genes <- intersect(top_genes, VariableFeatures(seurat_phase))


p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE, xnudge = 0, ynudge = 0)


# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "HTO_maxID",
        split.by = "HTO_maxID")

tail(seurat_phase)
