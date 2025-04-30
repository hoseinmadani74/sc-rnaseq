library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

# Increase allowed memory for future
options(future.globals.maxSize = 2 * 1024^3)  # 2GB 

data.dir = "C:/Users/hosei/Desktop/thesis/no-introns/sample_5/outs/filtered_feature_bc_matrix"
seurat_data <- Read10X(data.dir)
names(seurat_data)
hto_data <- seurat_data[["Antibody Capture"]]  # Modify based on actual structure

seurat_data$`Gene Expression` <- seurat_data$`Gene Expression`[, Matrix::colSums(seurat_data$`Gene Expression`) > 0]
seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`, project = "Sample_5")

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-") # For mouse data
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)

seurat_obj$ribosomalRatio <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl") / 100  # Ribosomal ratio

dim(seurat_obj)

seurat_obj_f4 <- subset(
  seurat_obj,
  subset = 
    nFeature_RNA > 1000 &      # avoid empty droplets
    nCount_RNA > 500 &       # Minimum UMI 
    percent.mt < 15 &         # dying cells
    ribosomalRatio < 0.3 &    
    log10GenesPerUMI > 0.8  # Filters cells with poor complexity
)

dim(seurat_obj_f4)
#common cells for htos to do demux

common_cells <- colnames(seurat_obj_f4)
hto_data <- hto_data[, common_cells]

counts <- GetAssayData(object = seurat_obj_f4, slot = "counts")  

nonzero <- counts > 0

# Keep genes expressed in at least 10 cells
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Subset  only these genes
filtered_counts <- counts[keep_genes, ]

seurat_obj_f4 <- CreateSeuratObject(counts = filtered_counts, meta.data = seurat_obj_f4@meta.data)

hto_data <- as(hto_data, "dgCMatrix")

# Group Mouse1 + Mouse2 into Group1, and so on
grouped_hto_data <- rbind(
  Group1 = Matrix::colMeans(hto_data[c("mouse 1", "mouse 2"), , drop = FALSE]),
  Group2 = Matrix::colMeans(hto_data[c("mouse 3", "mouse 4"), , drop = FALSE]),
  Group3 = Matrix::colMeans(hto_data[c("mouse 5", "mouse 6"), , drop = FALSE])
)


seurat_obj_f4[["HTO"]] <- CreateAssayObject(counts = grouped_hto_data)

rownames(grouped_hto_data)
seurat_obj_f4_norm <- NormalizeData(seurat_obj_f4, assay = "HTO", normalization.method = "CLR")
seurat_obj_f4_demuxed <- HTODemux(seurat_obj_f4_norm, assay = "HTO", positive.quantile = 0.95)

table(seurat_obj_f4_demuxed$HTO_classification.global)

num_singlet <- sum(seurat_obj_f4_demuxed$HTO_classification.global == "Singlet")
total_cells <- length(seurat_obj_f4_demuxed$HTO_classification.global)
percent_singlet <- (num_singlet / total_cells) * 100

print(paste("Percentage of Singlets:", round(percent_singlet, 2), "%"))

# Keep only singlets
singlets_5 <- subset(seurat_obj_f4_demuxed, subset = HTO_classification.global == "Singlet")

singlets_1$sample_id <- "Sample_1"
singlets_2$sample_id <- "Sample_2"
singlets_3$sample_id <- "Sample_3"
singlets_4$sample_id <- "Sample_4"
singlets_5$sample_id <- "Sample_5"

allsample <- merge(singlets_1, y = list(singlets_2, singlets_3, singlets_4, singlets_5), add.cell.ids = c("S1", "S2", "S3", "S4", "S5"))

table(pooled_singlets$sample_id)

pooled_singlets <- allsample

 allsample
# Normalize using SCTransform
options(future.globals.maxSize = 4 * 1024^3)

pooled_singlets <- SCTransform(pooled_singlets, verbose = FALSE)

# PCA, clustering, UMAP
pooled_singlets <- RunPCA(pooled_singlets, verbose = FALSE)
pooled_singlets <- FindNeighbors(pooled_singlets, dims = 1:30)
pooled_singlets <- FindClusters(pooled_singlets, resolution = 0.5)
pooled_singlets <- RunUMAP(pooled_singlets, dims = 1:30)

# Plot grouped by sample origin
DimPlot(pooled_singlets, group.by = "HTO_maxID", reduction = "umap")

