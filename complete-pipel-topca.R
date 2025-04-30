library(Seurat)
library(ggplot2)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)

seurat_list <- list()

singlet_list <- list()


for (i in 1:5) {
  s <- paste0("sample_", i)
  message("Processing of ", s)
  
  data.dir <- paste0("C:/Users/hosei/Desktop/thesis/no-introns/", s, "/outs/filtered_feature_bc_matrix")
  
  # Read data
  seurat_data <- Read10X(data.dir)
  
  # Adjust for HTO or other features
  hto_data <- seurat_data[["Antibody Capture"]]  # optional, depending on structure
  
  # Filter zero-column cells
  seurat_data$`Gene Expression` <- seurat_data$`Gene Expression`[, Matrix::colSums(seurat_data$`Gene Expression`) > 0]
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = seurat_data$`Gene Expression`, project = s)
  
  # QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")  # Mouse mitochondrial genes
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  seurat_obj$ribosomalRatio <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl") / 100  # Ribosomal
  
  # Optionally save or store for later
  saveRDS(seurat_obj, file = paste0(s, "_seurat_obj.rds"))
  
  
  #df <- seurat_obj@meta.data
  
  seurat_obj_f4 <- subset(
    seurat_obj,
    subset = 
      nFeature_RNA > 500 &      # avoid empty droplets
      nCount_RNA > 500 &       # Minimum UMI 
      percent.mt < 15 &         # dying cells
      ribosomalRatio < 0.3 &    
      log10GenesPerUMI > 0.8  # Filters cells with poor complexity
  )
  
  dim(seurat_obj_f4)
  #common cells for htos to do demu
  
  seurat_list[[s]] <- seurat_obj_f4
  
  df <- seurat_obj_f4@meta.data
  
  #1 ncell

  p0 <- ggplot(df,(aes(x=s, fill=s)))  +
    geom_bar(fill = "black") +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle(s,"Numebr of cells")
  
  print(p0)
  
  ggsave((paste0(s, "_NCells.png")), plot = p0, width = 6, height = 4, dpi = 300)
  ggsave((paste0(s, "_NCells.pdf")), plot = p0, width = 6, height = 4, dpi = 300)
  
  
  #2 nFeature , nGene
  
  p1 <- ggplot(df, aes(nFeature_RNA)) +
    geom_histogram(bins = 100, fill = "purple") +
    ggtitle(s,"Genes per Cell")
  
  print(p1)
  ggsave((paste0(s, "_nGene.png")), plot = p1, width = 6, height = 4, dpi = 300)
  ggsave((paste0(s, "_nGene.pdf")), plot = p1, width = 6, height = 4, dpi = 300)
  
  
  #3 UMIs, ncount
  
  p2 <- df %>%
    ggplot(aes(nCount_RNA)) +
    geom_histogram(bins = 100, fill = "darkgreen") +
    ggtitle(s, "UMIs per Cell")
  
  print(p2)
  ggsave((paste0(s, "_nUMIs.png")), plot = p2, width = 6, height = 4, dpi = 300)
  ggsave((paste0(s, "_nUMIs.pdf")), plot = p2, width = 6, height = 4, dpi = 300)
  
  
  
  #4 complexity
  
  
  p4 <- df %>%
    ggplot(aes(x = log10GenesPerUMI, color = s, fill = s)) +
    geom_histogram(aes(y = ..density..), bins = 100, alpha = 0.3, position = "identity") +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 0.8, linetype = "dashed", color = "darkred") 
  
  print(p4)   
  ggsave((paste0(s, "_Complexity_density.png")), plot = p4, width = 6, height = 4, dpi = 300)
  ggsave((paste0(s, "_Complexity_density.pdf")), plot = p4, width = 6, height = 4, dpi = 300)
  
  
  
  
  #5 Mito
  
  p5 <- ggplot(df, aes(percent.mt)) +
    geom_histogram(bins = 100, fill = "brown") +
    ggtitle(s, "Mitochondrial Ratio (%)")
  
  print(p5)   
  
  ggsave((paste0(s, "_MitoRatio.png")), plot = p5, width = 6, height = 4, dpi = 300)
  ggsave((paste0(s, "_MitoRatio.pdf")), plot = p5, width = 6, height = 4, dpi = 300)
  
  
}


