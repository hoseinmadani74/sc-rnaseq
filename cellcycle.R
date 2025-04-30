#cell cycle maus
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Now install AnnotationHub and ensembldb
BiocManager::install(c("AnnotationHub", "ensembldb"))

library(rCurl)
library(AnnotationHub)
library(ensembldb)


cc_file <- "C:/Users/hosei/Desktop/thesis/qc/mus.csv" 
cell_cycle_genes <- read.csv(cc_file)

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))
# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

#data after qc, demuxed
#filtered_seurat <- allsample

# Create new merged RNA assay 
counts_matrix <- GetAssayData(seurat_phase, assay = "RNA", layer = "counts.1")
seurat_phase[["RNA"]] <- CreateAssayObject(counts = counts_matrix)

# Normalize again
seurat_phase <- NormalizeData(seurat_phase, verbose = FALSE)


# Trim the lists?
g2m_genes <- g2m_genes[g2m_genes %in% rownames(seurat_phase)]
s_genes <- s_genes[s_genes %in% rownames(seurat_phase)]

# Perform cell cycle scoring
seurat_phase <- CellCycleScoring(
  seurat_phase,
  g2m.features = g2m_genes,
  s.features = s_genes,
  assay = "RNA"
)


View(seurat_phase@meta.data) 


# to know about the most variables better view for future analysis
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA and color by cell cycle phase
seurat_phase <- RunPCA(seurat_phase)

# Visualize the PCA, grouping by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")