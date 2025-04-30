#plot


library(Seurat)
library(ggplot2)
library(patchwork)


#some plots before filtering

#FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


#VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI", "ribosomalRatio"), pt.size = 0.1, ncol = 3)

#FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle("Genes vs UMIs")

#FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle("UMIs vs % Mito")

#FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "log10GenesPerUMI") + ggtitle("UMIs vs Complexity")

df <- seurat_obj@meta.data



df <- allsample@meta.data


ggplot(df, aes(nFeature_RNA)) +
  geom_histogram(bins = 100, fill = "steelblue") +
  ggtitle("Genes per Cell")

ggplot(df, aes(nCount_RNA)) +
  geom_histogram(bins = 100, fill = "darkgreen") +
  theme_minimal() +
  ggtitle("UMIs per Cell")

ggplot(df, aes(percent.mt)) +
  geom_histogram(bins = 50, fill = "firebrick") +
  theme_minimal() +
  ggtitle("Mitochondrial Ratio (%)")

#basic qc plots

metadata %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")


VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0.001) +
  ggtitle("Total UMI Counts per Cell (Count Depth)")


FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("Gene Count vs UMI Count")

#1-histogram-count depth


ggplot(seurat_obj@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Histogram of Count Depth (nCount_RNA)",
       x = "Total UMIs per Cell", y = "Number of Cells") +
  scale_x_log10()  #  log scale

#2-violin

VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
                                 "log10GenesPerUMI", "ribosomalRatio"),
        ncol = 3, pt.size = 0.1) +
  theme(
    axis.text = element_text(size = 12),        # Axis tick labels
    axis.title = element_text(size = 14),       # Axis titles
    plot.title = element_text(size = 16),       # Plot title
    legend.text = element_text(size = 12),      # Legend labels
    legend.title = element_text(size = 13)      # Legend title
  )


#3. Scatter plots relationships 

FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("Gene vs UMI Count")

FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("Mitochondrial Content vs UMI")

FeatureScatter(seurat_obj, feature1 = "log10GenesPerUMI", feature2 = "nFeature_RNA") +
  ggtitle("Complexity vs Gene Count")

#adjust the sizes of plot details

FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("Mitochondrial Content vs UMI") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )


RidgePlot(seurat_obj_f4, assay = "HTO", features = rownames(seurat_obj_f4[["HTO"]]))
HTOHeatmap(seurat_obj_f4_demuxed, assay = "HTO", classification = "HTO_classification.global")

#VlnPlot(seurat_obj_f4, features = "nCount_HTO", assay = "HTO", pt.size = 0.1)

#summary(seurat_obj_f4$nCount_HTO)
VlnPlot(seurat_obj_f4_demuxed, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "log10GenesPerUMI"),         pt.size = 0.02, ncol = 3)



