data <- Read10X("/data1/zj_wt/outs/filtered_feature_bc_matrix/")
wt <- CreateSeuratObject(counts = data, project = "wt", min.cells = 3, min.features = 200)
wt[["percent.mt"]] <- PercentageFeatureSet(object = wt, pattern = "^mt-")
VlnPlot(wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(wt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

wt <- SCTransform(wt, method = "glmGamPoi",verbose = FALSE)
wt <- RunPCA(wt, npcs = 50)

wt <- RunUMAP(wt,
              reduction = "pca",
              dims = 1:pc_threshold)
wt <- FindNeighbors(wt,
                    reduction = "pca",
                    dims = 1:pc_threshold)
wt <- FindClusters(wt,
                   resolution = 0.4,
                   algorithm = 1,
)

sweep.mue14 <- paramSweep_v3(wt, PCs = 1:pc_threshold, sct = TRUE)
sweep.stats_pfc77 <- summarizeSweep(sweep.mue14 , GT = FALSE)
bcmvn <- find.pK(sweep.stats_pfc77)
homotypic.prop <- modelHomotypic(wt@meta.data$seurat_clusters) 
DoubletRate = ncol(wt)*8*1e-6 ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <-  round(DoubletRate*ncol(wt))   ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #提取最佳pk值
# 计算双细胞比例
wt <- findouble(wt)
plot8<-UMAPPlot(wt, group.by="DF_hi.lo",label.color=c("black","gold","red"))
print(plot8)
dev.off()
library(SeuratDisk) # Interfaces for HDF5-Based Single Cell File Formats
SaveH5Seurat(wt, filename = "wt.h5Seurat")
Convert("wt.h5Seurat", dest = "h5ad")
system('rm wt.h5Seurat')


SaveH5Seurat(mut, filename = "mut.h5Seurat")
Convert("mut.h5Seurat", dest = "h5ad")
system('rm mut.h5Seurat')

data <- Read10X("/data1/zj_wt/outs/filtered_feature_bc_matrix/")
scewt <- readRDS('/data2/zj_sc/rna/scewt.Rdata')
cellwt <- rownames(data.frame(data=scewt@metadata$decontX$estimates$all_cells$contamination) %>% filter(data<0.7))
cellwt <-str_split(cellwt,"_",simplify = T)[,3]
data <- data[,cellwt]
wt <- CreateSeuratObject(counts = data, project = "wt", min.cells = 3, min.features = 200)
wt[["percent.mt"]] <- PercentageFeatureSet(object = wt, pattern = "^mt-")
VlnPlot(wt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(wt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(wt, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

wt <- subset(wt, subset = nCount_RNA<60000&nFeature_RNA > 500 & percent.mt <5)


data <- Read10X("/data1/zj_mut/outs/filtered_feature_bc_matrix/")
scemut <- readRDS('/data2/zj_sc/rna/scemut.Rdata')
plotDecontXContamination(scemut)
ggsave(paste(paste0("/data2/zj_sc/rna/",'decoxmut'),"pdf",sep = "."), height = 15, width = 13)
cellmut <- rownames(data.frame(data=scemut@metadata$decontX$estimates$all_cells$contamination) %>% filter(data<0.7))
cellmut <-str_split(cellmut,"_",simplify = T)[,3]
data <- data[,cellmut]
mut <- CreateSeuratObject(counts = data, project = "mut", min.cells = 3, min.features = 200)
mut[["percent.mt"]] <- PercentageFeatureSet(object = mut, pattern = "^mt-")
VlnPlot(mut, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(mut, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mut, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

mut <- subset(mut, subset = nCount_RNA<60000&nFeature_RNA > 500 & percent.mt <5)


##integrate

library(Seurat) # Tools for Single Cell Genomics # Tools for Single Cell Genomics # Tools for Single Cell Genomics
library(SeuratData) # Install and Manage Seurat Datasets
library(patchwork) # The Composer of Plots
library(dplyr) # A Grammar of Data Manipulation # A Grammar of Data Manipulation # A Grammar of Data Manipulation # A Grammar of Data Manipulation
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics

moumerge <- merge(x = wt, y = mut,
                    add.cell.ids = c('wt', 'mut'))
moulist <- SplitObject(moumerge, split.by = "orig.ident")
mou.combined.sct <- IntegrateData_own(anchorset = moulist)
mou.combined.sct <- RunPCA(mou.combined.sct, verbose = FALSE)
ElbowPlot(object = mou.combined.sct, ndims = 30, reduction = "pca") 
mou.combined.sct <- RunUMAP(mou.combined.sct, reduction = "pca", dims = 1:15, verbose = FALSE)
mou.combined.sct <- FindNeighbors(mou.combined.sct, reduction = "pca", dims = 1:15)
mou.combined.sct <- FindClusters(mou.combined.sct,resolution = 0.1, verbose = FALSE,algorithm = 4)
p1 <- DimPlot(mou.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(mou.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
              repel = TRUE)
p1|p2
ggsave(paste(paste0("/data2/zj_sc/rna/",'celltype'),"pdf",sep = "."), height =15, width = 13)

diff2<- FindMarkers(mou.combined.sct2, assay = "SCT", ident.1 = "wt", ident.2 = "mut",
                    group='orig.ident',only.pos = TRUE,verbose = FALSE)

up.bp <- enrichGO(
  gene  = subset(gene,external_gene_name%in%rownames(diff))$ensembl_gene_id,
  keyType = "ENSEMBL",
  OrgDb   = org.Mm.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

dotplot(up.bp)
mou.combined.sct@meta.data[mou.combined.sct$seurat_clusters%in%c(11,6,8,2,5,9,10,12),]$seurat_clusters <- '2'