load('/data2/zj/atac/second.Rdata')

addArchRGenome("mm10")
set.seed(42)

addArchRThreads(threads = 30)

##
prog1wt <- proj1[proj1$Sample%in%('WT-pcdhaq'),]
prog1mut <- proj1[proj1$Sample%in%('MUT-pcdhaq'),]

mouwt <- subset(mou.combined.sct,orig.ident%in%'wt')
mouwt <- CreateSeuratObject(mouwt@assays$RNA@counts,'wt')
mouwt@meta.data <- subset(mou.combined.sct,orig.ident%in%'wt')@meta.data
mouwt$seurat_clusters <- as.factor(as.numeric(mouwt$seurat_clusters))
mouwt <- NormalizeData(mouwt, normalization.method = 'LogNormalize') %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
# mouwt@assays[["SCT"]]@SCTModel.list <- list(model1=mouwt@assays[["SCT"]]@SCTModel.list$model1.1)
prog1wt <- addGeneIntegrationMatrix(
  ArchRProj = prog1wt, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = mouwt,
  addToArrow = T,
  force= TRUE,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

p.hist <- ggplot() + 
  geom_histogram(aes(prog1wt$predictedScore_Un), 
                 fill="grey", color="black", binwidth = 0.05) +
  geom_vline(xintercept = 0.8, color="red") +
  theme_classic()

plotPDF(p.hist, 
        name = "Plot-CellPredictedScorewt.pdf", 
        ArchRProj = prog1wt,
        addDOC = FALSE, width = 6, height = 6)

p.inter <- plotEmbedding(ArchRProj = prog1wt,
                         colorBy = "cellColData", 
                         name = "predictedGroup_Un", 
                         embedding = "UMAP")

plotPDF(p.inter, 
        name = "Plot-UMAP-RNA-Integrationwt.pdf", 
        ArchRProj = prog1wt,
        addDOC = FALSE, width = 5, height = 5)
prog1wt <- addPeak2GeneLinks(
  ArchRProj = prog1wt,
  reducedDims = "IterativeLSI",
  maxDist = 250000,
  predictionCutoff = 0.4
)
peak_count<-getMatrixFromProject(prog1wt,useMatrix = "PeakMatrix")

moumut <- subset(mou.combined.sct,orig.ident%in%'mut')
moumut <- CreateSeuratObject(moumut@assays$RNA@counts,'mut')
moumut@meta.data <- subset(mou.combined.sct,orig.ident%in%'mut')@meta.data
moumut$seurat_clusters <- as.factor(as.numeric(moumut$seurat_clusters))
moumut <- NormalizeData(moumut, normalization.method = 'LogNormalize') %>% FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
# mouwt@assays[["SCT"]]@SCTModel.list <- list(model1=mouwt@assays[["SCT"]]@SCTModel.list$model1.1)
prog1mut <- addGeneIntegrationMatrix(
  ArchRProj = prog1mut, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = mouwt,
  addToArrow = T,
  force= TRUE,
  groupRNA = "seurat_clusters",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

p.hist <- ggplot() + 
  geom_histogram(aes(prog1mut$predictedScore_Un), 
                 fill="grey", color="black", binwidth = 0.05) +
  geom_vline(xintercept = 0.8, color="red") +
  theme_classic()

plotPDF(p.hist, 
        name = "Plot-CellPredictedScoremut.pdf", 
        ArchRProj = prog1mut,
        addDOC = FALSE, width = 6, height = 6)

p.inter <- plotEmbedding(ArchRProj = prog1mut,
                         colorBy = "cellColData", 
                         name = "predictedGroup_Un", 
                         embedding = "UMAP")

plotPDF(p.inter, 
        name = "Plot-UMAP-RNA-Integrationmut.pdf", 
        ArchRProj = prog1wt,
        addDOC = FALSE, width = 5, height = 5)
prog1mut <- addPeak2GeneLinks(
  ArchRProj = prog1mut,
  reducedDims = "IterativeLSI",
  maxDist = 250000,
  predictionCutoff = 0.4
)
peak_count<-getMatrixFromProject(prog1mut,useMatrix = "PeakMatrix")

cell_group <- list(ATAC = split(prog1mut$cellNames, prog1mut$Clusters),
                   RNA = split(names(moumut$seurat_clusters), moumut$seurat_clusters))

load("/data3/index/single_cell/atac_fuction.Rdata")