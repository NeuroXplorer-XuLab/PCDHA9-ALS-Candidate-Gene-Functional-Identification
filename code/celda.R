sce <- import(sampleDirs = c("/data1/zj_mut/"))
sce.raw <- import(sampleDirs = c("/data1/zj_mut/"), dataType = "raw")
scemut <- decontX(sce, background=sce.raw)
umap <- reducedDim(scemut, "decontX_UMAP")
plotDimReduceCluster(x = scemut$decontX_clusters,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(scemut)
try <- decontXcounts(scemut,val)
plotDecontXResults(scemut, reducedDimName = "decontX_UMAP")
scewt <- readRDS('/data2/zj_sc/rna/scewt.Rdata')
umap <- reducedDim(scewt, "decontX_UMAP")
plotDimReduceCluster(x = scewt$decontX_clusters,
                     dim1 = umap[, 1], dim2 = umap[, 2])
plotDecontXContamination(scewt)
