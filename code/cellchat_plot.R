searchPair(signaling = 'CX3C', pairLR.use = moucellchat$mut@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)

cellchatgene <- function(data1,signal,data2){
  print(signal)
  ligand <- unique(str_match(str_split(netAnalysis_contribution(data1, signaling = signal,return.data =T)$LR.contribution$name,"-",simplify = T)[,1],'[A-Z]*[a-z]*\\w+'))
  receptor <- unique(str_match(str_split(netAnalysis_contribution(data1, signaling = signal,return.data =T)$LR.contribution$name,"-",simplify = T)[,2],'[A-Z]*[a-z]*\\w+'))
  print(paste0('ligad expression'))
  data <- as.data.frame(AverageExpression(data2,assay = 'SCT',ligand,group.by = 'orig.ident'))
  print(data)
  print(data[,1]/data[,2])
  print(paste0('receptor expression'))
  data <- as.data.frame(AverageExpression(data2,assay = 'SCT',receptor,group.by = 'orig.ident'))
  print(data)
  print(data[,1]/data[,2])  
  p2 <-FeaturePlot(data2,c(ligand,receptor),reduction = 'tsne',split.by  = 'orig.ident')
  print(p2)
  print(c(ligand,receptor))
}


library(ArchR) # Analyzing single-cell regulatory chromatin in R. # Analyzing single-cell regulatory chromatin in R. # Analyzing single-cell regulatory chromatin in R. # Analyzing single-cell regulatory chromatin in R.
p.track <- plotBrowserTrack(
  ArchRProj = prog1wt,
  groupBy = "predictedGroup_Un",
  geneSymbol = generxra,
  tileSize = 250,
  upstream = 2E5,
  downstream = 2E5,
  loops = getPeak2GeneLinks(prog1wt)
)

plotPDF(p.track,
        name = paste0("Plot-Tracks-peak2geneWTRXRA.pdf"),
        ArchRProj = proj1,
        addDOC = FALSE, width = 12, height = 6)

p.track <- plotBrowserTrack(
  ArchRProj = prog1mut,
  groupBy = "predictedGroup_Un",
  geneSymbol = generxra,
  tileSize = 250,
  upstream = 2E5,
  downstream = 2E5,
  loops = getPeak2GeneLinks(prog1mut)
)

plotPDF(p.track,
        name = paste0("Plot-Tracks-peak2geneMUTRXRA.pdf"),
        ArchRProj = proj1,
        addDOC = FALSE, width = 12, height = 6)

ArchRBrowser(proj1)


geneok <- stringr::str_to_title(c('kcnc2','cabin1','mapkbp1','chd6'
)) 
p.track <- plotBrowserTrack(
  ArchRProj = prog1wt,
  groupBy = "predictedGroup_Un",
  geneSymbol = geneok,
  tileSize = 250,
  upstream = 2E5,
  downstream = 2E5,
  loops = getPeak2GeneLinks(prog1wt)
)

plotPDF(p.track,
        name = paste0("Plot-Tracks-peak2geneokwt.pdf"),
        ArchRProj = proj1,
        addDOC = FALSE, width = 12, height = 6)

p.track <- plotBrowserTrack(
  ArchRProj = prog1mut,
  groupBy = "predictedGroup_Un",
  geneSymbol = geneok,
  tileSize = 250,
  upstream = 2E5,
  downstream = 2E5,
  loops = getPeak2GeneLinks(prog1mut)
)

plotPDF(p.track,
        name = paste0("Plot-Tracks-peak2geneokmut.pdf"),
        ArchRProj = proj1,
        addDOC = FALSE, width = 12, height = 6)


generemain <- intersect(subset(grnmut.final,Upstream_TF%in%c('Rxra'))$Downstream_Gene,diff_motor$gene)
p.track <- plotBrowserTrack(
  ArchRProj = prog1wt,
  groupBy = "predictedGroup_Un",
  geneSymbol = generemain,
  tileSize = 250,
  upstream = 2E5,
  downstream = 2E5,
  loops = getPeak2GeneLinks(prog1wt)
)

plotPDF(p.track,
        name = paste0("Plot-Tracks-peak2geneokwt.pdf"),
        ArchRProj = proj1,
        addDOC = FALSE, width = 12, height = 6)

p.track <- plotBrowserTrack(
  ArchRProj = prog1mut,
  groupBy = "predictedGroup_Un",
  geneSymbol = generemain,
  tileSize = 250,
  upstream = 2E5,
  downstream = 2E5,
  loops = getPeak2GeneLinks(prog1mut)
)

plotPDF(p.track,
        name = paste0("Plot-Tracks-peak2genemut.pdf"),
        ArchRProj = proj1,
        addDOC = FALSE, width = 12, height = 6)



up.bp <- enrichGO(
  gene  = subset(gene,external_gene_name%in%(subset(grnmut.final,Upstream_TF%in%c('Rxra'))$Downstream_Gene))$ensembl_gene_id,
  keyType = "ENSEMBL",
  OrgDb   = org.Mm.eg.db,
  ont     = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

dotplot(up.bp)

genego <- stri_sort(genego)
p.track <- plotBrowserTrack(
  ArchRProj = prog1wt,
  groupBy = "predictedGroup_Un",
  geneSymbol = genego,
  tileSize = 250,
  upstream = 2E5,
  downstream = 2E5,
  loops = getPeak2GeneLinks(prog1wt)
)

plotPDF(p.track,
        name = paste0("Plot-Tracks-peak2geneokwt.pdf"),
        ArchRProj = proj1,
        addDOC = FALSE, width = 12, height = 6)

p.track <- plotBrowserTrack(
  ArchRProj = prog1mut,
  groupBy = "predictedGroup_Un",
  geneSymbol = genego,
  tileSize = 250,
  upstream = 2E5,
  downstream = 2E5,
  loops = getPeak2GeneLinks(prog1mut)
)

plotPDF(p.track,
        name = paste0("Plot-Tracks-peak2geneokmut.pdf"),
        ArchRProj = proj1,
        addDOC = FALSE, width = 12, height = 6)
