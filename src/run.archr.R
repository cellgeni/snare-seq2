#!/software/R-4.1.0/bin/Rscript

library(visutils)
library(ArchR)
library(Rsamtools)
library(Seurat)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
#args = c('/lustre/scratch126/cellgen/cellgeni/pasham/raw/2110.SNARE-seq/v04.fq/pipeline.out/test_C3_msBrain_Nu','mouse','16','/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/barcodes/bc.white.list.txt')
dir = args[1]
species = args[2]
cpu = as.integer(args[3])
WL = readLines(args[4])

setwd(dir)

addArchRThreads(threads = cpu) 
if(species=='human'){
  addArchRGenome("hg38")
}else if(species=='mouse'){
  addArchRGenome("mm10")
}else{
  stop("species should be human or mouse")
}
addArchRChrPrefix(chrPrefix = TRUE)

bam = 'ATAC/out.sorted.bam'
name=basename(getwd())

bamparam = list(isPaired = TRUE,isProperPair = TRUE,isSecondaryAlignment=FALSE,isNotPassingQualityControls=FALSE)
ArrowFiles = createArrowFiles(inputFiles = bam,
                              sampleNames = name,
                              minTSS = 0,
                              minFrags = 0,
                              maxFrags = Inf,
                              addTileMat = TRUE,
                              addGeneScoreMat = TRUE,
                              cleanTmp=T,
                              bcTag='CB',
                              validBarcodes = WL,
                              bamFlag = bamparam,
                              force = TRUE)


p = ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "archR",
  copyArrows = TRUE#This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

# subset by GEX
rna = Seurat::Read10X('cDNA/output/GeneFull/raw/')
rna = colSums(rna)
cells = readLines('cDNA/output/GeneFull/filtered/barcodes.tsv')


p@cellColData$barcode = sapply(strsplit(rownames(p@cellColData),'#'),'[',2)
p@cellColData$cdna_cell = p@cellColData$barcode %in% cells
p@cellColData$cdna_total_umi = rna[p@cellColData$barcode]
saveRDS(p,'whole.archr.prj.rds')

# continue with archer ######
pf = subsetArchRProject(p,rownames(p@cellColData)[p@cellColData$cdna_cell],force = TRUE)


pf <- addIterativeLSI(
  ArchRProj = pf,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

pf <- addClusters(
  input = pf,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1
)

pf <- addUMAP(
  ArchRProj = pf, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

saveRDS(pf,'filtered.archr.prj.rds')

# process GEX
sspipeline = function(d,mdata=NULL,nvargenes=3000,ndim=30){
  d = CreateSeuratObject(d,meta.data = mdata)
  d[["percent.mt"]] <- PercentageFeatureSet(d, pattern = "^MT-")
  d$percent.mt[is.na(d$percent.mt)] = 0
  d = NormalizeData(d)
  
  d = FindVariableFeatures(d, selection.method = "vst", nfeatures = nvargenes)
  
  d <- ScaleData(d, features = rownames(d))
  d <- RunPCA(d, features = VariableFeatures(object = d))
  
  d <- FindNeighbors(d, dims = 1:ndim)
  d <- FindClusters(d, resolution = 0.5)
  d = RunUMAP(d, dims = 1:ndim)
  d
}

g=Seurat::Read10X('cDNA/output/GeneFull/filtered')
g=sspipeline(g)
saveRDS(g,'filtered.cdna.rds')
