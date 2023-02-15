#!/software/R-4.1.0/bin/Rscript
library(visutils)
library(ArchR)
library(Rsamtools)
library(Seurat)
library(Matrix)
set.seed(1234)
MINUMI = 500
MINFRAGS=1000


args = commandArgs(trailingOnly=TRUE)
# workingdir, human|mouse, nthreads, atacwhitelist, cdnawhitelist
#args = c('/lustre/scratch126/cellgen/cellgeni/pasham/raw/2110.SNARE-seq/v04.fq/pipeline.out/F2_msBrain_Nu','mouse','8','/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/barcodes/bc.white.list.txt','/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/barcodes/bc.white.list.txt')
#args = c('/lustre/scratch126/cellgen/cellgeni/pasham/raw/2110.SNARE-seq/v04.fq/pipeline.out/10x/E7_75_rep1','mouse','16','/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/barcodes/barcode.10x.737K-arc-v1.rev.txt','/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/barcodes/barcode.10x.gex.737K-arc-v1.txt')

dir = args[1]
species = args[2]
cpu = as.integer(args[3])
# alignes atac and cDNA barcodes. they are same for snare, but different for 10x; so it should be ALWAYS used when atac and cdna are analysed jointly
WL = data.frame(atac=readLines(args[4]),
                cdna=readLines(args[5]))

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

bam = 'ATAC/out.sorted_markdupl.bam'
name=basename(getwd())

# set isDuplicate to TRUE to skip deduplication
bamparam = list(isProperPair = TRUE,isSecondaryAlignment=FALSE,isNotPassingQualityControls=FALSE,isDuplicate=FALSE)

ArrowFiles = createArrowFiles(inputFiles = bam,
                              sampleNames = name,
                              minTSS = 0,
                              minFrags = 0,
                              maxFrags = Inf,
                              addTileMat = TRUE,
                              addGeneScoreMat = TRUE,
                              cleanTmp=T,
                              bcTag='CB',
                              validBarcodes = WL$atac,
                              bamFlag = bamparam,
                              force = TRUE)


p = ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "archR",
  copyArrows = TRUE
)

# F2_msBrain_Nu
# sum(p@cellColData$nFrags)
# 2630485 - no markdupl, no flags
# 2630485 - no markdupl, with flags
# 2630485 - markdupl no CB, no flags
#  750575 - markdupl no CB, with flags
#  771542 - markdupl wich CB, with flags

# subset by GEX
rna = Seurat::Read10X('cDNA/output/GeneFull/raw/')
rna = colSums(rna)
cells = readLines('cDNA/output/GeneFull/filtered/barcodes.tsv')


p@cellColData$barcode.atac = sapply(strsplit(rownames(p@cellColData),'#'),'[',2)
p@cellColData$barcode.gex = WL$cdna[match(p@cellColData$barcode.atac,WL$atac)]
p@cellColData$cdna_total_umi = rna[p@cellColData$barcode.gex]
p@cellColData$cdna_total_umi[is.na(p@cellColData$cdna_total_umi)] = 0
p@cellColData$cdna_cell = p@cellColData$barcode.gex %in% cells & p@cellColData$nFrags >= MINFRAGS & p@cellColData$cdna_total_umi >= MINUMI
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

#######################
### peak calling ######
gc()
pf <- addGroupCoverages(ArchRProj = pf, groupBy = "Clusters")
#pathToMacs2 = '/home/jovyan/my-conda-envs/ngs/bin/macs2'
pathToMacs2 <- findMacs2()

pf <- addReproduciblePeakSet(
  ArchRProj = pf, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2
)

pf = addPeakMatrix(pf)
saveRDS(pf,'filtered.archr.prj.rds')
######################


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

cdnaf=Seurat::Read10X('cDNA/output/GeneFull/filtered')
# filter by atac
cdnaf = cdnaf[,pf@cellColData$barcode.gex]

cdnaf=sspipeline(cdnaf)
cdnaf$barcode.atac = WL$atac[match(colnames(cdnaf),WL$cdna)]
saveRDS(cdnaf,'filtered.cdna.rds')
saveRDS(WL,'barcodes.rds')


# save QC #####
fragsize = plotFragmentSizes(ArchRProj = pf,returnDF = TRUE)
tssenrich = plotTSSEnrichment(ArchRProj = pf,returnDF = TRUE)

# prepare UMAPs and clusterings
cdna.umap = as.data.frame(cdnaf@reductions$umap@cell.embeddings)
cdna.umap = cbind(cdna.umap,cdnaf@meta.data)

atac.umap = pf@embeddings$UMAP$df
rownames(atac.umap) = sapply(strsplit(rownames(atac.umap),'#'),'[',2)
ccd = as.data.frame(pf@cellColData)
rownames(ccd) = sapply(strsplit(rownames(ccd),'#'),'[',2)
cmn = intersect(rownames(atac.umap),rownames(ccd)) # sometime smth weird happens here, so it is to ensure correct merge
atac.umap = cbind(atac.umap[cmn,],ccd[cmn,])

wl = WL[WL$atac %in% rownames(atac.umap) & WL$cdna %in% rownames(cdna.umap),]
atac.umap = atac.umap[wl$atac,]
cdna.umap = cdna.umap[wl$cdna,]



jsons = list.files('.',pattern = '*.json',full.names = TRUE)
atacstat = apply(do.call(rbind,lapply(jsons[grepl('atac.+bc',jsons)],function(f){unlist(jsonlite::read_json(f)$general)})),2,sum)
starsolo = read.csv('cDNA/output/GeneFull/Summary.csv',header = F,row.names = 1)
starsolog= read.csv('cDNA/output/Gene/Summary.csv',header = F,row.names = 1)
inxstat = read.table('ATAC/out.sorted_markdupl.idxstats.txt')
cdnaa = Seurat::Read10X('cDNA/output/GeneFull/raw')
peaks = getPeakSet(pf)
dupl.stat = read.table('ATAC/markdupl_metrics.txt',skip = 6,nrows = 1,header = T)

if(length(jsons[grepl('cdna.+bc',jsons)])>0){# SNARE-seq
  cdnastat = apply(do.call(rbind,lapply(jsons[grepl('cdna.+bc',jsons)],function(f){unlist(jsonlite::read_json(f)$general)})),2,sum)
}else{# for 10x I use STAR to fix barcodes
  cdnastat = c()
  cdnastat[1] = starsolo[1,1]
  cdnastat[3] = NA
  cdnastat[3] = round(starsolo[2,1] * starsolo[1,1])
}



info = list()
info$Sample = getwd()
info$Genome = pf@genomeAnnotation$genome
info$nCells = ncol(cdnaf)


cdna.metrics = list()
cdna.metrics[['Total reads']]                    = cdnastat[1]
cdna.metrics[['Valid barcodes no mismatches']]   = cdnastat[2]
cdna.metrics[['Valid barcodes <=1 mismatches']]  = cdnastat[3]
cdna.metrics[['Uniquely mapped reads']]          = starsolo['Reads Mapped to Genome: Unique',1]
cdna.metrics[['Sequencing Saturation']]          = starsolo['Sequencing Saturation',1] 
cdna.metrics[['GeneFull unique reads']]          = starsolo['Reads Mapped to GeneFull: Unique GeneFull',1]
cdna.metrics[['Gene unique reads']]              = starsolog['Reads Mapped to Gene: Unique Gene',1]
cdna.metrics[['Total UMI']]                      = sum(cdnaa)
cdna.metrics[['Fraction MT UMI']]                = sum(cdnaa[grep('^MT-',rownames(cdnaa)),])/sum(cdnaa)
cdna.metrics[['UMI in cells']]                   = sum(cdnaf$nCount_RNA)
cdna.metrics[['Median UMI per cell']]            = median(cdnaf$nCount_RNA)
cdna.metrics[['Median genes per cell']]          = median(cdnaf$nFeature_RNA)
cdna.metrics[['Total genes detected']]           = sum(rowSums(cdnaf@assays$RNA@counts)>0)

atac.metrics = list()
atac.metrics[['Read pairs']]                              = atacstat[1]
atac.metrics[['Valid barcodes no mismatches']]            = atacstat[2]
atac.metrics[['Valid barcodes <=1 mismatches']]           = atacstat[3]
atac.metrics[['Fraction MT reads']]                       = inxstat[inxstat[,1] == 'chrM',3]/sum(inxstat[,3])
atac.metrics[['Fraction non-chr reads']]                  = sum(inxstat[!startsWith(inxstat[,1],'chr'),3])/sum(inxstat[,3])
atac.metrics[['Fraction duplicates']]                     = dupl.stat$PERCENT_DUPLICATION[1]

atac.metrics[['Total fragments used']]                    = sum(p@cellColData$nFrags)
atac.metrics[['Fragments in cells']]                      = sum(pf@cellColData$nFrags)
atac.metrics[['Median fragments per cell']]               = median(pf@cellColData$nFrags)
atac.metrics[['Number of peaks']]                         = length(peaks)
atac.metrics[['Fraction of genome in peaks']]             = sum(width(peaks))/sum(inxstat[startsWith(inxstat[,1],'chr') & inxstat[,1] != 'chrM',2])
atac.metrics[['mean TSS enrichment score']]               = mean(pf@cellColData$TSSEnrichment)
atac.metrics[['Fraction of fragments overlapping TSS']]   = sum(pf@cellColData$ReadsInTSS)/sum(pf@cellColData$nFrags)/2 # seems these are ends, not the fragments; so i divide by 2
atac.metrics[['Fraction of fragments overlapping peaks']] = sum(pf@cellColData$ReadsInPeaks)/sum(pf@cellColData$nFrags)/2

qc = list(cdna.umap=cdna.umap,
          atac.umap=atac.umap,
          metrics=list(info=info,cdna.metrics=cdna.metrics,atac.metrics=atac.metrics),
          fragsize=fragsize,
          tssenrich=tssenrich
)
saveRDS(qc,'qc.rds')
