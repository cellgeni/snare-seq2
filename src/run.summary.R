library(visutils)
library(ArchR)
library(Rsamtools)
library(Seurat)
library(Matrix)

# rmarkdown doesn't work on farm (some issues with Cairo) 
# so let's run Summary from jhub

# make html per-sample summaries
sids = read.table('raw/v04.fq/pipeline.out/samples.txt')
for(s in sids$V1){
  print(s)
  rmarkdown::render('/nfs/cellgeni/pasham/projects/2110.SNARE-seq/src/pipeline/src/summary.Rmd',
                    output_file=paste0('/nfs/users/nfs_p/pm19/nfs/projects/2110.SNARE-seq/raw/v04.fq/pipeline.out/',s,'/summary.html'),
                    params=list(folder = paste0('/nfs/users/nfs_p/pm19/nfs/projects/2110.SNARE-seq/raw/v04.fq/pipeline.out/',s)),
                    clean = TRUE)
}

# try load archR project
p = readRDS('raw/v04.fq/pipeline.out/2_modSHARE_cells/filtered.archr.prj.rds')
p
plotTSSEnrichment(p)
# works!

# agglomerate summary ####
summs = lapply(sids$V1,function(s)readRDS(paste0('raw/v04.fq/pipeline.out/',s,'/qc.rds')))
names(summs) = sids$V1

for(s in names(summs)){
  print(s)
  summs[[s]]$ataca = readRDS(paste0('raw/v04.fq/pipeline.out/',s,'/whole.archr.prj.rds'))
}
names(summs)

summs = summs[c(1,3,4,5,2,6,7,8)]

metrics = do.call(rbind,lapply(summs,function(x)as.data.frame(x$metrics,check.names=FALSE)))

fragsizemax = max(unlist(lapply(summs,function(x)x$fragsize$fragmentPercent)))
tssenrichmax = max(unlist(lapply(summs,function(x)x$tssenrich$smoothValue)))
nFragsmax = max(unlist(lapply(summs,function(v)v$ataca@cellColData$nFrags[v$ataca@cellColData$cdna_cell])))
nFragsmin = min(unlist(lapply(summs,function(v)v$ataca@cellColData$nFrags[v$ataca@cellColData$cdna_cell])))
TSSEnrichmentmax = quantile(unlist(lapply(summs,function(v)v$ataca@cellColData$TSSEnrichment[v$ataca@cellColData$cdna_cell])),p=0.98)

#tssenrichmaxx = max(unlist(lapply(summs,function(x)x$tssenrich$x)))

pdf('raw/v04.fq/pipeline.out/summary.figs.pdf',w=9*3.5,h=8*3)
#png('raw/v04.fq/pipeline.out/summary.figs.png',w=9*3.5,h=8*3,unit='in',res=300)
par(mfrow=c(8,9),bty='n',tcl=-0.2,mgp=c(1.1,0.2,0),mar=c(3,3,1.5,0))
for(s in names(summs)){
  v = summs[[s]]
  
  cells = v$ataca@cellColData$cdna_cell
  plot(v$ataca@cellColData$nFrags+1,v$ataca@cellColData$cdna_total_umi+1,log='xy',pch=16,bty='n',xlab='ATAC fragments',ylab='cDNA UMIs',col=cells+1,cex=0.5,main=s)
  legend('topleft',bty='n',col=1:2,legend=paste0(c('Non-Cells','Cells'),' (',c(sum(!cells),sum(cells)),')'),pch=16)
  
  plot(v$fragsize$fragmentSize,v$fragsize$fragmentPercent,xlab='ATAC-seq Fragment Size (bp)',ylab='Percentage of Fragments',type='l',lwd=2,main=s,ylim=c(0,fragsizemax))
  plot(v$tssenrich$x,v$tssenrich$smoothValue,xlab='Distance From Center (bp)',ylab='Normalized Intensity Profile',type='l',lwd=2,main=s,ylim=c(0,tssenrichmax))
  
  p = v$ataca@cellColData
  p = p[p$cdna_cell,]
  d = pointKde2d(p$nFrags,p$TSSEnrichment,approx = F)
  plot(p$nFrags,p$TSSEnrichment,pch=16,log='x',col=num2col(d),xlab='nFrags',ylab='TSS Enrichment',xlim=c(nFragsmin,nFragsmax),ylim=c(0,TSSEnrichmentmax))
  
  par(mar=c(0,0,1,10))
  plotVisium(v$cdna.umap,v$cdna.umap$clusters,main='cDNA, clusters by cDNA',t='xy',pch=16,show.cluster.sizes = TRUE,label.clusters = TRUE)
  plotVisium(v$cdna.umap,v$atac.umap$clusters,main='cDNA, clusters by ATAC',t='xy',pch=16,show.cluster.sizes = TRUE,label.clusters = TRUE)
  
  plotVisium(v$atac.umap,v$cdna.umap$clusters,main='ATAC, clusters by cDNA',t='xy',pch=16,show.cluster.sizes = TRUE,label.clusters = TRUE)
  plotVisium(v$atac.umap,v$atac.umap$clusters,main='ATAC, clusters by ATAC',t='xy',pch=16,show.cluster.sizes = TRUE,label.clusters = TRUE)
  
  par(mar=c(3,3,1.5,0))
  t = table(v$atac.umap$clusters,v$cdna.umap$clusters)
  imageWithText(t[,order(apply(t,2,which.max))],xlab='clusters by ATAC',ylab='clusters by cDNA')
}
dev.off()


# save metrics
openxlsx::write.xlsx(metrics,'raw/v04.fq/pipeline.out/metrics.xlsx')
write.csv(metrics,'raw/v04.fq/pipeline.out/metrics.csv')

dim(metrics)
pdf('raw/v04.fq/pipeline.out/metrics.pdf',w=6*4.5,h=3*4)
par(mfrow=c(3,6),bty='n',tcl=-0.2,mgp=c(1.1,0.2,0),mar=c(10,3,1.5,0))
for(m in colnames(metrics)[-1:-2]){
  x = as.numeric(gsub('%','',metrics[,m]))
  if(any(grepl('%',metrics[,m])))
    m = paste0(m,' (%)')
  barplot(x,names.arg =  metrics$Sample,las=3,main=m,ylab=m)
}
dev.off()

