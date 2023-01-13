sids = read.table('raw/v04.fq/sample.names.txt')
sids$V2 = gsub('modSHAREnuclei','modSHARE_nuclei',sids$V2)
sids$type = 'cDNA'
sids$type[grep('AC',sids$V2)] = 'AC'
sids$sample_id = gsub('AC|cDNA|_cDNA|_AC','',sids$V2)
sids = do.call(rbind,lapply(split(sids,sids$sample_id),function(x){
  a = x[x$type=='AC',]
  d = x[x$type=='cDNA',]
  data.frame(sample_id=a$sample_id,
             cDNA_sid=d$V1,
             AC_sid=a$V1)
}))


sids
fqpath='/nfs/users/nfs_p/pm19/nfs/projects/2110.SNARE-seq/raw/v04.fq/'
h=data.frame(TAG=sids$sample_id,
             ATACDIR=paste0(fqpath,sids$AC_sid),
             CDNADIR=paste0(fqpath,sids$cDNA_sid),
             ATACTAG=sids$AC_sid,
             CDNATAG=sids$cDNA_sid,
             SPECIES='human')

fqpath='/lustre/scratch117/cellgen/cellgeni/tmp/Zinah_T319/221205_msBrain_MiSeq_ATAC-RNA/'

msids = c('C3_msBrain_Nu','D2_msBrain_Nu','F2_msBrain_Nu')
m=data.frame(TAG=msids,
             ATACDIR=fqpath,
             CDNADIR=fqpath,
             ATACTAG=paste0('AC-',msids),
             CDNATAG=paste0('cDNA-',msids),
             SPECIES='mouse')

dir.create('raw/v04.fq/pipeline.out')
write.table(rbind(h,m),file = 'raw/v04.fq/pipeline.out/samples.txt',sep = ' ',quote = F,row.names = F,col.names = F)
