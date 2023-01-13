# snare-seq2
The pipeline for processing SNARE-seq cDNA-ATAC multiome. The pipeline is based on STARsolo, BWA and archR.

The entry point is snareseq2.sh.

It is very first version of the pipeline, it depends on enviroment (all tools have to be installed, paths to it and to the indexes are hardcoded in snareseq2.sh).

# TODO
1. Call peaks
2. Implement joint cell calling
3. Count % of mapped reads from ATAC (BWA doesn't report it)
4. Remove temporary files
5. Containerize
6. Post-analysis (?)
7. Change to nextflow (?)
