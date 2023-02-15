# snare-seq2
The pipeline for processing SNARE-seq cDNA-ATAC multiome. The pipeline is based on STARsolo, BWA and archR.

The entry point is snareseq2.sh for SNARE-seq and 10xmulti.sh for 10x multiome.

It is very first version of the pipeline, it depends on enviroment (all tools have to be installed, paths to it and to the indexes are hardcoded in snareseq2.sh).

# TODO
1. Implement joint cell calling
2. Containerize
3. Post-analysis (?)
4. Change to nextflow (?)
