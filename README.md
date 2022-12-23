#### * This repo contains the commands used in differential expression analysis of the bulk RNA-SEQ data obtained from the EMT6 nad 4T1 cancer cells co-cultured with the mouse primary DRG neurons.

#### * Before the analysis in R enviroment, Kallisto pseudo-aligner v0.46.1 was installed using debian distribution of the software.

#### * Index of transcriptome (cDNA) file was generated using kallisto index function.

#### * Used mouse references: "Mus_musculus.GRCm39.108.gtf" and "Mus_musculus.GRCm39.cdna.all.fa" (downloaded at 13.11.2022 from ensembl ftp).

#### * Quality check of the fastq files was performed with fastqc and not any trimming was performed before pseudo-alignments.

#### * Tximport and Deseq2 packages were used in the downstream analyzes.
