# tRF target prediction
A pipeline for the prediction of 22nt 3' tRNA fragment targets in the mouse transcriptome, with a focus on target sites derived from the primer binding site of LTR-retrotransposons. To this end, the GENCODE annotation is complemented with transcripts assembled with RNA-seq data from the early embryo, when LTR-retrotransposons are released from silencing and tRF3s are highly abundant. Target prediction is performed with miRanda with seed-weighting removed.

Scripts for transcriptome assembly are adapted from Modzelewski et al 2021 (https://epigenome.wustl.edu/TE_Transcript_Assembly/tool.html). Scripts for SMART-seq+5' read processing are adapted from Oomen et al. 2025 (https://github.com/meoomen/Smartseq5).

Except for one-off tasks such as index generation and tRF sequence extraction, individual scripts should be executed with the corresponding wrapper under "wrappers". Editing may be necessary for compatibility with your HPC.

## Transcriptome assembly

1.1 If necessary, download public RNA-seq data with "wrappers/download_SRA_wrapper.sh". Metadata or URL inputs can be downloaded from https://sra-explorer.info.
1.2 Remove adapaters, polyA tails, and low quality bases with "wrappers/trim_wrapper.sh". For SMART-seq+5' data from Oomen et al. 2025, use "wrappers/Oomen_trim_wrapper.sh" followed by "wrappers/trim_wrapper.sh" with adapter trimming set to "F".
1.3 First pass alignment with "wrappers/STAR_wrapper.sh". If necessary, first generate an index with "transcriptome_assembly/generate_STAR_index.sh".
1.4 Update the STAR index with new splice junctions with "transcriptome_assembly/generate_STAR_index_second_pass.sh". This requires gathering all SJ.out.tab files from the first alignment into a single directory.
1.5 Second pass alignment with "wrappers/STAR_wrapper.sh", using the updated STAR index from step 4.
1.6 Stringtie assembly with "wrappers/stringtie_wrapper.sh". Strandedness must be determined manually (e.g. with "check_strandedness", https://github.com/signalbash/how_are_we_stranded_here)
1.7 Merge stringtie assembly into a single GTF file and filter for strand information with "transcriptome_assembly/stringtie_merge.sh". This requires moving all .gtf files from step 6 into a single directory. 
1.8 Extract fasta sequences from the transcriptome GTF generated in step 7 with "transcriptome_assembly/gtf_to_fasta.sh". Individual transcripts are split into separate fasta files for faster miRanda processing.

## ORF prediction

2.1
2.2

## tRF target prediction

3.1 Generate a fasta of tRF3b sequences from mature tRNA sequenes from GtRNAdb (https://gtrnadb.ucsc.edu) with "miranda/generate_tRF_fasta.R".
3.2 Run miRanda with either default miRNA or custom tRF parameters with "wrappers/miranda_wrapper.sh". Requires transcriptome input from step 1.8 and tRF input from step 3.1. Note this takes a long time given the size of the transcriptome (100,000 transcripts +).

