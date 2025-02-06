# tRF_target_prediction
Prediction of 22nt 3' tRNA fragment targets in the mouse transcriptome, with a focus on target sites derived from the primer binding site of LTR-retrotransposons. 

Scripts for transcriptome assembly are heavily adapted from Modzelewski et al 2021. Originals are available from https://epigenome.wustl.edu/TE_Transcript_Assembly/tool.html. 

Except for one-off tasks such as index generation, individual scripts should be executed with the corresponding wrapper under "wrappers". Editing may be necessary for compatibility with your HPC.

Transcriptome assembly

1. If necessary, download public RNA-seq data with "wrappers/download_SRA_wrapper.sh". Metadata or URL inputs can be downloaded from https://sra-explorer.info.
2. Remove adapaters, polyA tails, and low quality bases with "wrappers/trim_wrapper.sh". For SMART-seq+5' data from Oomen et al. 2025, use "wrappers/Oomen_trim_wrapper.sh" followed by "wrappers/trim_wrapper.sh" with adapter trimming set to "F".
3. First pass alignment with "wrappers/STAR_wrapper.sh". If necessary, first generate an index with "transcriptome_assembly/generate_STAR_index.sh".
4. Update the STAR index with new splice junctions with "transcriptome_assembly/generate_STAR_index_second_pass.sh". This requires gathering all SJ.out.tab files from the first alignment into a single directory.
5. Second pass alignment with "wrappers/STAR_wrapper.sh", using the updated STAR index from step 4.
6. Stringtie assembly with "wrappers/stringtie_wrapper.sh". Strandedness must be determined manually (e.g. with "check_strandedness", https://github.com/signalbash/how_are_we_stranded_here)
7. Merge stringtie assembly into a single GTF file and filter for strand information with "transcriptome_assembly/stringtie_merge.sh"). 

ORF prediction

tRF target prediction
