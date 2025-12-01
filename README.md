# tRF target prediction
A pipeline for the prediction of 22nt 3' tRNA fragment (3'-tRF) targets in the mouse genome, with a focus on target sites derived from the primer binding site of LTR-retrotransposons.

To recover target sites in novel LTR-retrotransposon derived transcripts, a transcriptome is assembled from RNA-seq data from the early embryo. Scripts for transcriptome assembly are adapted from Modzelewski et al 2021 (https://epigenome.wustl.edu/TE_Transcript_Assembly/tool.html).

Except for one-off tasks such as index generation and tRF sequence extraction, individual scripts should be executed with the corresponding wrapper under "wrappers". Editing may be necessary for compatibility with your HPC.

## 3'-tRF target prediction

1. miranda/generate_tRF_fasta.R: generate a fasta of tRF3b sequences from mature tRNA sequences downladed from GtRNAdb (https://gtrnadb.ucsc.edu).
2. miranda/create_sliding_windows.sh: divide the input genome into windows of 10,000 bp, overlapping by 50 bp, to accommodate memory limitations of miRanda. 
3. miranda/miranda_genome.sh: run miRanda using either miRNA settings or custom "tRF" settings that unweight the seed.
4. miranda/gather_summaries.sh: gather summary files into one output folder.
5. R_scripts/miranda_to_bed.sh: process summmary files to .csv and .bed output formats.

## Transcriptome assembly

1. If necessary, download public RNA-seq data with "wrappers/download_SRA_wrapper.sh". Metadata or URL inputs can be downloaded from https://sra-explorer.info.
2. Remove adapaters, polyA tails, and low quality bases with "wrappers/trim_wrapper.sh".
3. First pass alignment with "wrappers/STAR_wrapper.sh". If necessary, first generate an index with "transcriptome_assembly/generate_STAR_index.sh".
4. Update the STAR index with new splice junctions with "transcriptome_assembly/generate_STAR_index_second_pass.sh". This requires gathering all SJ.out.tab files from the first alignment into a single directory.
5. Second pass alignment with "wrappers/STAR_wrapper.sh", using the updated STAR index from step 4.
6. Stringtie assembly with "wrappers/stringtie_wrapper.sh". Strandedness must be determined manually (e.g. with "check_strandedness", https://github.com/signalbash/how_are_we_stranded_here)
7. Merge stringtie assembly into a single GTF file and filter for strand information with "transcriptome_assembly/stringtie_merge.sh". This requires moving all .gtf files from step 6 into a single directory. 

## ORF prediction

1. ORF_prediction/find_orfs.sh: find open reading frames in the stringtie-assembled transcriptome.
2. ORF_prediction/blast_orfs.sh: blast predicted open reading frames against the Refseq annotation. Retains only the top scoring hit.

## Annotate target prediction

1. R_scripts/annotate_targets.R: annnotate features of predicted target sites, including their location relative to LTRs and GENCODE/StringTie transcripts.
2. R_scripts/permutation_analysis.sh: calculates the Z-score for overlap between predicted target sites and LTRs at different alignment score thresholds.
