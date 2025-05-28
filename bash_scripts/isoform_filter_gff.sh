#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 2:30:00
#SBATCH -J filter_isoforms
#SBATCH -o filter_isoforms_native.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

#filter gff files for the longest isoforms

module load bioinfo-tools AGAT/1.3.2

# ANNOT_GTF=/proj/naiss2023-6-65/Milena/annotation_pipeline/all_proteinrefs_annotation/annotation_species/C_maculatus/braker/braker.gtf
# echo $ANNOT_GTF
# FILTERED_GTF="${ANNOT_GTF%.*}_isoform_filtered.gff" # originally gtf but keep_longest_isoform.pl automatically returns gff version 3 
# 
# perl /proj/naiss2023-6-65/Milena/gene_family_analysis/filter_longest_isoform/agat_sp_keep_longest_isoform.pl -gff $ANNOT_GTF -o $FILTERED_GTF
# echo $(ls -lh $FILTERED_GTF)

W_DIR=$(pwd)

# ANNOT_DIRS=/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Kaufmann_2023_comparison/RNA*
ANNOT_DIRS=/proj/naiss2023-6-65/Milena/gene_family_analysis/native_annotations_gff/*gff

for ANNOT_GTF in $(echo $ANNOT_DIRS)
do 
    echo $ANNOT_GTF
    FILTERED_GTF="${ANNOT_GTF%.*}_isoform_filtered.gff" # originally gtf but keep_longest_isoform.pl automatically returns gff version 3 
    rm $FILTERED_GTF
    OVERLAP_FILTERED_GTF="${ANNOT_GTF%.*}_overlap_filtered.gff" # originally gtf but keep_longest_isoform.pl automatically returns gff version 3 
    rm $OVERLAP_FILTERED_GTF
    perl /proj/naiss2023-6-65/Milena/gene_family_analysis/filter_longest_isoform/agat_sp_fix_overlaping_genes.pl -f $ANNOT_GTF  -o $OVERLAP_FILTERED_GTF
    perl /proj/naiss2023-6-65/Milena/gene_family_analysis/filter_longest_isoform/agat_sp_keep_longest_isoform.pl -gff $OVERLAP_FILTERED_GTF -o $FILTERED_GTF
    rm $OVERLAP_FILTERED_GTF
    echo $(ls -lh $FILTERED_GTF)

done



# extract some filtering stats:

# all filtered isoforms here:
# ls /proj/naiss2023-6-65/Milena/coleoptera_sequences/sequence_downloads/*ncbi_download/ncbi_dataset/data/*/*isoform_filtered* > filtered_isoforms_filepaths.txt

# basic stats from the log files:
# grep "=> N.*: \d*" /proj/naiss2023-6-65/Milena/coleoptera_sequences/sequence_downloads/*ncbi_download/ncbi_dataset/data/*/genomic.agat.log > unfiltered_isoforms_parsing_stats.txt




