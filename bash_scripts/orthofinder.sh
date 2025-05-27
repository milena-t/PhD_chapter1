#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 2:00:00
#SBATCH -J orthofinder_native
#SBATCH -o orthofinder_native.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools OrthoFinder/2.5.2

# the output directory will be created by orthofinder and should NOT already exist!
# OUT_DIR=/proj/naiss2023-6-65/Milena/gene_family_analysis/orthofinder_orthoDB_overlap_filtered/results_overlap_filtered
OUT_DIR=/proj/naiss2023-6-65/Milena/gene_family_analysis/orthofinder_coleoptera/properly_isoform_filtered
# OUT_DIR=/proj/naiss2023-6-65/Milena/gene_family_analysis/orthofinder_orthoDB_TE_filtered
# orthofinder runs with protein sequences here! nucleotide sequences require the -d flag
# IN_FASTA_DIR=/proj/naiss2023-6-65/Milena/gene_family_analysis/orthofinder_orthoDB_overlap_filtered/overlap_filtered_proteinseqs
IN_FASTA_DIR=/proj/naiss2023-6-65/Milena/gene_family_analysis/proteinseqs_native_annotations/isoform_filtered_for_sure_proteins
# IN_FASTA_DIR=/proj/naiss2023-6-65/Milena/gene_family_analysis/orthoDB_TE_filtered_proteinseqs

# remove it if it does exist, this overwrites all preexisting results
rm -r $OUT_DIR

orthofinder.py -t 20 -a 20 -f $IN_FASTA_DIR -o $OUT_DIR
