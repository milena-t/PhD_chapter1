#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 1:00:00
#SBATCH -J blast_repeat_library
#SBATCH -o blast_repeat_library.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools blast/2.13.0+

LIBRARIES_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/repeatmasking/all_repeatmasked/

A_obtectus_REPEATS=${LIBRARIES_DIR}A_obtectus_repeats-families.fa
A_verrucosus_REPEATS=${LIBRARIES_DIR}A_verrucosus_repeats-families.fa
B_siliquastri_REPEATS=${LIBRARIES_DIR}B_siliquastri_repeats-families.fa
C_analis_REPEATS=${LIBRARIES_DIR}C_analis_repeats-families.fa
C_chinensis_REPEATS=${LIBRARIES_DIR}C_chinensis_repeats-families.fa
C_maculatus_REPEATS=${LIBRARIES_DIR}C_maculatus_repeats-families.fa
C_maculatus_superscaffold_REPEATS=${LIBRARIES_DIR}C_maculatus_superscaffold_repeats-families.fa
C_septempunctata_REPEATS=${LIBRARIES_DIR}C_septempunctata_repeats-families.fa
D_melanogaster_REPEATS=${LIBRARIES_DIR}D_melanogaster_repeats-families.fa
D_ponderosae_REPEATS=${LIBRARIES_DIR}D_ponderosae_repeats-families.fa
I_luminosus_REPEATS=${LIBRARIES_DIR}I_luminosus_repeats-families.fa
P_pyralis_REPEATS=${LIBRARIES_DIR}P_pyralis_repeats-families.fa
T_castaneum_REPEATS=${LIBRARIES_DIR}T_castaneum_repeats-families.fa
T_molitor_REPEATS=${LIBRARIES_DIR}T_molitor_repeats-families.fa
Z_morio_REPEATS=${LIBRARIES_DIR}Z_morio_repeats-families.fa

A_obtectus_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/A_obtectus/isoform_filtered_no_bad_proteins.faa
A_verrucosus_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/A_verrucosus/isoform_filtered_no_bad_proteins.faa
B_siliquastri_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/B_siliquastri/isoform_filtered_no_bad_proteins.faa
C_analis_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_analis/isoform_filtered_no_bad_proteins.faa
C_chinensis_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_chinensis/isoform_filtered_no_bad_proteins.faa
C_maculatus_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus/isoform_filtered_no_bad_proteins.faa
C_septempunctata_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_septempunctata/isoform_filtered_no_bad_proteins.faa
D_melanogaster_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/D_melanogaster/isoform_filtered_no_bad_proteins.faa
D_ponderosae_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/D_ponderosae/isoform_filtered_no_bad_proteins.faa
I_luminosus_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/I_luminosus/isoform_filtered_no_bad_proteins.faa
P_pyralis_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/P_pyralis/isoform_filtered_no_bad_proteins.faa
R_ferrugineus_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/R_ferrugineus/isoform_filtered_no_bad_proteins.faa
T_castaneum_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/T_castaneum/isoform_filtered_no_bad_proteins.faa
T_molitor_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/T_molitor/isoform_filtered_no_bad_proteins.faa
Z_morio_PROTEINS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/Z_morio/isoform_filtered_no_bad_proteins.faa

# blastx :  compares a protein query sequence against a nucleotide sequence database dynamically translated in all six reading frames (both strands).

makeblastdb -in $A_obtectus_PROTEINS -dbtype prot
makeblastdb -in $A_verrucosus_PROTEINS -dbtype prot
makeblastdb -in $B_siliquastri_PROTEINS -dbtype prot
makeblastdb -in $C_analis_PROTEINS -dbtype prot
makeblastdb -in $C_chinensis_PROTEINS -dbtype prot
makeblastdb -in $C_maculatus_PROTEINS -dbtype prot
makeblastdb -in $C_septempunctata_PROTEINS -dbtype prot
makeblastdb -in $D_melanogaster_PROTEINS -dbtype prot
makeblastdb -in $D_ponderosae_PROTEINS -dbtype prot
makeblastdb -in $I_luminosus_PROTEINS -dbtype prot
makeblastdb -in $P_pyralis_PROTEINS -dbtype prot
makeblastdb -in $R_ferrugineus_PROTEINS -dbtype prot
makeblastdb -in $T_castaneum_PROTEINS -dbtype prot
makeblastdb -in $T_molitor_PROTEINS -dbtype prot
makeblastdb -in $Z_morio_PROTEINS -dbtype prot

# I will use the evalue filtering here immediately out of blast and then do the bit score filtering in the post processing
# 

blastx -query $A_obtectus_REPEATS -db $A_obtectus_PROTEINS -out A_obtectus_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $A_verrucosus_REPEATS -db $A_verrucosus_PROTEINS -out A_verrucosus_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $B_siliquastri_REPEATS -db $B_siliquastri_PROTEINS -out B_siliquastri_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $C_analis_REPEATS -db $C_analis_PROTEINS -out C_analis_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $C_chinensis_REPEATS -db $C_chinensis_PROTEINS -out C_chinensis_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $C_maculatus_REPEATS -db $C_maculatus_PROTEINS -out C_maculatus_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $C_maculatus_superscaffold_REPEATS -db $C_septempunctata_PROTEINS -out C_septempunctata_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $C_septempunctata_REPEATS -db $D_melanogaster_PROTEINS -out D_melanogaster_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $D_melanogaster_REPEATS -db $D_ponderosae_PROTEINS -out D_ponderosae_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $D_ponderosae_REPEATS -db $I_luminosus_PROTEINS -out I_luminosus_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $I_luminosus_REPEATS -db $P_pyralis_PROTEINS -out P_pyralis_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $P_pyralis_REPEATS -db $R_ferrugineus_PROTEINS -out R_ferrugineus_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $T_castaneum_REPEATS -db $T_castaneum_PROTEINS -out T_castaneum_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $T_molitor_REPEATS -db $T_molitor_PROTEINS -out T_molitor_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10
blastx -query $Z_morio_REPEATS -db $Z_morio_PROTEINS -out Z_morio_blast_repeat_family.out -outfmt 6 -num_threads 20 -evalue 1e-10

