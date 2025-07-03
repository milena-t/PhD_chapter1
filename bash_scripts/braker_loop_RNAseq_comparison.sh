#!/bin/sh

ORTHODB_ARTHROPODA=/proj/naiss2023-6-65/Milena/annotation_pipeline/annotation_protein_data/OrthoDB_Arthropoda_v11.fa

# # Lu2024
# LU_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lu2024_simple
# cd $LU_DIR
# 
# sbatch --job-name="Cmac_scaff_Lu2024" --output="Cmac_scaff_Lu2024.out" -t 5-00:00:00 \
# /proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh Cmac_Lu2024 ${LU_DIR}/assembly_genomic.fna.masked \
# $ORTHODB_ARTHROPODA SRR27182120,SRR27182119


# South India line (diverse set for Differential expression analysis from Sayadi 2016)
SI_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_SI_diverse
cd $SI_DIR
sbatch --job-name="Cmac_scaff_SI_DE" --output="Cmac_scaff_SI_diverse.out" -t 5-00:00:00 \
/proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh Cmac_SI_diverse ${SI_DIR}/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA ## TODO get right numbers SRR27182120,SRR27182119

# Lome line (diverse set for differential expression from Kaufmann 2024, same population as assembly)
SI_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lome_diverse
cd $SI_DIR
sbatch --job-name="Cmac_scaff_SI_DE" --output="Cmac_scaff_SI_diverse.out" -t 5-00:00:00 \
/proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh Cmac_SI_diverse ${SI_DIR}/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA ## TODO get right numbers SRR27182120,SRR27182119