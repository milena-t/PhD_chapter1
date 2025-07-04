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
$ORTHODB_ARTHROPODA SRX1541697,SRX1541641,SRX1541696,SRX1541643,SRX1541689,SRX1541646,SRX1541693,SRX1541649,SRX1541651,SRX1541640,SRX1541639

# Lome line (diverse set for differential expression from Kaufmann 2024, same population as assembly)
LOME_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lome_diverse
cd $LOME_DIR
sbatch --job-name="Cmac_scaff_LOME_DE" --output="Cmac_scaff_LOME_diverse.out" -t 5-00:00:00 \
/proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh Cmac_SI_diverse ${LOME_DIR}/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA ERX11759692,ERX11759666,ERX11759673,ERX11759707,ERX11759712,ERX11759732,ERX11759689,ERX11759708,ERX11759670,ERX11759684