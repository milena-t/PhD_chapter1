#!/bin/sh

ORTHODB_ARTHROPODA=/proj/naiss2023-6-65/Milena/annotation_pipeline/annotation_protein_data/OrthoDB_Arthropoda_v11.fa
# Lu2024
LU_DIR=/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_Lu2024_simple
cd $LU_DIR

sbatch --job-name="Cmac_scaff_Lu2024" --output="Cmac_scaff_Lu2024.out" -t 5-00:00:00 \
/proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_with_RNAseq_in_SNIC_TMP.sh Cmac_Lu2024 ${LU_DIR}/assembly_genomic.fna.masked \
$ORTHODB_ARTHROPODA SRR27182120,SRR27182119



# cd /proj/naiss2023-6-65/Milena/annotation_pipeline/RNA_proteinrefs_annotation/annotation_species/superscaffolded_comparison/Cmac_SI_diverse
# cd /proj/naiss2023-6-65/Milena/annotation_pipeline/RNA_proteinrefs_annotation/annotation_species/superscaffolded_comparison/Cmac_Lome_diverse