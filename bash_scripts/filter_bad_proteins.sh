#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 3:00
#SBATCH -J filter_bad_proteins.sh
#SBATCH -o filter_bad_proteins.log

module load bioinfo-tools augustus/3.3.3-CGP Biopython/1.73-foss-2019a

ANNOT_DIRS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation

for ANNOT_DIR in $(echo $ANNOT_DIRS/*) # add species names to only do select species, e.g D_pond
do 
    cd $ANNOT_DIR
    ANNOT_PROTEINS=isoform_filtered_proteins.faa
    FILTERED_PROTEINS="${ANNOT_PROTEINS%_proteins.*}_no_bad_proteins.faa"
    
    SPECIES="${ANNOT_DIR#$ANNOT_DIRS/}"
    echo $SPECIES
    
    python3 src/filter_bad_proteins.py $ANNOT_PROTEINS $FILTERED_PROTEINS

    # add the species name to the protein fasta headers to be able to tell the difference in orthofinder later
    sed -i "s/>/>${SPECIES}_/g" $FILTERED_PROTEINS

done
 