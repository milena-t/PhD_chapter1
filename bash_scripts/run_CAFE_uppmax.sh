#!/bin/bash -l

#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH -J CAFE_%j
#SBATCH -o CAFE_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

CAFE_path=/proj/naiss2023-6-65/Milena/CAFE/CAFE5-5.0/bin/cafe5

ORTHODB_INFILE=/Users/milena/work/chapter1_final_postprocessing/orthofinder_uniform/CAFE_input_orthoDB_from_N0.tsv
TREE_ORTHODB=/Users/milena/work/chapter1_final_postprocessing/orthofinder_uniform/SpeciesTree_uniform.nw
TREE_ODB_MOD="${TREE_ORTHODB%.*}_only_species_names.nw"
sed 's/_filtered_proteinfasta_TE_filtered//g' $TREE_ORTHODB > $TREE_ODB_MOD
CAFE_path -i $ORTHODB_INFILE -t $TREE_ORTHODB -c 2 -o /Users/milena/work/chapter1_final_postprocessing/CAFE_uniform > CAFE_out_uniform.log &