#!/bin/bash -l

#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 10:00:00
#SBATCH -J CAFE_uniform
#SBATCH -o CAFE_uniform.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

CAFE_path=/proj/naiss2023-6-65/Milena/CAFE/CAFE5-5.0/bin/cafe5

ORTHODB_INFILE=/proj/naiss2023-6-65/Milena/CAFE/updated_input/CAFE_input_orthoDB_from_N0_filtered.tsv
TREE_ULTRAMETRIC=/proj/naiss2023-6-65/Milena/CAFE/updated_input/native_orthofinder_tree_ultrametrick.nw

# may need to remove the filename suffix of the tree leaf names
# TREE_ODB_MOD="${TREE_ORTHODB%.*}_only_species_names.nw"
# sed 's/_filtered_proteinfasta_TE_filtered//g' $TREE_ORTHODB > $TREE_ODB_MOD

$CAFE_path -i $ORTHODB_INFILE -t $TREE_ULTRAMETRIC -c 16 -o /proj/naiss2023-6-65/Milena/CAFE