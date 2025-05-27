#!/bin/bash -l

NATIVE_INFILE=/Users/milena/work/chapter1_final_postprocessing/orthofinder_native/CAFE_input_native_from_N0.tsv
TREE_NATIVE=/Users/milena/work/chapter1_final_postprocessing/orthofinder_native/SpeciesTree_native.nw
TREE_NAT_MOD="${TREE_NATIVE%.*}_only_species_names.nw"
sed 's/_native_isoform_filtered_proteins//g' $TREE_NATIVE > $TREE_NAT_MOD
cafe5 -i $NATIVE_INFILE -t $TREE_NAT_MOD -c 2 -o /Users/milena/work/chapter1_final_postprocessing/CAFE_native > CAFE_out_native.log &

ORTHODB_INFILE=/Users/milena/work/chapter1_final_postprocessing/orthofinder_uniform/CAFE_input_orthoDB_from_N0.tsv
TREE_ORTHODB=/Users/milena/work/chapter1_final_postprocessing/orthofinder_uniform/SpeciesTree_uniform.nw
TREE_ODB_MOD="${TREE_ORTHODB%.*}_only_species_names.nw"
sed 's/_filtered_proteinfasta_TE_filtered//g' $TREE_ORTHODB > $TREE_ODB_MOD
cafe5 -i $ORTHODB_INFILE -t $TREE_ORTHODB -c 2 -o /Users/milena/work/chapter1_final_postprocessing/CAFE_uniform > CAFE_out_uniform.log &