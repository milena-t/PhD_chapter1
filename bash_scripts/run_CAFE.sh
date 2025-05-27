#!/bin/bash -l

NATIVE_INFILE=/Users/milena/work/chapter1_final_postprocessing/orthofinder_native/CAFE_input_native_from_N0.tsv
TREE_NATIVE=/Users/milena/work/chapter1_final_postprocessing/orthofinder_native/SpeciesTree_native.nw
cafe5 -i $NATIVE_INFILE -t $TREE_NATIVE -c 2 -o /Users/milena/work/chapter1_final_postprocessing/CAFE_native > CAFE_out_native.log &

ORTHODB_INFILE=/Users/milena/work/chapter1_final_postprocessing/orthofinder_uniform/CAFE_input_orthoDB_from_N0.tsv
TREE_ORTHODB=/Users/milena/work/chapter1_final_postprocessing/orthofinder_uniform/SpeciesTree_uniform.nw
cafe5 -i $ORTHODB_INFILE -t $TREE_ORTHODB -c 2 -o /Users/milena/work/chapter1_final_postprocessing/CAFE_uniform > CAFE_out_uniform.log