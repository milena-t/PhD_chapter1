#!/bin/bash

## Run ReVis in transcript surroundings mode (https://github.com/milena-t/ReVis)

## TODO check the CAFE5 input file to see what difference it makes!! 
CAFE_RES=/Users/milena/work/PhD_code/ReVis/example_data/CAFE5_Base_family_results.txt
OF_RES=/Users/milena/work/PhD_code/ReVis/example_data/N0.tsv

U_NAME=milena

# python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
#     --compute_tables \
#     --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/B_siliquastri_assembly_genomic.fna.out \
#     --annotation_gff /Users/${U_NAME}/work/orthoDB_annotations/B_siliquastri_orthoDB_filtered.gff \
#     --orthogroups $OF_RES \
#     --CAFE5_results $CAFE_RES \
#     --species_name B_siliquastri \
#     --bp 10000 \
#     --GF_size_percentile 90 \
#     --verbose

python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
    --compute_tables \
    --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
    --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/C_mac_superscaffolded_assembly_genomic.fasta.out \
    --annotation_gff /Users/${U_NAME}/work/orthoDB_annotations/C_maculatus_superscaffolded_orthoDB_filtered.gff \
    --orthogroups $OF_RES \
    --CAFE5_results $CAFE_RES \
    --species_name C_maculatus \
    --bp 10000 \
    --GF_size_percentile 90 \
    --verbose

# python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
#     --compute_tables \
#     --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/A_obtectus_assembly_genomic.fna.out \
#     --annotation_gff /Users/${U_NAME}/work/orthoDB_annotations/A_obtectus_orthoDB_filtered.gff \
#     --orthogroups $OF_RES \
#     --CAFE5_results $CAFE_RES \
#     --species_name A_obtectus \
#     --bp 10000 \
#     --GF_size_percentile 90 \
#     --verbose

# python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis.py \
#     --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/C_maculatus_superscaffolded_masking.ori.out \
#     --masker_out_gff /Users/${U_NAME}/work/repeatmasking/repeat_gffs/C_maculatus_superscaffolded_assembly_genomic.fna.out.gff \
#     --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --species_name C_maculatus \
#     --window_length 3e6 \
#     --plot_overlap_filtered \
#     --verbose \
#     --plot

# python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis.py \
#     --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/B_siliquastri_masking.ori.out \
#     --masker_out_gff /Users/${U_NAME}/work/repeatmasking/repeat_gffs/B_siliquastri_assembly_genomic.fna.out.gff \
#     --annotation_gff /Users/${U_NAME}/work/orthoDB_annotations/B_siliquastri_orthoDB_filtered.gff \
#     --merge_gene_windows 2 \
#     --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --species_name B_siliquastri \
#     --window_length 1e6 \
#     --plot_overlap_filtered \
#     --verbose \
#     --plot