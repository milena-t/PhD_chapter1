#!/bin/bash

## Run ReVis in transcript surroundings mode (https://github.com/milena-t/ReVis)

## TODO check the CAFE5 input file to see what difference it makes!! 

# python3 /Users/miltr339/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
#     --compute_tables \
#     --out_dir /Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --masker_outfile /Users/miltr339/work/repeatmasking/repeat_gffs/B_siliquastri_masking.ori.out \
#     --annotation_gff /Users/miltr339/work/orthoDB_annotations/B_siliquastri_braker_isoform_filtered.gff \
#     --orthogroups /Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv \
#     --CAFE5_results /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence/run1/Base_family_results.txt \
#     --species_name B_siliquastri \
#     --bp 10000 \
#     --GF_size_percentile 90 \
#     --verbose

python3 /Users/miltr339/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
    --compute_tables \
    --out_dir /Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables \
    --masker_outfile /Users/miltr339/work/repeatmasking/repeat_gffs/C_maculatus_repeats.out \
    --annotation_gff /Users/miltr339/work/orthoDB_annotations/C_maculatus_braker_isoform_filtered.gff \
    --orthogroups /Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv \
    --CAFE5_results /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence/run1/Base_family_results.txt \
    --species_name C_maculatus \
    --bp 10000 \
    --GF_size_percentile 90 \
    --verbose

# python3 /Users/miltr339/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
#     --compute_tables \
#     --out_dir /Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --masker_outfile /Users/miltr339/work/repeatmasking/repeat_gffs/A_obtectus_masking.ori.out \
#     --annotation_gff /Users/miltr339/work/orthoDB_annotations/A_obtectus_braker_isoform_filtered.gff \
#     --orthogroups /Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv \
#     --CAFE5_results /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence/run1/Base_family_results.txt \
#     --species_name A_obtectus \
#     --bp 10000 \
#     --GF_size_percentile 90 \
#     --verbose

# python3 /Users/miltr339/work/PhD_code/ReVis/src/ReVis/ReVis.py \
#     --masker_outfile /Users/miltr339/work/repeatmasking/repeat_gffs/C_maculatus_superscaffolded_masking.ori.out \
#     --masker_out_gff /Users/miltr339/work/repeatmasking/repeat_gffs/C_maculatus_superscaffolded_repeats.out.gff \
#     --out_dir /Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --species_name C_maculatus \
#     --window_length 3e6 \
#     --plot_overlap_filtered \
#     --verbose \
#     --plot

# python3 /Users/miltr339/work/PhD_code/ReVis/src/ReVis/ReVis.py \
#     --masker_outfile /Users/miltr339/work/repeatmasking/repeat_gffs/B_siliquastri_masking.ori.out \
#     --masker_out_gff /Users/miltr339/work/repeatmasking/repeat_gffs/B_siliquastri_repeats.out.gff \
#     --annotation_gff /Users/miltr339/work/orthoDB_annotations/B_siliquastri_braker_isoform_filtered.gff \
#     --merge_gene_windows 2 \
#     --out_dir /Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --species_name B_siliquastri \
#     --window_length 1e6 \
#     --plot_overlap_filtered \
#     --verbose \
#     --plot