#!/bin/bash

## Run ReVis in transcript surroundings mode (https://github.com/milena-t/ReVis)

## TODO check the CAFE5 input file to see what difference it makes!! 
CAFE_RES=/Users/milena/work/PhD_code/ReVis/example_data/CAFE5_Base_family_results.txt
OF_RES=/Users/milena/work/PhD_code/ReVis/example_data/N0.tsv

U_NAME=miltr339

# python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
#     --compute_tables_from_list \
#     --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/B_siliquastri_masking.ori.out \
#     --annotation_gff /Users/${U_NAME}/work/orthoDB_annotations/B_siliquastri_braker_isoform_filtered.gff \
#     --all_list /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/overlap_all_transcripts_B_siliquastri.txt \
#     --sig_list /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/overlap_sig_transcripts_B_siliquastri.txt \
#     --species_name B_siliquastri \
#     --bp 10000 \
#     --GF_size_percentile 90 \
#     --verbose

# python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
#     --compute_tables_from_list \
#     --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
#     --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/C_maculatus_superscaffolded_masking.ori.out \
#     --annotation_gff /Users/${U_NAME}/work/orthoDB_annotations/C_maculatus_superscaffolded_annotation_isoform_filtered.gff \
#     --all_list /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/overlap_all_transcripts_C_maculatus.txt \
#     --sig_list /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/overlap_sig_transcripts_C_maculatus.txt \
#     --species_name C_maculatus \
#     --bp 10000 \
#     --GF_size_percentile 90 \
#     --verbose

python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
    --compute_tables_from_list \
    --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
    --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/A_obtectus_masking.ori.out \
    --annotation_gff /Users/${U_NAME}/work/orthoDB_annotations/A_obtectus_braker_isoform_filtered.gff \
    --all_list /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/overlap_all_transcripts_A_obtectus.txt \
    --sig_list /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/overlap_sig_transcripts_A_obtectus.txt \
    --species_name A_obtectus \
    --bp 10000 \
    --GF_size_percentile 90 \
    --verbose

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