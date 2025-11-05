#!/bin/bash

## Run ReVis in transcript surroundings mode (https://github.com/milena-t/ReVis)

CAFE_RES=/Users/milena/work/PhD_code/ReVis/example_data/CAFE5_Base_family_results.txt
OF_RES=/Users/milena/work/PhD_code/ReVis/example_data/N0.tsv

U_NAME=miltr339

# for SPECIES in A_obtectus # A_verrucosus C_chinensis C_maculatus C_septempunctata D_melanogaster D_ponderosae I_luminosus P_pyralis R_ferrugineus T_castaneum T_molitor Z_morio # B_siliquastri 
# do
#     echo ""
#     echo ""
#     echo "--------------------------------------------------------------------------------------------------------------"
#     echo "-------------------------------------------${SPECIES}-------------------------------------------------------"
#     echo "--------------------------------------------------------------------------------------------------------------"
# 
#     python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
#         --compute_tables_from_list \
#         --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
#         --masker_outfile /Users/${U_NAME}/work/repeatmasking/repeat_gffs/${SPECIES}_repeats.out \
#         --annotation_gff /Users/${U_NAME}/work/orthoDB_annotations/${SPECIES}_braker_isoform_filtered.gff \
#         --all_list /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/overlap_all_transcripts_${SPECIES}.txt \
#         --sig_list /Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/overlap_sig_transcripts_${SPECIES}.txt \
#         --species_name ${SPECIES} \
#         --bp 10000 \
#         --verbose
# done   

## do the plotting based on existing tables

REP_TABLES=/Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables
SPECIES=C_maculatus
for SPECIES in A_obtectus A_verrucosus B_siliquastri C_chinensis C_maculatus C_septempunctata D_melanogaster D_ponderosae I_luminosus P_pyralis R_ferrugineus T_castaneum T_molitor Z_morio # B_siliquastri 
do
    python3 /Users/${U_NAME}/work/PhD_code/ReVis/src/ReVis/ReVis_transcript_surroundings.py \
        --plot \
        --out_dir /Users/${U_NAME}/work/PhD_code/PhD_chapter1/data/repeats_tables \
        --all_before_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_before_all_transcripts.txt \
        --sig_before_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt \
        --all_after_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_after_all_transcripts.txt \
        --sig_after_table ${REP_TABLES}/${SPECIES}_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt \
        --num_transcripts ${REP_TABLES}/${SPECIES}_transcript_numbers.txt \
        --species_name $SPECIES \
        --verbose
done
##


### I can't do the repeat histograms in a loop because they all require different window lengths
