#!/bin/sh
#SBATCH -A naiss2024-5-135
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 20:00
#SBATCH -J single_exon_stats
#SBATCH -o single_exon_stats.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools python3

ANN_PATH=/Users/miltr339/work
SRC_PATH=${ANN_PATH}/PhD_code/PhD_chapter1/src/

# all orthoDB annotations

python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/A_obtectus_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_A_obtectus_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/A_verrucosus_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_A_verrucousus_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/B_siliquastri_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_B_siliquastri_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/C_analis_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_C_analis_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/C_chinensis_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_C_chinensis_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/C_maculatus_superscaffolded_annotation_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_C_maculatus_superscaffolded_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/C_septempunctata_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_C_septempunctata_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/D_ponderosae_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_D_ponderosae_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/D_melanogaster_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_D_melanogaster_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/I_luminosus_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_I_luminosus_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/P_pyralis_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_P_pyralis_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/R_ferrugineus_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_R_ferrugineus_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/T_molitor_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_T_molitor_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/T_castaneum_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_T_castaneum_braker_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/Z_morio_braker_isoform_filtered.gff  True > ${ANN_PATH}/single_exon_genes/orthoDB/orthoDB_Z_morio_braker_single_exon_stats_with_transcript_list.txt

# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/A_obtectus_isoform_filtered.gff  False > orthoDB_acanthoscelides_obtectus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/A_verrucousus_isoform_filtered.gff  False > orthoDB_asbolus_verrucosus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/B_siliquastri_isoform_filtered.gff  False > orthoDB_bruchidius_siliquastri_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/C_analis_isoform_filtered.gff  False > orthoDB_callosobruchus_analis_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/C_chinensis_isoform_filtered.gff  False > orthoDB_callosobruchus_chinensis_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/C_maculatus_isoform_filtered.gff  False > orthoDB_callosobruchus_maculatus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/C_septempunctata_isoform_filtered.gff  False > orthoDB_coccinella_septempunctata_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/D_ponderosae_isoform_filtered.gff  False > orthoDB_dendroctonus_ponderosae_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/D_melanogaster_isoform_filtered.gff  False > orthoDB_drosophila_melanogaster_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/I_luminosus_isoform_filtered.gff  False > orthoDB_ignelater_luminosus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/P_pyralis_isoform_filtered.gff  False > orthoDB_photinus_pyralis_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/R_ferrugineus_isoform_filtered.gff  False > orthoDB_rhynchophorus_ferrugineus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/T_molitor_isoform_filtered.gff  False > orthoDB_tenebrio_molitor_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/T_castaneum_isoform_filtered.gff  False > orthoDB_tribolium_castaneum_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/orthoDB_annotations/Z_morio_isoform_filtered.gff  False > orthoDB_zophobas_morio_single_exon_stats_no_transcript_list.txt


# do all native annotations


echo "---> A_obtectus"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/A_obtectus_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_A_obtectus_single_exon_stats_with_transcript_list.txt
echo "---> A_verrucosus"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/A_verrucosus_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_A_verrucosus_single_exon_stats_with_transcript_list.txt
echo "---> B_siliquastri"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/B_siliquastri_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_B_siliquastri_single_exon_stats_with_transcript_list.txt
echo "---> C_analis"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/C_analis_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_C_analis_single_exon_stats_with_transcript_list.txt
echo "---> C_chinensis"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/C_chinensis_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_C_chinensis_single_exon_stats_with_transcript_list.txt
echo "---> C_maculatus_superscaffolded"
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/C_maculatus_superscaffolded_liftover_annotation.gff True > ${ANN_PATH}/single_exon_genes/native/native_C_maculatus_supersc_single_exon_stats_with_transcript_list.txt
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/liftover_annotation_cmac/C_mac_Y_s_annotation.gff True > ${ANN_PATH}/single_exon_genes/native/native_C_maculatus_supersc_single_exon_stats_with_transcript_list.txt
echo "---> C_septempunctata"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/C_septempunctata_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_C_septempunctata_single_exon_stats_with_transcript_list.txt
echo "---> D_melanogaster"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/D_melanogaster_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_D_melanogaster_single_exon_stats_with_transcript_list.txt
echo "---> D_ponderosae"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/D_ponderosae_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_D_ponderosae_single_exon_stats_with_transcript_list.txt
echo "---> I_luminosus"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/I_luminosus_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_I_luminosus_single_exon_stats_with_transcript_list.txt
echo "---> P_pyralis"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/P_pyralis_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_P_pyralis_single_exon_stats_with_transcript_list.txt
echo "---> R_ferrugineus"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/R_ferrugineus_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_R_ferrugineus_single_exon_stats_with_transcript_list.txt
echo "---> T_castaneum"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/T_castaneum_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_T_castaneum_single_exon_stats_with_transcript_list.txt
echo "---> T_molitor"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/T_molitor_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_T_molitor_single_exon_stats_with_transcript_list.txt
echo "---> Z_morio"
python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/Z_morio_annotation_isoform_filtered.gff True > ${ANN_PATH}/single_exon_genes/native/native_Z_morio_single_exon_stats_with_transcript_list.txt

# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/acanthoscelides_obtectus.gff.isoform_filtered False > native_acanthoscelides_obtectus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/asbolus_verrucosus.gff.isoform_filtered False > native_asbolus_verrucosus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/bruchidius_siliquastri.gff.isoform_filtered False > native_bruchidius_siliquastri_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/callosobruchus_analis.gff.isoform_filtered False > native_callosobruchus_analis_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/callosobruchus_chinensis.gff.isoform_filtered False > native_callosobruchus_chinensis_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/callosobruchus_maculatus.gff.isoform_filtered False > native_callosobruchus_maculatus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/coccinella_septempunctata.gff.isoform_filtered False > native_coccinella_septempunctata_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/dendroctonus_ponderosae.gff.isoform_filtered False > native_dendroctonus_ponderosae_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/drosophila_melanogaster.gff.isoform_filtered False > native_drosophila_melanogaster_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/ignelater_luminosus.gff.isoform_filtered False > native_ignelater_luminosus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/photinus_pyralis.gff.isoform_filtered False > native_photinus_pyralis_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/rhynchophorus_ferrugineus.gff.isoform_filtered False > native_rhynchophorus_ferrugineus_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/tenebrio_molitor.gff.isoform_filtered False > native_tenebrio_molitor_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/tribolium_castaneum.gff.isoform_filtered False > native_tribolium_castaneum_single_exon_stats_no_transcript_list.txt
# python3 ${SRC_PATH}calculate_single_exon_stats.py ${ANN_PATH}/native_annotations/all_native_annot/zophobas_morio.gff.isoform_filtered False > native_zophobas_morio_single_exon_stats_no_transcript_list.txt










