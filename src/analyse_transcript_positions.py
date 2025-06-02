"""
Analyse the relative positions of all transcripts in an orthogroup to try and infer the transcript mechanism.
Assumed that all paralogs that are in close proximity might be tandem duplicated, and dispersed paralogs are
duplicated through different mechanisms
"""

import parse_gff as gff
import parse_orthogroups as OGs


def filepaths_native():
    native_annot_dir = "/Users/miltr339/work/native_annotations/all_native_annot/"
    native_annotations = {
        "A_obtectus" : f"{native_annot_dir}A_obtectus_annotation_isoform_filtered.gff",
        "A_verrucosus" : f"{native_annot_dir}A_verrucosus_annotation_isoform_filtered.gff",
        "B_siliquastri" : f"{native_annot_dir}B_siliquastri_annotation_isoform_filtered.gff",
        "C_analis" : f"{native_annot_dir}C_analis_annotation_isoform_filtered.gff",
        "C_chinensis" : f"{native_annot_dir}C_chinensis_annotation_isoform_filtered.gff",
        "C_maculatus" : f"{native_annot_dir}C_maculatus_superscaffolded_liftover_annotation.gff",
        "C_septempunctata" : f"{native_annot_dir}C_septempunctata_annotation_isoform_filtered.gff",
        "D_melanogaster" : f"{native_annot_dir}D_melanogaster_annotation_isoform_filtered.gff",
        "D_ponderosae" : f"{native_annot_dir}D_ponderosae_annotation_isoform_filtered.gff",
        "I_luminosus" : f"{native_annot_dir}I_luminosus_annotation_isoform_filtered.gff",
        "P_pyralis" : f"{native_annot_dir}P_pyralis_annotation_isoform_filtered.gff",
        "R_ferrugineus" : f"{native_annot_dir}R_ferrugineus_annotation_isoform_filtered.gff",
        "T_castaneum" : f"{native_annot_dir}T_castaneum_annotation_isoform_filtered.gff",
        "T_molitor" : f"{native_annot_dir}T_molitor_annotation_isoform_filtered.gff",
        "Z_morio" : f"{native_annot_dir}Z_morio_annotation_isoform_filtered.gff",
    }
    
    native_proteinseqs_dir = "/Users/miltr339/work/native_proteinseqs/"
    native_proteinseqs={
        "A_obtectus" : f"{native_proteinseqs_dir}A_obtectus.faa",
        "A_verrucosus" : f"{native_proteinseqs_dir}A_verrucosus.faa",
        "B_siliquastri" : f"{native_proteinseqs_dir}B_siliquastri.faa",
        "C_chinensis" : f"{native_proteinseqs_dir}C_chinensis.faa",
        "C_maculatus" : f"{native_proteinseqs_dir}C_maculatus.faa",
        "C_septempunctata" : f"{native_proteinseqs_dir}C_septempunctata.faa",
        "D_melanogaster" : f"{native_proteinseqs_dir}D_melanogaster.faa",
        "D_ponderosae" : f"{native_proteinseqs_dir}D_ponderosae.faa",
        "I_luminosus" : f"{native_proteinseqs_dir}I_luminosus.faa",
        "P_pyralis" : f"{native_proteinseqs_dir}P_pyralis.faa",
        "R_ferrugineus" : f"{native_proteinseqs_dir}R_ferrugineus.faa",
        "T_castaneum" : f"{native_proteinseqs_dir}T_castaneum.faa",
        "T_molitor" : f"{native_proteinseqs_dir}T_molitor.faa",
        "Z_morio" : f"{native_proteinseqs_dir}Z_morio.faa",
    }

    # orthogroups_native = "/Users/miltr339/work/orthofinder/native_orthogroups/N0.tsv"
    # orthogroups_native = "/Users/milena/work/orthofinder/native_orthogroups/N0.tsv"
    # sig_native = "/Users/miltr339/Box Sync/code/CAFE/native_from_N0_Base_family_results.txt"
    orthogroups_native = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/N0.tsv"
    sig_native = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_native_Base_Family_results.txt"

    return native_annotations,orthogroups_native,sig_native,native_proteinseqs

def filepaths_orthoDB():
    orthoDB_annot_dir = "/Users/miltr339/work/orthoDB_annotations/"
    orthoDB_annotations = {
        "A_obtectus" : f"{orthoDB_annot_dir}A_obtectus_braker_isoform_filtered.gff",
        "A_verrucosus" : f"{orthoDB_annot_dir}A_verrucosus_braker_isoform_filtered.gff",
        "B_siliquastri" : f"{orthoDB_annot_dir}B_siliquastri_braker_isoform_filtered.gff",
        "C_analis" : f"{orthoDB_annot_dir}C_analis_braker_isoform_filtered.gff",
        "C_chinensis" : f"{orthoDB_annot_dir}C_chinensis_braker_isoform_filtered.gff",
        "C_maculatus" : f"{orthoDB_annot_dir}C_maculatus_superscaffolded_annotation_isoform_filtered.gff",
        "C_septempunctata" : f"{orthoDB_annot_dir}C_septempunctata_braker_isoform_filtered.gff",
        "D_melanogaster" : f"{orthoDB_annot_dir}D_melanogaster_braker_isoform_filtered.gff",
        "D_ponderosae" : f"{orthoDB_annot_dir}D_ponderosae_braker_isoform_filtered.gff",
        "I_luminosus" : f"{orthoDB_annot_dir}I_luminosus_braker_isoform_filtered.gff",
        "P_pyralis" : f"{orthoDB_annot_dir}P_pyralis_braker_isoform_filtered.gff",
        "R_ferrugineus" : f"{orthoDB_annot_dir}R_ferrugineus_braker_isoform_filtered.gff",
        "T_castaneum" : f"{orthoDB_annot_dir}T_castaneum_braker_isoform_filtered.gff",
        "T_molitor" : f"{orthoDB_annot_dir}T_molitor_braker_isoform_filtered.gff",
        "Z_morio" : f"{orthoDB_annot_dir}Z_morio_braker_isoform_filtered.gff",
    }

    orthoDB_proteinseqs_dir = "/Users/miltr339/work/orthoDB_proteinseqs_TE_filtered/"
    orthoDB_proteinseqs = {
        "A_obtectus" : f"{orthoDB_proteinseqs_dir}A_obtectus_filtered_proteinfasta_TE_filtered.fa",
        "A_verrucosus" : f"{orthoDB_proteinseqs_dir}A_verrucosus_filtered_proteinfasta_TE_filtered.fa",
        "B_siliquastri" : f"{orthoDB_proteinseqs_dir}B_siliquastri_filtered_proteinfasta_TE_filtered.fa",
        "C_analis" : f"{orthoDB_proteinseqs_dir}C_analis_filtered_proteinfasta_TE_filtered.fa",
        "C_chinensis" : f"{orthoDB_proteinseqs_dir}C_chinensis_filtered_proteinfasta_TE_filtered.fa",
        "C_maculatus" : f"{orthoDB_proteinseqs_dir}C_maculatus_filtered_proteinfasta_TE_filtered.fa",
        "C_septempunctata" : f"{orthoDB_proteinseqs_dir}C_septempunctata_filtered_proteinfasta_TE_filtered.fa",
        "D_melanogaster" : f"{orthoDB_proteinseqs_dir}D_melanogaster_filtered_proteinfasta_TE_filtered.fa",
        "D_ponderosae" : f"{orthoDB_proteinseqs_dir}D_ponderosae_filtered_proteinfasta_TE_filtered.fa",
        "I_luminosus" : f"{orthoDB_proteinseqs_dir}I_luminosus_filtered_proteinfasta_TE_filtered.fa",
        "P_pyralis" : f"{orthoDB_proteinseqs_dir}P_pyralis_filtered_proteinfasta_TE_filtered.fa",
        "R_ferrugineus" : f"{orthoDB_proteinseqs_dir}R_ferrugineus_filtered_proteinfasta_TE_filtered.fa",
        "T_castaneum" : f"{orthoDB_proteinseqs_dir}T_castaneum_filtered_proteinfasta_TE_filtered.fa",
        "T_molitor" : f"{orthoDB_proteinseqs_dir}T_molitor_filtered_proteinfasta_TE_filtered.fa",
        "Z_morio" : f"{orthoDB_proteinseqs_dir}Z_morio_filtered_proteinfasta_TE_filtered.fa",
    }

    # orthogroups_orthoDB = "/Users/miltr339/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # orthogroups_orthoDB = "/Users/milena/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # sig_orthoDB = "/Users/miltr339/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"
    
    orthogroups_orthoDB = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv"
    sig_orthoDB = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_uniform_Base_Family_results.txt"

    return orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs





if __name__ == "__main__":
    
    
    orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs = filepaths_orthoDB()
    native_annotations, orthogroups_native, sig_native, native_proteinseqs = filepaths_native()
    dmel_unfiltered_annot = "/Users/miltr339/work/native_annotations/d_melanogaster_NOT_isoform_filtered.gff"    
