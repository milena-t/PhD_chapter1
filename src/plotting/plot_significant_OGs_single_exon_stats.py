### some significant transcripts have surroundings associated with LINEs.
### check how many significant transcripts are single-exon vs non-significant transcripts

import parse_gff as gff
import parse_orthogroups as OGs
import plot_basics as my_plotting

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

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
        # "C_analis" : f"{orthoDB_annot_dir}C_analis_braker_isoform_filtered.gff",
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
        # "C_analis" : f"{orthoDB_proteinseqs_dir}C_analis_filtered_proteinfasta_TE_filtered.fa",
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


def get_single_exon_proportion(OGs_dict, annot_dict, species):
    """
    get number of all transcripts and the number of which are single exon for one species
    """
    all_transcripts = 0
    single_exon = 0
    not_found = 0
    for orthogroup, by_species_list in OGs_dict.items():
        transcript_list = by_species_list[species]
        for transcript in transcript_list:
            if len(transcript)==0:
                continue
            transcript=transcript[:-2] # remove tailing "_1"
            all_transcripts += 1
            try:
                annot_transcript = annot_dict[transcript]
                exons = [child_feature for child_feature in annot_transcript.child_ids_list if annot_dict[child_feature].category == gff.FeatureCategory.Exon]
                if len(exons) == 1:
                    single_exon += 1
            except:
                not_found += 1

    if not_found != 0:
        print(f"transcripts not found in annotation: {not_found}")
    return all_transcripts, single_exon


def plot_percentages(perc_sig_dict, perc_all_dict, tree_path, filename = "single_exon_percentages.png"):
    fs = 13 # set font size
    # plot each column in the dataframe as a line in the same plot thorugh a for-loop

    # speciesnames = gff.make_species_order_from_tree(species_tree)
    fig, (ax_data, ax_tree) = plt.subplots(2, 1, figsize=(10, 15), gridspec_kw={'height_ratios': [1, 2]}, constrained_layout=True)
    species_names_unsorted = my_plotting.plot_tree_manually(tree_path, ax_tree)
    # get species order from plotted tree
    species_coords_sorted = sorted(list(species_names_unsorted.keys()))
    speciesnames = [species_names_unsorted[species_coord] for species_coord in species_coords_sorted]
    
    ylab="percentage of single-exon genes"

    perc_all_vec = [perc_all_dict[species] for species in speciesnames]
    perc_sig_vec = [perc_sig_dict[species] for species in speciesnames]
    ax_data.plot(speciesnames, perc_sig_vec, label = f"significant transcripts", color = "#ED7D3A")
    ax_data.plot(speciesnames, perc_all_vec, linestyle=':', label = f"all CAFE transcripts", color = "#ED7D3A")

    ax_data.legend(loc='upper center')
    # set grid only for X axis ticks 
    ax_data.grid(True)
    ax_data.yaxis.grid(False)
    ax_data.set_ylabel(ylab, fontsize = fs)
    ax_data.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x < 0 else f'{int(x)} %'))
    ax_data.set_xticklabels([species.replace("_", ". ") for species in speciesnames], rotation=90, fontsize=fs)

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)



if __name__ == "__main__":
    
    orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs = filepaths_orthoDB()
    native_annotations, orthogroups_native, sig_native, native_proteinseqs = filepaths_native()

    tree = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_only_species_names.nw"

    orthoDB_sig_list, orthoDB_all_list =OGs.get_sig_orthogroups(sig_orthoDB)
    orthoDB_sig_OGs_dict = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_list = orthoDB_sig_list)
    orthoDB_all_OGs_dict = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_list = orthoDB_all_list)

    sig_single_exon_percent = {}
    all_single_exon_percent = {}

    for species in orthoDB_annotations.keys():
        print(f"\n\n\t --> {species}")
        annot_dict = gff.parse_gff3_general(orthoDB_annotations[species], verbose=False)
        
        ## all CAFE genes
        all_transcripts, single_exon = get_single_exon_proportion(orthoDB_all_OGs_dict, annot_dict, species)
        proportion = single_exon/all_transcripts
        percent = proportion*100
        all_single_exon_percent[species] = percent
        print(f"all transcripts single exon proportion: {percent:.2f}% ({single_exon}/{all_transcripts})")

        ## sig CAFE genes
        all_transcripts, single_exon = get_single_exon_proportion(orthoDB_sig_OGs_dict, annot_dict, species)
        proportion = single_exon/all_transcripts
        percent = proportion*100
        sig_single_exon_percent[species] = percent
        print(f"sig transcripts single exon proportion: {percent:.2f}% ({single_exon}/{all_transcripts})")

    print("\n\n")
    for species in sig_single_exon_percent.keys():
        sig_percent = sig_single_exon_percent[species]
        all_percent = all_single_exon_percent[species]
        print(f"{species} single exon proportions: sig. transcripts and all transcripts -->\t {sig_percent:.2f}% // {all_percent:.2f}%")

    data_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data"
    plot_percentages(perc_sig_dict= sig_single_exon_percent, perc_all_dict= all_single_exon_percent, tree_path= tree, filename=f"{data_dir}/orthoDB_single_exon_percentages")

