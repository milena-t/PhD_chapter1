# plot the batch summary statistics of busco analysis from native and orthoDB annotations

import re
import parse_gff as gff
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

def make_species_order_from_tree(newick_tree_path):
    # Regular expression to extract leaf names
    # This matches strings between commas, parentheses, and before colons.
    leaf_pattern = r'(?<=\(|,)([a-zA-Z0-9_]+)(?=:)'
    with open(newick_tree_path, "r") as newick_tree_file:
        newick_tree_string = newick_tree_file.readlines()[0]
        # print(newick_tree_string)
        leaf_names = re.findall(leaf_pattern, newick_tree_string)
        # leaf names are like "A_obtectus_filtered_proteinfasta" but we only care about the species names in the beginning
        species_names = [gff.split_at_second_occurrence(leaf, "_") for leaf in leaf_names]
    return species_names



def get_tree_name_from_header(header_str, species_names):
    header_found = False
    # check for full species name 
    for species in species_names:
        if species in header_str:
            header_found = True
            return(species)
    # check only for second part of the species name (obtectus instead of A_obtectus)
    if not header_found:
        for species in species_names:
            species = species.split("_")[-1]
            if species in header_str:
                header_found = True
                return(species)

    if not header_found:
        return(header_str)

def read_busco_summary(busco_batch_summary_path, species_list, name_association = {}):
    # make nested dictionary: { species : { Dataset	Complete : float , Single : float ,  Duplicated : float , Fragmented : float , Missing : float , n_markers : int } }
    busco_stats_dict = {}
    with open(busco_batch_summary_path, "r") as batch_summary:
        header = [heading.strip() for heading in batch_summary.readline().strip().split("\t")] 
        
        for species_line_str in batch_summary:

            line_dict = {}     

            species_line = species_line_str.strip().split("\t")
            if len(name_association) == 0:
                species = get_tree_name_from_header(species_line[0].strip(), species_list)
            else:
                try:
                    species = get_tree_name_from_header(name_association[species_line[0].strip()], species_list)

                except:
                    print(species_line[0])
                    species = species_line[0]

            if "failed" in species_line_str:
                for i in list(range(2,len(header)-1)): # exclude filename and Dataset headings
                    line_dict[header[i]] = 0.0     
            else:
                for i in list(range(2,len(header)-1)): # exclude filename and Dataset headings
                    line_dict[header[i]] = float(species_line[i].strip())

            busco_stats_dict[species] = line_dict
        
        return busco_stats_dict



def plot_busco_stats(native_busco, orthoDB_busco = {}, species_names = [], outfile_name = ""):
    print(f" plotting for these {len(species_names)} species: \n{species_names}")

    colors_orthoDB = {
        "Complete" : "#D36D0D",
        "Single" : "#F2933A" , # this is the same shade as uniform masking yellow in all the other plots
        "Duplicated" : "#F4B352",
        "Fragmented" : "#F8CF8B",
        "Missing" : "#828282"
    }

    colors_native = {
        "Complete" : "#582829",
        "Single" : "#863c3d" , # this is the same shade as native red in all the other plots
        "Duplicated" : "#cd5b5b",
        "Fragmented" : "#ff7277",
        "Missing" : "#828282"
    }

    # X coordinates for the groups
    x = np.arange(len(species_names))
    # Width of the bars
    long_figure = False
    if len(orthoDB_busco) == 0:
        width = 0.8 # (this is a fraction of the standardized 1 unit of space between axis ticks)
    else:
        width = 0.35
        long_figure = True
    # position bars right
    if len(orthoDB_busco)>0:
        annotation_label = " (native)"
        x_subtr = width/2 + width/15
    else:
        annotation_label = ""
        x_subtr = 0


    legend_columns = 2
    fs = 50 # fontsize is scaled with the dpi somehow which i have to do extra because i change the aspect ratio manually below

    # set figure aspect ratio
    if long_figure:
        aspect_ratio = 30 / 12
    else:
        aspect_ratio = 17 / 12

    height_pixels = 2000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig = plt.figure(figsize=(width_pixels / 100, height_pixels / 100), dpi=100, frameon = False)
    ax = fig.add_subplot(111)

    native_numbers = {
        # 'Complete' : np.array([native_busco.get(species, {}).get('Complete', 0) for species in species_names]),
        "Single" : np.array([native_busco.get(species, {}).get("Single", 0) for species in species_names]),
        "Duplicated" : np.array([native_busco.get(species, {}).get("Duplicated", 0) for species in species_names]),
        "Fragmented" : np.array([native_busco.get(species, {}).get("Fragmented", 0) for species in species_names]),
        "Missing" : np.array([native_busco.get(species, {}).get("Missing", 0) for species in species_names])
    }
        
    bottom = np.zeros(len(species_names))
    for busco_category, busco_numbers in native_numbers.items():
        p = ax.bar(x - x_subtr, busco_numbers, width, label=busco_category + annotation_label, color= colors_native[busco_category], bottom=bottom)
        bottom += busco_numbers

    if len(orthoDB_busco)>0:
        annotation_label = " (uniform)"
        legend_columns = 4

        orthoDB_numbers = {
            # 'Complete' : np.array([orthoDB_busco.get(species, {}).get('Complete', 0) for species in species_names]),
            "Single" : np.array([orthoDB_busco.get(species, {}).get("Single", 0) for species in species_names]),
            "Duplicated" : np.array([orthoDB_busco.get(species, {}).get("Duplicated", 0) for species in species_names]),
            "Fragmented" : np.array([orthoDB_busco.get(species, {}).get("Fragmented", 0) for species in species_names]),
            "Missing" : np.array([orthoDB_busco.get(species, {}).get("Missing", 0) for species in species_names])
        }

        missing_species = [species for species in species_names if species not in orthoDB_busco]
        print(f"Missing species in orthoDB: {missing_species}")

        bottom = np.zeros(len(species_names))
        for busco_category, busco_numbers in orthoDB_numbers.items():
            #p = ax.bar(x + x_subtr, busco_numbers, width, label=busco_category, color= colors[busco_category], bottom=bottom)
            try:
                #p = ax.bar(x + x_subtr, busco_numbers, width, label=busco_category, color= colors[busco_category], bottom=bottom)
                p = ax.bar(x + x_subtr, busco_numbers, width, label=busco_category + annotation_label, color= colors_orthoDB[busco_category], bottom=bottom)
                bottom += busco_numbers
            except:
                print(busco_category)
                print(orthoDB_numbers)

                raise RuntimeError

    ### set up labels and axes etc. ###

    # y axis ticks
    ax.set_ylabel('percentage of BUSCOs', fontsize=fs+4)
    # if len(orthoDB_busco)==0:
    #     ax.set_title('busco completeness', fontsize=fs+4)
    # else:
    #     ax.set_title('busco completeness (left: native, right: orthoDB)', fontsize=fs+4)
    ax.set_xticks(x)
    ax.set_xlabel('', fontsize=fs+4)
    xtick_labels = [species.replace("_", ". ") for species in species_names]
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=fs)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 100 else f'{int(x)}%'))
    ax.tick_params(axis='y', labelsize=fs)
    # ax.set_yticklabels([int(tick) for tick in ax.get_yticks() if tick <101], fontsize=fs)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize=fs*0.8, ncol=legend_columns, loc='upper center')
    # add space at the top of the plot for the legend
    ymax = 100*1.2
    ax.set_ylim(0, int(ymax))
    ax.set_xlim(-0.5, len(xtick_labels)-0.5)

    # plt.tight_layout()

    plt.savefig(outfile_name, dpi = 300, transparent = True, bbox_inches='tight')
    print("Figure saved in the current working directory directory as: "+outfile_name)
    





if __name__ == '__main__':

    native_name_association = {
        "Coccinella_septempunctata_ncbi_GCF_907165205.1.fna_filtered_bad_transcripts_filtered.faa" : "C_septempunctata",
        "Drosophila_melanogaster_ncbi_GCF_000001215.4.fna_filtered_bad_transcripts_filtered.faa" : "D_melanogaster",
        "GCA_027724725.1_ASM2772472v1__transcripts_isoform_filtered_bad_transcripts_filtered.faa" : "Z_morio",
        "GCA_011009095.1_Ilumi1.2__transcripts_isoform_filtered_bad_transcripts_filtered.faa" : "I_luminosus",
        "Bruchidius_siliquastri-GCA_949316355.1.fna_filtered_bad_transcripts_filtered.faa" : "B_siliquastri",
        "Callosobruchus_chinensis.fna_filtered_bad_transcripts_filtered.faa" : "C_chinensis",
        "GCA_008802855.1_Ppyr1.3__transcripts_isoform_filtered_bad_transcripts_filtered.faa" : "P_pyralis",
        "GCA_014462685.1_ASM1446268v1__transcripts_isoform_filtered_bad_transcripts_filtered.faa" : "R_ferrugineus",
        "acanthoscelides_obtectus.fna_filtered_bad_transcripts_filtered_unique_headers.faa" : "A_obtectus",
        # "acanthoscelides_obtectus.fna_filtered_bad_transcripts_filtered.faa" : "A_obtectus",
        "Callosobruchus_analis.fna_filtered_bad_transcripts_filtered.faa" : "C_analis",
        "GCA_020466635.2_Dpon_M_20191212v2__transcripts_isoform_filtered_bad_transcripts_filtered.faa" : "D_ponderosae",
        "GCA_004193795.1_BDFB_1.0__transcripts_isoform_filtered_bad_transcripts_filtered.faa" : "A_verrucosus",
        "Tribolium_castaneum_NCBI_GCF_031307605.1.fna_filtered_bad_transcripts_filtered.faa" : "T_castaneum",
        "GCA_027725215.1_ASM2772521v1__transcripts_isoform_filtered_bad_transcripts_filtered.faa" : "T_molitor",
        "callosobruchus_maculatus.fna_filtered_bad_transcripts_filtered.faa" : "C_maculatus"
    }

    species_tree = "/Users/miltr339/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw"
    species_names = make_species_order_from_tree(species_tree)

    orthoDB_batch_uppmax = "/proj/naiss2023-6-65/Milena/annotation_pipeline/busco/busco/busco2/BUSCO_orthoDB_annotation/batch_summary.txt"
    native_batch_uppmax = "/proj/naiss2023-6-65/Milena/annotation_pipeline/busco/busco/busco2/BUSCO_native_annotation/batch_summary.txt"

    orthoDB_batch = "/Users/miltr339/work/busco/orthoDB_batch_summary.txt"    
    native_batch = "/Users/miltr339/work/busco/native_batch_summary.txt"

    orthoDB_busco_stats_dict = read_busco_summary(orthoDB_batch, species_names)
    # print(orthoDB_busco_stats_dict)
    # for key, value in orthoDB_busco_stats_dict.items():
    #     print(key)
    #     print(value)

    native_busco_stats_dict = read_busco_summary(native_batch, species_names, name_association=native_name_association)
    # print(native_busco_stats_dict["A_obtectus"])

    # plot_busco_stats(native_busco_stats_dict, orthoDB_busco_stats_dict, species_names)
    data = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/"
    plot_busco_stats(native_busco_stats_dict, orthoDB_busco = orthoDB_busco_stats_dict, species_names = species_names, outfile_name=f"{data}busco_stats_uniform_repeatmasking.png")
    # plot_busco_stats(native_busco_stats_dict, species_names = species_names, outfile_name="busco_stats_test.png")