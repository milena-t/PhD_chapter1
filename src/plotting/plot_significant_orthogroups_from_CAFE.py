# plot the orthogroup size per species for all orthogroups, and highlight the significant ones

import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from statistics import mean
import numpy as np
import pandas as pd
import src.parse_gff as gff


def make_species_order_from_tree(newick_tree_path):
    """
    Takes a newick tree and extracts a list of all leaf names in order of occurence
    Makes it possible to align x-axes in later plots with trees
    """
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

def read_orthogroups_input(filepath):
    """
    read the input tsv with the orthogroup numbers into a dict of dicts
    
    {
        orthogroup : {
            species1 : num,
            species2 : num,
            ...
        }
        orthogroup2 : {
            species1 : num,
            species2 : num,
            ...
        }
    }
    """
    with open(filepath, "r") as orthogroups_file:
        orthogroups_lines = orthogroups_file.readlines()
        species = orthogroups_lines[0].strip().split("\t")[2:]
        output = {}
        for orthogroup_line in orthogroups_lines[1:]:
            orthogroups = orthogroup_line.strip().split("\t")
            output[orthogroups[1]] = { species[i] : int(orthogroups[i+2])for i in range(len(species))}
    return(output)


def get_sig_orthogroups(filepath, p_sig = 0.05):
    """ 
    get a list of significant orthogroup IDs
    """
    sig_list = []
    cafe_list = []
    with open(filepath, "r") as file:
        next(file) # skip first line
        for line in file:
            orthogroup, p_value, sig_bool = line.strip().split("\t")
            if float(p_value)<p_sig:
                sig_list.append(orthogroup)
            cafe_list.append(orthogroup)

    return(sig_list, cafe_list)


def get_means(orthogroups_dict, sig_list, all_cafe_list, species_names = []):
    """
    get mean orthogroup sizes in all species for both significant and non-significant orthogroups
    """

    if len(species_names)==0:
        OGs = list(orthogroups_dict.keys())
        species_names = list(orthogroups_dict[OGs[0]].keys())

    all_unsig = {species : [] for species in species_names}
    all_sig = {species : [] for species in species_names}

    for orthogroup in sig_list:
        for species in species_names:
            all_sig[species].append(orthogroups_dict[orthogroup][species])
    
    for orthogroup in all_cafe_list:
        if orthogroup in sig_list:
            continue
        for species in species_names:
            all_unsig[species].append(orthogroups_dict[orthogroup][species])       
    
    ### print mean line for significant and non-significant orthogroups
    # not very aussagekräftig, unfortunately

    all_unsig = {species : mean(all_unsig[species]) for species in species_names}
    all_sig = {species : mean(all_sig[species]) for species in species_names}

    return(all_unsig, all_sig)


def read_whole_genome_stats(filepath):
    stats_dict = {}
    with open(filepath, "r") as file:
        file = file.readlines()
        data_headers = file[0].strip().split(",")[1:]
        for line in file[1:]:
            line = line.strip().split(",")
            # for i in range(len(data_headers)):
            #     print(f"{i} --> {data_headers[i]} : {line[i+1]}")
            stats_dict[line[0]] = {data_headers[i] : float(line[i+1]) for i in range(len(data_headers))}
    return(stats_dict)


def plot_gene_counts(orthogroups_dict, sig_list, all_cafe_list, species_names, annotation = "native", filename = "significant_orthogroups_from_CAFE.png"):
    fs = 15 # set font size
    # plot each column in the dataframe as a line in the same plot thorugh a for-loop
    fig = plt.figure(figsize=(10,8))
    
    ax = fig.add_subplot(1, 1, 1)
    
    ylab="number of orthogroup members"
    # get a list of lists with [native, orthoDB] number of gene families per species

    legend_labels = []
    ncol_legend = 0

    colors = {
        "native" : "#b82946", # red
        "orthoDB": "#F2933A", # orange
        "background" : "#838383"#030303" # black
    }

    ymax = 0

    # calculate the mean line
    all_unsig = {species : [] for species in species_names}
    all_sig = {species : [] for species in species_names}

    for orthogroup in sig_list:
        gene_family_members = [orthogroups_dict[orthogroup][species] for species in species_names]
        for species in species_names:
            all_unsig[species].append(orthogroups_dict[orthogroup][species])
        if max(gene_family_members) > ymax:
            ymax = max(gene_family_members)
        ax.plot(species_names, gene_family_members, color = colors[annotation], alpha = 0.4, linewidth = 0.8)
    
    for orthogroup in all_cafe_list:
        if orthogroup in sig_list:
            continue
        for species in species_names:
            all_sig[species].append(orthogroups_dict[orthogroup][species])
        gene_family_members = [orthogroups_dict[orthogroup][species] for species in species_names]
        if max(gene_family_members) > ymax:
            ymax = max(gene_family_members)
        ax.plot(species_names, gene_family_members, color = colors["background"], alpha = 0.25, linewidth = 0.8)
    
    ### print mean line for significant and non-significant orthogroups
    # not very aussagekräftig, unfortunately

    # all_unsig_list = [mean(all_unsig[species]) for species in species_names]
    # print(all_unsig_list)
    # ax.plot(species_names, all_unsig_list, color = colors["background"])
    # all_sig_list = [mean(all_sig[species]) for species in species_names]
    # print(all_sig_list)
    # ax.plot(species_names, all_sig_list, color = colors[annotation])

    ax.set_ylabel(ylab, fontsize = fs)
    plt.xticks(labels=[species.replace("_", ". ") for species in species_names], ticks=species_names, rotation = 90, fontsize = fs)
    
    handles = []
    labels = [] 
    handles.append(mpatches.Patch(color=colors[annotation]))
    labels.append("significant")
    handles.append(mpatches.Patch(color=colors["background"]))
    labels.append("non-significant")

    ax.legend(handles, labels, fontsize = fs, loc='upper right', title_fontsize = fs)

    #ax.legend(fontsize = fs)
    # set grid only for X axis ticks 
    ax.grid(True)
    ax.yaxis.grid(False)

    ymax = ymax*1.1
    ax.set_ylim(-5,ymax)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)



def plot_means(native, orthoDB, whole_genome_stats, species_names, x_category = "", filename = "mean_orthogroups_from_CAFE.png", return_table = True):
    fs = 15 # set font size
    # plot each column in the dataframe as a line in the same plot thorugh a for-loop
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1, 1, 1)
    
    ylab="mean number of orthogroup members"
    # get a list of lists with [native, orthoDB] number of gene families per species

    legend_labels = []
    ncol_legend = 0

    colors = {
        "native" : "#b82946", # red
        "native_unsignificant" : "#DE6880", #light red
        "orthoDB": "#F2933A", # orange
        "orthoDB_unsignificant" : "#F6B679", # light orange
        "background" : "#838383" # grey
    }

    native_unsigfnicant = [native["unsignificant"][species] for species in species_names]
    native_sigfnicant = [native["significant"][species] for species in species_names]    
    orthoDB_unsigfnicant = [orthoDB["unsignificant"][species] for species in species_names]
    orthoDB_sigfnicant = [orthoDB["significant"][species] for species in species_names]

    
    if len(x_category) == 0:
        ax.plot(species_names, native_sigfnicant, color = colors["native"], label = "significant (native)")
        ax.plot(species_names, native_unsigfnicant, color = colors["native"], label = "unsignificant (native)", linestyle = ":")
        ax.plot(species_names, orthoDB_sigfnicant, color = colors["orthoDB"], label = "significant (orthoDB)")
        ax.plot(species_names, orthoDB_unsigfnicant, color = colors["orthoDB"], label = "unsignificant (orthoDB)", linestyle = ":")
        plt.xticks(labels=[species.replace("_", ". ") for species in species_names], ticks=species_names, rotation = 90, fontsize = fs)
        ax.grid(True)
        ax.yaxis.grid(False)
    
    elif len(x_category) > 0:
        try:
            x_values = [whole_genome_stats[species][x_category] for species in species_names]
        except:
            for key, value in whole_genome_stats[species_names[0]]:
                print(key,value)
            raise RuntimeError(f"{x_category} is an invalid category for x axis.")
        ax.scatter(x_values, native_sigfnicant, color = colors["native"], label = "significant (native)")
        ax.scatter(x_values, native_unsigfnicant, color = colors["native"], label = "unsignificant (native)", marker = "o", facecolors = "none")
        ax.scatter(x_values, orthoDB_sigfnicant, color = colors["orthoDB"], label = "significant (orthoDB)")
        ax.scatter(x_values, orthoDB_unsigfnicant, color = colors["orthoDB"], label = "unsignificant (orthoDB)", marker = "o", facecolors = "none")
        # make regression lines

        m_nat, b_nat = np.polyfit(x_values, native_sigfnicant, 1)
        m_odb, b_odb = np.polyfit(x_values, orthoDB_sigfnicant, 1)
        # print(f"native incline: {m_nat}")
        # print(f"orthoDB incline: {m_odb}")
        ax.plot(x_values, [m_nat*x_value+b_nat for x_value in x_values], color = colors["native"], linewidth = 1 , label = f"reg. line (native) incline: {m_nat:.3f}")
        ax.plot(x_values, [m_odb*x_value+b_odb for x_value in x_values], color = colors["orthoDB"], linewidth = 1 , label = f"reg. line (orthoDB) incline: {m_odb:.3f}")

        x_header = x_category.replace("_", " ")
        if "repeat" in x_category:
            ax.set_xlabel(f"{x_header} in the genome", fontsize = fs)
        if "size" in x_category:
            ax.set_xlabel(f"{x_header} in Mb", fontsize = fs)

        if return_table:
            table_df = pd.DataFrame({
                "native_unsigfnicant_mean_OG_size" : native_unsigfnicant,
                "native_sigfnicant_mean_OG_size" : native_sigfnicant,
                "orthoDB_unsigfnicant_mean_OG_size" : orthoDB_unsigfnicant,
                "orthoDB_sigfnicant_mean_OG_size" : orthoDB_sigfnicant,
                x_category : x_values
            })
            csv_filename = filename.replace(".png",".csv")
            table_df.to_csv(csv_filename, index=False) 
            print(f"saved data as file in: {csv_filename}")

    ax.set_ylabel(ylab, fontsize = fs)

    ax.legend(fontsize = fs, loc='upper left', title_fontsize = fs)

    ymax = max([max(native_unsigfnicant), max(native_sigfnicant), max(orthoDB_unsigfnicant), max(orthoDB_sigfnicant)])
    ymax = ymax*1.25
    ax.set_ylim(0,ymax)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = False)
    print("Figure saved in the current working directory directory as: "+filename)



if __name__ == "__main__":

    tree = "/Users/miltr339/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw"
    species_names = make_species_order_from_tree(tree)

    orthogroups_native = "/Users/miltr339/Box Sync/code/CAFE/CAFE_input_native_from_N0.tsv"
    orthogroups_orthoDB = "/Users/miltr339/Box Sync/code/CAFE/CAFE_input_orthoDB_TE_filtered.tsv"
    sig_native = "/Users/miltr339/Box Sync/code/CAFE/native_from_N0_Base_family_results.txt"
    sig_orthoDB = "/Users/miltr339/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"

    genome_sizes_dict = {"D_melanogaster" : 180,
                    "I_luminosus" : 842,
                    "P_pyralis" : 471,
                    "C_septempunctata" : 399,
                    "A_verrucosus" : 250,
                    "T_castaneum" : 204,
                    "T_molitor" : 258,
                    "Z_morio" : 461,
                    "R_ferrugineus" : 589,
                    "D_ponderosae" : 223,
                    "A_obtectus" : 949,
                    "B_siliquastri" : 375,
                    "C_chinensis" : 701,
                    "C_analis" : 971,
                    "C_maculatus" : 1202 
                    }

    whole_genome_stats_filepath = "/Users/miltr339/Box Sync/code/annotation_pipeline/repeatmasking_eval/eval_stats_3_annots_with_genome_size.csv"
    whole_genome_stats = read_whole_genome_stats(whole_genome_stats_filepath)

    native_dict = read_orthogroups_input(orthogroups_native)
    native_sig_list, native_cafe_list = get_sig_orthogroups(sig_native)
    print(f"{len(native_cafe_list)} of {len(native_dict)} native orthogroups considered in CAFE, of which {len(native_sig_list)} are significantly changing")
    all_unsig_native, all_sig_native = get_means(native_dict, native_sig_list, native_cafe_list, species_names)
    # plot_gene_counts(native_dict, native_sig_list, native_cafe_list, species_names, annotation = "native", filename = "native_significant_orthogroups_from_CAFE.png")

    orthoDB_dict = read_orthogroups_input(orthogroups_orthoDB)
    orthoDB_sig_list, orthoDB_cafe_list = get_sig_orthogroups(sig_orthoDB)
    print(f"{len(orthoDB_cafe_list)} of {len(orthoDB_dict)} orthoDB orthogroups considered in CAFE, of which {len(orthoDB_sig_list)} are significantly changing")
    all_unsig_orthoDB, all_sig_orthoDB = get_means(orthoDB_dict, orthoDB_sig_list, orthoDB_cafe_list, species_names)
    # plot_gene_counts(orthoDB_dict, orthoDB_sig_list, orthoDB_cafe_list, species_names, annotation = "orthoDB", filename = "orthoDB_significant_orthogroups_from_CAFE.png")

    native_means = {
        "unsignificant" : all_unsig_native,
        "significant" : all_sig_native
    }
    orthoDB_means = {
        "unsignificant" : all_unsig_orthoDB, 
        "significant" : all_sig_orthoDB
    }

    # plot_means(native_means, orthoDB_means, whole_genome_stats, species_names)
    plot_means(native_means, orthoDB_means, whole_genome_stats, species_names, x_category="genome_size", filename="mean_orthogroups_from_CAFE_vs_genome_size.png")
    plot_means(native_means, orthoDB_means, whole_genome_stats, species_names, x_category="repeat_percentage", filename="mean_orthogroups_from_CAFE_vs_repeats.png")