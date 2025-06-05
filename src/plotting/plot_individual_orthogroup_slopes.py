"""
for all orthogroups of choice:
    compute the slope of the orthogroup (14 datapoints, one for each species with [GF size, genome size])
plot the slope of all of those together, and highlight OGs with high slopes
"""

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
import parse_gff as gff
import parse_orthogroups as OGs


def plot_slopes(GF_sizes_dict, species_list, exp_dict, x_label, filename = "sig_OGs_inclines.png", color_category = "orthoDB", percentile = 99):
    """
    Plot fitted linear regression for each significant orthogroup.
    exp_dict is the dictionary with the x-axis variables, like genome size or repeat content
    returns a dictionary with { orthogroupID : incline }
    the percentile option lets you show only upper and lower percentile of inclines
    """
    if "/" in filename and color_category not in filename:
        filename_base = filename.split("/")[-1]
        filename_path = "/".join(filename.split("/")[:-1])
        filename = f"{filename_path}/{color_category}_{filename_base}"
    else:
        filename = f"{color_category}_{filename}"

    fs = 15 # set font size
    # plot each column in the dataframe as a line in the same plot thorugh a for-loop
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    
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
    
    inclines = {}
    intercepts = {}
    OG_sizes = {}

    for orthogroup, GF_sizes in tqdm(GF_sizes_dict.items()):
        GF_sizes_vec = [GF_sizes[species] if species in GF_sizes else 0 for species in species_list ]
        x_axis_vec = [exp_dict[species] for species in species_list]

        m, b = np.polyfit(x_axis_vec, GF_sizes_vec, 1)
        inclines[orthogroup] = m
        intercepts[orthogroup] = b

        OG_sizes[orthogroup] = sum(list(GF_sizes.values()))

    inclines_list = list(inclines.values())
    percentile_upper = np.percentile(inclines_list, q = percentile)
    percentile_lower = np.percentile(inclines_list, q = 100-percentile)
    print(f"max incline: {max(inclines_list):.5f},\t\tmin_incline: {min(inclines_list):.5f} \n95th percentile: {percentile_upper:.5f},\t5th percentile: {percentile_lower:.5f}")

    not_plotted_lines = 0
    for orthogroup, incline in tqdm(inclines.items()):
        if incline > percentile_upper or incline < percentile_lower:
            # print(f"{orthogroup} : incline {incline:.4f}")
            ax.plot(x_axis_vec, [incline*x_value+intercepts[orthogroup] for x_value in x_axis_vec], color = colors[color_category], linewidth = 1)
        else:
            not_plotted_lines+=1
            # ax.plot(x_axis_vec, [incline*x_value+intercepts[orthogroup] for x_value in x_axis_vec], color = colors["background"], linewidth = 0.5)

    ylab = f"gene family size (only {100-percentile}th and {percentile}th percentile shown, \n{len(inclines)-not_plotted_lines} of {len(inclines)} orthogroups)"
    ax.set_ylabel(ylab, fontsize = fs)
    ax.set_xlabel(x_label, fontsize = fs)
    title_ = x_label.split(" in")[0]
    title = f"regression lines of gene family size vs. {title_}"
    plt.title(title, fontsize=fs*1.2)

    plt.tight_layout()
    plt.savefig(filename, dpi = 300, transparent = True)
    print(f"{len(inclines)-not_plotted_lines} of {len(inclines)} orthogroups not in upper or lower {percentile}th percentile were not plotted")
    print("Figure with regression lines saved in the current working directory directory as: "+filename)

    ################################
    # make scatterplot of regression inclines and orthogroup sizes
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params(axis='both', which='major', labelsize=fs)

    orthogroups_list = list(intercepts.keys())
    inclines_list = [inclines[orthogroup] for orthogroup in orthogroups_list]
    OG_sizes_list = [OG_sizes[orthogroup] for orthogroup in orthogroups_list]
    colors_list = [colors[color_category] if incline > percentile_upper or incline < percentile_lower else colors["background"] for incline in inclines_list]

    ax.scatter(OG_sizes_list, inclines_list, color = colors_list)
    ylab = f"regression slope (color by {100-percentile}th and {percentile}th percentile, \n{len(inclines)-not_plotted_lines} of {len(inclines)} orthogroups)"
    ax.set_ylabel(ylab, fontsize = fs)
    ax.set_xlabel("orthogroup size", fontsize = fs)
    title = f"regression slopes of gene family size vs. {title_}"
    plt.title(title, fontsize=fs*1.2)

    plt.tight_layout()

    filename = filename.split(".png")[0]
    filename = f"{filename}_vs_OG_size.png"
    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure with slopes and OG sizes saved in the current working directory directory as: "+filename)
    
    return inclines



if __name__ == "__main__":
    tree = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_only_species_names.nw"
    species_names = gff.make_species_order_from_tree(tree)

    data_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/"
    orthogroups_native_filepath = f"{data_dir}orthofinder_native/N0.tsv"
    orthogroups_orthoDB_filepath = f"{data_dir}orthofinder_uniform/N0.tsv"
    sig_native = f"{data_dir}CAFE_native_Base_Family_results.txt"
    sig_orthoDB = f"{data_dir}CAFE_uniform_Base_Family_results.txt"

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
    repeat_percentages = {"D_melanogaster" : 22.17,
                        "I_luminosus" : 41.85,
                        "P_pyralis" : 49.90,
                        "C_septempunctata" : 61.54,
                        "A_verrucosus" : 26.63,
                        "T_castaneum" : 55.20,
                        "T_molitor" : 40.12,
                        "Z_morio" : 51.94,
                        "R_ferrugineus" : 49.64,
                        "D_ponderosae" : 20.36,
                        "A_obtectus" : 72.39,
                        "B_siliquastri" : 45.20,
                        "C_chinensis" : 59.33,
                        "C_analis" : 69.29,
                        "C_maculatus" : 73.49
                        }
    

    print(f"\n\tnative")
    native_dict_lists = OGs.parse_orthogroups_dict(orthogroups_native_filepath)
    native_dict = OGs.get_GF_sizes(native_dict_lists)
    native_sig_list, native_cafe_list = OGs.get_sig_orthogroups(sig_native)

    plot_slopes(GF_sizes_dict=native_dict, species_list = species_names, exp_dict=genome_sizes_dict, x_label = "Genome size in Mb", color_category="native", filename = f"{data_dir}sig_OGs_vs_GS_inclines.png")
    plot_slopes(GF_sizes_dict=native_dict, species_list = species_names, exp_dict=repeat_percentages, x_label = "Repeat content in percent", color_category="native", filename = f"{data_dir}sig_OGs_vs_reps_inclines.png")


    print(f"\n\torthoDB")
    orthoDB_dict_lists = OGs.parse_orthogroups_dict(orthogroups_orthoDB_filepath)
    orthoDB_dict = OGs.get_GF_sizes(orthoDB_dict_lists)
    orthoDB_sig_list, orthoDB_cafe_list = OGs.get_sig_orthogroups(sig_orthoDB)

    plot_slopes(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=genome_sizes_dict, x_label = "Genome size in Mb", filename = f"{data_dir}sig_OGs_vs_GS_inclines.png")
    plot_slopes(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=repeat_percentages, x_label = "Repeat content in percent", filename = f"{data_dir}sig_OGs_vs_reps_inclines.png")