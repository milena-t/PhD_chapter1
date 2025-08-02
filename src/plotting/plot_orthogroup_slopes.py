"""
for all orthogroups of choice:
    compute the slope of the orthogroup (14 datapoints, one for each species with [GF size, genome size])
plot the slope of all of those together, and highlight OGs with high slopes
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import parse_gff as gff
import analyze_multiple_CAFE_runs as CAFE
import parse_orthogroups as OGs
import compute_PIC as PIC


def plot_slopes(GF_sizes_dict, species_list, exp_dict, x_label, tree_path, filename = "sig_OGs_inclines.png", color_category = "orthoDB", percentile = 99, sig_list = [], log10_GF=False, log2_GF=True, correct_bh = False):
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

    inclines = {}
    intercepts = {}
    p_values = {}
    std_errs = {}
    return_dict = {}

    OG_sizes = {}
    
    for orthogroup, GF_sizes in tqdm(GF_sizes_dict.items()):
        
        ## log-transform the GF sizes
        if log10_GF:
            GF_sizes_dict = {species : np.log10(GF_sizes[species]) if species in GF_sizes else 0 for species in species_list }
        elif log2_GF:
            GF_sizes_dict = {species : np.log2(GF_sizes[species]) if species in GF_sizes else 0 for species in species_list }
        else:
            GF_sizes_dict = {species : GF_sizes[species] if species in GF_sizes else 0 for species in species_list }

        ## calculate PICs
        PICs_GF_sizes = PIC.calculate_PIC(tree_path=tree_path, trait_values=GF_sizes_dict)
        x_axis_vec = PIC.calculate_PIC(tree_path=tree_path, trait_values=exp_dict)

        # m, b = np.polyfit(x_axis_vec, GF_sizes_vec, 1)
        result = scipy.stats.linregress(x_axis_vec, PICs_GF_sizes)
        inclines[orthogroup] = result.slope
        intercepts[orthogroup] = result.intercept
        p_values[orthogroup] = result.pvalue
        std_errs[orthogroup] = result.stderr
        return_dict[orthogroup] = [result.slope, result.pvalue, "x"]
        
        OG_size = sum(list(GF_sizes.values()))
        OG_sizes[orthogroup] = OG_size

        # if OG_size>200:
        #     print(f" --> {orthogroup} : size {OG_size}, p-value {result.pvalue:.2f}, \n\t GF sizes : {GF_sizes_dict[orthogroup]}")

    if sig_list==[]:
        inclines_list = list(inclines.values())
        percentile_upper = np.percentile(inclines_list, q = percentile)
        percentile_lower = np.percentile(inclines_list, q = 100-percentile)
        print(f"max incline: {max(inclines_list):.5f},\t\tmin_incline: {min(inclines_list):.5f} \n{percentile}th percentile: {percentile_upper:.5f},\t{100-percentile}th percentile: {percentile_lower:.5f}")

    fs = 22 # set font size
    
    ylab="mean number of orthogroup members"
    # get a list of lists with [native, orthoDB] number of gene families per species

    colors = {
        "native" : "#b82946", # red
        "native_unsignificant" : "#DE6880", #light red
        "native_multiple_testing_sig" : "#861D32", #dark orange
        "orthoDB": "#F2933A", # orange
        "orthoDB_unsignificant" : "#F6B679", # light orange
        "orthoDB_multiple_testing_sig" : "#C0630C", #dark orange
        "background" : "#838383" # grey
    }

    not_significant_lines = 0
    significant_lines = 0
    if sig_list==[]:
        for orthogroup, incline in tqdm(inclines.items()):
            if incline > percentile_upper or incline < percentile_lower:
                significant_lines+=1
            else:
                not_significant_lines+=1
                # ax.plot(x_axis_vec, [incline*x_value+intercepts[orthogroup] for x_value in x_axis_vec], color = colors["background"], linewidth = 0.5)
    else:
        for orthogroup, incline in tqdm(inclines.items()):
            if orthogroup in sig_list:
                significant_lines+=1
            else:
                not_significant_lines+=1

    ################################
    # make scatterplot of regression inclines and orthogroup sizes
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params(axis='both', which='major', labelsize=fs)

    if sig_list == []:
        orthogroups_list = list(intercepts.keys())
        inclines_list = [inclines[orthogroup] for orthogroup in orthogroups_list]
        OG_sizes_list = [OG_sizes[orthogroup] for orthogroup in orthogroups_list]
        colors_list = [colors[color_category] if incline > percentile_upper or incline < percentile_lower else colors["background"] for incline in inclines_list]
        if len(inclines_list) != len(OG_sizes_list):
            raise RuntimeError(f"slopes list: {len(inclines_list)}\nOG sizes list: {len(OG_sizes_list)}")
        ax.scatter(OG_sizes_list, inclines_list, color = colors_list, s=75)
    else:
        # inclines_list = [inclines[orthogroup] for orthogroup in sig_list]
        inclines_sig_list = []
        OG_sizes_sig_list = []
        # OG_sizes_list = [OG_sizes[orthogroup] for orthogroup in sig_list]
        inclines_unsig_list = []
        OG_sizes_unsig_list = []

        if correct_bh == False:
            for orthogroup in sig_list:
                p_val = p_values[orthogroup]
                if p_val<0.05:
                    inclines_sig_list.append(inclines[orthogroup])
                    OG_sizes_sig_list.append(OG_sizes[orthogroup])
                else:   
                    inclines_unsig_list.append(inclines[orthogroup])
                    OG_sizes_unsig_list.append(OG_sizes[orthogroup])
        
            ax.scatter(OG_sizes_unsig_list, inclines_unsig_list, color = colors[f"{color_category}_unsignificant"], s=30, marker = "x", label = "unsignificant")# with marker="o" use facecolors = "none" to make an un-filled circle
            ax.scatter(OG_sizes_sig_list, inclines_sig_list, color = colors[color_category], s=75, label = "significant")

        elif correct_bh:
            inclines_bh_cor_sig_list = []
            OG_sizes_bh_cor_sig_list = []
            p_values_list = [p_values[orthogroup] for orthogroup in sig_list]
            reject, p_values_bh, _, _ = multipletests(p_values_list, alpha=0.05, method='fdr_bh')
            # p_values_bh = scipy.stats.false_discovery_control(p_values_list) # default benjamini-hochberg correction

            for i, orthogroup in enumerate(sig_list):
                p_val = p_values_list[i]
                p_cor = p_values_bh[i]
                if p_cor < 0.05:
                    print(f"\t\t-- {orthogroup}")
                    inclines_bh_cor_sig_list.append(inclines[orthogroup])
                    OG_sizes_bh_cor_sig_list.append(OG_sizes[orthogroup])
                    return_dict[orthogroup][-1] = "y"
                elif p_val < 0.05:
                    inclines_sig_list.append(inclines[orthogroup])
                    OG_sizes_sig_list.append(OG_sizes[orthogroup])
                    return_dict[orthogroup][-1] = "n"
                else:   
                    inclines_unsig_list.append(inclines[orthogroup])
                    OG_sizes_unsig_list.append(OG_sizes[orthogroup])
                    return_dict[orthogroup][-1] = "n"
            
            ax.scatter(OG_sizes_unsig_list, inclines_unsig_list, color = colors[f"{color_category}_unsignificant"], s=30, marker = "x", label = "unsignificant")# with marker="o" use facecolors = "none" to make an un-filled circle
            ax.scatter(OG_sizes_sig_list, inclines_sig_list, color = colors[color_category], s=75, label = "significant")
            ax.scatter(OG_sizes_bh_cor_sig_list, inclines_bh_cor_sig_list, color = colors[f"{color_category}_multiple_testing_sig"], s=75, marker="v", label = "B.H. corrected")

    if sig_list==[]:
        x_text_coord = max(OG_sizes_list)
        x_text_coord = 0.8*x_text_coord
        ax.axhline(percentile_upper, linestyle='--', color = colors["background"])
        ax.text(x_text_coord, percentile_upper+percentile_upper*1.25, f'{percentile_upper:.3f}', va='top', ha="center", fontsize=fs, color=colors["background"])
        ax.axhline(percentile_lower, linestyle='--', color = colors["background"])
        ax.text(x_text_coord, percentile_lower-percentile_upper*0.25, f'{percentile_lower:.3f}', va='top', ha="center", fontsize=fs, color=colors["background"])
    
    if sig_list==[]:
        ylab = f"regression slopes of individual orthogroups \n(color by {100-percentile}th and {percentile}th percentile, {significant_lines} of \n{len(inclines)} orthogroups outside percentile bounds)"
    else:
        ylab = f"regression slopes of individual orthogroups \n(only {len(sig_list)} significant orthogroups shown)"
    ax.set_ylabel(ylab, fontsize = fs)
    ax.set_xlabel("orthogroup size", fontsize = fs)
    title_ = x_label.split(" in")[0]
    if log10_GF:
        title = f"log10(Gene family size) vs. {title_}"
    elif log2_GF:
        title = f"log2(Gene family size) vs. {title_}"
    else:
        title = f"Gene family size vs. {title_}"
    plt.title(title, fontsize=fs*1.2)
    ax.legend(fontsize = fs, loc='lower right', title_fontsize = fs)

    filename = filename.split(".png")[0]
    filename = f"{filename}_vs_OG_size"
    if sig_list==[]:
        filename = f"{filename}_{percentile}th_percentile_colors.png"
    else:
        filename = f"{filename}_sig_OGs_colors.png"

    plt.savefig(filename, dpi = 300, transparent = True, bbox_inches='tight')
    print("Figure with slopes and OG sizes saved in the current working directory directory as: "+filename)
    
    return return_dict



if __name__ == "__main__":

    try:
        tree = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_ultrametric.nw"
        species_names = gff.make_species_order_from_tree(tree)
        data_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/"
        CAFE_runs_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence"
    except:
        tree = "/Users/milena/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_ultrametric.nw"
        species_names = gff.make_species_order_from_tree(tree)
        data_dir = "/Users/milena/work/PhD_code/PhD_chapter1/data/"
        CAFE_runs_dir = "/Users/milena/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence"

    orthogroups_native_filepath = f"{data_dir}orthofinder_native/N0.tsv"
    orthogroups_orthoDB_filepath = f"{data_dir}orthofinder_uniform/N0.tsv"
    # sig_native = f"{data_dir}CAFE_native_Base_Family_results.txt"
    # sig_orthoDB = f"{data_dir}CAFE_uniform_Base_Family_results.txt"

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
                    # "C_analis" : 971,
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
                        # "C_analis" : 69.29,
                        "C_maculatus" : 73.49
                        }
    
    if False:
        print(f"\n\tnative")
        native_dict_lists = OGs.parse_orthogroups_dict(orthogroups_native_filepath)
        native_dict = OGs.get_GF_sizes(native_dict_lists)
        native_sig_list, native_cafe_list = OGs.get_sig_orthogroups(sig_native)

        plot_slopes(GF_sizes_dict=native_dict, species_list = species_names, exp_dict=genome_sizes_dict, x_label = "Genome size in Mb", color_category="native", filename = f"{data_dir}sig_OGs_vs_GS_inclines.png")
        plot_slopes(GF_sizes_dict=native_dict, species_list = species_names, exp_dict=repeat_percentages, x_label = "Repeat content in percent", color_category="native", filename = f"{data_dir}sig_OGs_vs_reps_inclines.png")


    print(f"\n\torthoDB")
    # orthoDB_sig_list, orthoDB_cafe_list = OGs.get_sig_orthogroups(sig_orthoDB)
    orthoDB_sig_list, orthoDB_cafe_list = CAFE.get_overlap_OG_sig_list(CAFE_runs_dir)
    print(f"{len(orthoDB_sig_list)} significant orthogroups out of {len(orthoDB_cafe_list)} in total")
    orthoDB_dict_lists = OGs.parse_orthogroups_dict(orthogroups_orthoDB_filepath, orthoDB_cafe_list)
    orthoDB_dict = OGs.get_GF_sizes(orthoDB_dict_lists)

    print(f"\n\t\t * Genome size")
    # plot_slopes(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=genome_sizes_dict, x_label = "Genome size in Mb", filename = f"{data_dir}sig_OGs_vs_GS_inclines.png")
    GS_inclines = plot_slopes(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=genome_sizes_dict, x_label = "Genome size in Mb", tree_path = tree,  filename = f"{data_dir}sig_OGs_vs_GS_inclines_bh_corrected_PIC.png", sig_list=orthoDB_sig_list, correct_bh=True)
    gff.write_dict_to_file(GS_inclines, f"{data_dir}sig_OGs_vs_GS_inclines_pvalues.tsv", header=f"OG\tslope\tp-value\tsig_after_multiple_testing", separator="\t")

    print(f"\n\t\t * repeat content")
    # plot_slopes(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=repeat_percentages, x_label = "Repeat content in percent", filename = f"{data_dir}sig_OGs_vs_reps_inclines.png")
    TE_inclines = plot_slopes(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=repeat_percentages, x_label = "Repeat content in percent", tree_path = tree,  filename = f"{data_dir}sig_OGs_vs_reps_inclines_bh_corrected_PIC.png", sig_list=orthoDB_sig_list, correct_bh=True)
    gff.write_dict_to_file(TE_inclines, f"{data_dir}sig_OGs_vs_reps_inclines_pvalues.tsv", header=f"OG\tslope\tp-value\tsig_after_multiple_testing", separator="\t")

    ## last column of the sig_OGs_[...]_pvalues.tsv lists has one of three:
    #  * x: the orthogorup is not significant according to CAFE
    #  * n: the orthogroup is not significantly correlated after multiple testing correction
    #  * y: the orthogroup is significantly correlated after multiple testing correction