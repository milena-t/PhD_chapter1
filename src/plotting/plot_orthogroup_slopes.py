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


def calculate_OG_lin_reg(GF_sizes:dict, exp_dict:dict, tree_path:str, species_list:list, log10_GF=False, log2_GF=True):
    """
    Calculate the linear regression of an orthogroup, represented as a dictionary { species : gene family members }
    I need the tree to calculate phylogenetically independent contrasts
    """
    
    log_possible = True

    ## log-transform the GF sizes
    if log10_GF:
        GF_sizes_dict = {species : np.log10(GF_sizes[species]) if species in GF_sizes else 0 for species in species_list }
    elif log2_GF:
        GF_sizes_dict = {species : np.log2(GF_sizes[species]) if species in GF_sizes else 0 for species in species_list }
    else:
        GF_sizes_dict = {species : GF_sizes[species] if species in GF_sizes else 0 for species in species_list }
    
    # if any value is 0 the log doesn't work
    if log10_GF or log2_GF:
        if len([exp_val for exp_val in exp_dict.values() if exp_val== 0.0])>0:
            # print(f"can't log the explanatory variable because one of the values is 0!")
            log10_GF = False
            log2_GF = False
            log_possible=False

    ## log transform explanatory variables
    if log10_GF:
        exp_dict = {species : np.log10(exp_dict[species]) if species in exp_dict else 0 for species in species_list }
    elif log2_GF:
        exp_dict = {species : np.log2(exp_dict[species]) if species in exp_dict else 0 for species in species_list }
    else:
        exp_dict = {species : exp_dict[species] if species in exp_dict else 0 for species in species_list }

    ## calculate PICs
    PICs_GF_sizes = PIC.calculate_PIC(tree_path=tree_path, trait_values=GF_sizes_dict)
    x_axis_vec = PIC.calculate_PIC(tree_path=tree_path, trait_values=exp_dict)

    ## linear regression
    # abs(min(x_axis_vec)-max(x_axis_vec))
    # print(f"{}")
    result = scipy.stats.linregress(x_axis_vec, PICs_GF_sizes)
    
    return result,PICs_GF_sizes,x_axis_vec,log_possible


def test_normality_of_residuals(result,PICs_GF_sizes,x_axis_vec):
    """
    test the normality of residuals from a linear regression created with scipy.stats.linregress
    """
    ## test normality of residuals
    def predict(x):
        pred_PIC = x*result.slope + result.intercept
        return(pred_PIC)

    residuals = [PICs_GF_sizes[i] - predict(x_axis_vec[i]) for i in range(len(x_axis_vec))]
    stat, p_value = scipy.stats.shapiro(residuals)
    
    return stat,p_value
        



def get_plot_values(GF_sizes_dict, species_list, exp_dict, sig_list, tree_path, log10_GF=False, log2_GF=True):
    """
    calculate fitted linear regression for each significant orthogroup.
    exp_dict is the dictionary with the x-axis variables, like genome size or repeat content
    returns a dictionary with { orthogroupID : incline }
    !! includes FDR multiple testing correction !!
    """

    inclines = {}
    intercepts = {}
    p_values = {}
    std_errs = {}
    return_dict = {}
    OG_sizes = {}

    excluded_OGs = []
    included_OGs = []
    
    # for orthogroup in tqdm(sig_list):
    for orthogroup in sig_list:
        
        GF_sizes = GF_sizes_dict[orthogroup]
        result,PICs_GF_sizes,x_axis_vec,log_possible = calculate_OG_lin_reg(GF_sizes = GF_sizes, exp_dict= exp_dict, tree_path = tree_path, species_list = species_list, log10_GF=log10_GF, log2_GF=log2_GF)


        stat,p_value = test_normality_of_residuals(result,PICs_GF_sizes,x_axis_vec)
        
        ## if residuals not normal don't include this orthogroup in the analysis
        if p_value < 0.05:
            excluded_OGs.append(orthogroup)
        
        else:
            included_OGs.append(orthogroup)

            inclines[orthogroup] = result.slope
            intercepts[orthogroup] = result.intercept
            p_values[orthogroup] = result.pvalue
            std_errs[orthogroup] = result.stderr
            return_dict[orthogroup] = [result.slope, result.pvalue, "x"]

            OG_sizes[orthogroup] = sum([GF_sizes[species] if species in GF_sizes else 0 for species in species_list ])
        
    ## DO multiple testing correction
    p_values_list = [p_values[orthogroup] for orthogroup in included_OGs]
    reject, p_values_bh, _, _ = multipletests(p_values_list, alpha=0.05, method='fdr_bh')

    p_values_BH = {}
    for i, orthogroup in enumerate(included_OGs):
            p_values_BH[orthogroup] = p_values_bh[i]

    return inclines,intercepts,p_values,p_values_BH,std_errs,return_dict,OG_sizes,included_OGs,log_possible



def read_repeat_categories(path:str):
    """
    read the csv table for the repeat categories
    {
        repeat_category : { species  : percentage }
    }
    """

    with open(path, "r") as infile:
        lines = infile.readlines()
        header = lines[0].strip().split(",")
        categories = header[1:]
        repeat_categories = { category : {} for category in categories}

        for line in lines[1:]:
            line = line.strip().split(",")
            species = line[0]
            category_values = line[1:]
            for i,category_value in enumerate(category_values):
                repeat_categories[categories[i]][species] = float(category_value)

    return repeat_categories




def plot_slopes(inclines,intercepts,p_values,p_values_bh,std_errs,return_dict,OG_sizes, sig_list ,x_label, filename = "sig_OGs_inclines.png", color_category = "orthoDB", percentile = 99, log10_GF=False, log2_GF=True, log_possible=True, svg = False):

    ### PLOT 

    if "/" in filename and color_category not in filename:
        filename_base = filename.split("/")[-1]
        filename_path = "/".join(filename.split("/")[:-1])
        filename = f"{filename_path}/{color_category}_{filename_base}"
    else:
        filename = f"{color_category}_{filename}"

    not_significant_lines = 0
    significant_lines = 0

    for orthogroup, incline in tqdm(inclines.items()):
        if orthogroup in sig_list:
            significant_lines+=1
        else:
            not_significant_lines+=1

    ### PLOT parameters
    fs = 22 # set font size
    
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1, 1, 1)
    ax.tick_params(axis='both', which='major', labelsize=fs)

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


    # inclines_list = [inclines[orthogroup] for orthogroup in sig_list]
    inclines_sig_list = []
    OG_sizes_sig_list = []
    # OG_sizes_list = [OG_sizes[orthogroup] for orthogroup in sig_list]
    inclines_unsig_list = []
    OG_sizes_unsig_list = []
    # multiple testing correction    
    inclines_bh_cor_sig_list = []
    OG_sizes_bh_cor_sig_list = []

    for i, orthogroup in enumerate(sig_list):
        p_val = p_values[orthogroup]
        p_cor = p_values_bh[orthogroup]
        if p_cor < 0.05:
            print(f"\t\t-- {orthogroup} significant after FDR correction")
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

    ylab = f"regression slopes of individual orthogroups \n(only {len(sig_list)} significant orthogroups shown)"
    ax.set_ylabel(ylab, fontsize = fs)
    ax.set_xlabel("orthogroup size", fontsize = fs)
    title_ = x_label.split(" in")[0]
    if log_possible:
        title_ = f"log2({title_})"
    if log10_GF:
        title = f"lin. reg. {title_} vs. \nlog10(Gene family size) "
    elif log2_GF:
        title = f"lin. reg. {title_} vs. \nlog2(Gene family size) "
    else:
        title = f"lin. reg. {title_} vs. \nGene family size"
    plt.title(title, fontsize=fs*1.2)
    ax.legend(fontsize = fs, loc='lower right', title_fontsize = fs)

    filename = filename.split(".png")[0]
    filename = f"{filename}_vs_OG_size"
    if sig_list==[]:
        filename = f"{filename}_{percentile}th_percentile_colors.png"
    else:
        filename = f"{filename}_sig_OGs_colors.png"

    if svg:
        filename = filename.replace(".png", ".svg")

    # plt.show()
    plt.savefig(filename, dpi = 200, transparent = True, bbox_inches='tight')
    print("Figure with slopes and OG sizes saved in the current working directory directory as: "+filename)
    
    return return_dict





if __name__ == "__main__":

    try:
        tree = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_ultrametric.nw"
        species_names = gff.make_species_order_from_tree(tree)
        data_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/"
        CAFE_runs_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence"
        repeat_categories_in_species = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/whole_species_repeat_categories.csv"
    except:
        tree = "/Users/milena/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_ultrametric.nw"
        species_names = gff.make_species_order_from_tree(tree)
        data_dir = "/Users/milena/work/PhD_code/PhD_chapter1/data/"
        CAFE_runs_dir = "/Users/milena/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence"
        repeat_categories_in_species = "/Users/milena/work/PhD_code/PhD_chapter1/data/whole_species_repeat_categories.csv"

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


    if True:

        print(f"\n\torthoDB")
        # orthoDB_sig_list, orthoDB_cafe_list = OGs.get_sig_orthogroups(sig_orthoDB)
        orthoDB_sig_list, orthoDB_cafe_list = CAFE.get_overlap_OG_sig_list(CAFE_runs_dir)
        print(f"{len(orthoDB_sig_list)} significant orthogroups out of {len(orthoDB_cafe_list)} in total")
        orthoDB_dict_lists = OGs.parse_orthogroups_dict(orthogroups_orthoDB_filepath, orthoDB_cafe_list)
        orthoDB_dict = OGs.get_GF_sizes(orthoDB_dict_lists)

    if True:

        svg_bool = True

        species_names.remove("D_melanogaster")

        if False:
            print(f"\n\t\t * Genome size")
            inclines,intercepts,p_values,p_values_BH,std_errs,return_dict,OG_sizes,included_OGs,log_possible = get_plot_values(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=genome_sizes_dict, sig_list=orthoDB_sig_list, tree_path=tree)
            print(f"{len(included_OGs)} (of {len(orthoDB_sig_list)}) orthogroups included because the LR residuals are not normally distributed")
            GS_inclines = plot_slopes(inclines,intercepts,p_values,p_values_BH,std_errs,return_dict,OG_sizes, x_label = "Genome size in Mb",  filename = f"{data_dir}correlations/sig_OGs_vs_GS_inclines_bh_corrected_PIC.png", sig_list=included_OGs,log_possible=log_possible, svg=svg_bool)
            # gff.write_dict_to_file(GS_inclines, f"{data_dir}sig_OGs_vs_GS_inclines_pvalues.tsv", header=f"OG\tslope\tp-value\tsig_after_multiple_testing", separator="\t")

            print(f"\n\t\t * repeat content")
            inclines,intercepts,p_values,p_values_BH,std_errs,return_dict,OG_sizes,included_OGs,log_possible = get_plot_values(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=repeat_percentages, sig_list=orthoDB_sig_list, tree_path=tree)
            print(f"{len(included_OGs)} (of {len(orthoDB_sig_list)}) orthogroups included because the LR residuals are not normally distributed")
            TE_inclines = plot_slopes(inclines,intercepts,p_values,p_values_BH,std_errs,return_dict,OG_sizes, x_label = "Repeat content in percent",  filename = f"{data_dir}correlations/sig_OGs_vs_reps_inclines_bh_corrected_PIC.png", sig_list=included_OGs,log_possible=log_possible, svg=svg_bool)
            # gff.write_dict_to_file(TE_inclines, f"{data_dir}sig_OGs_vs_reps_inclines_pvalues.tsv", header=f"OG\tslope\tp-value\tsig_after_multiple_testing", separator="\t")

        ## do the individual repeat categories
        repeats_categories_dict = read_repeat_categories(repeat_categories_in_species)

        for repeat_category in repeats_categories_dict.keys():
            repeat_percentages = repeats_categories_dict[repeat_category]
            print(f"\n\t\t * repeat category: {repeat_category}")
            
            inclines,intercepts,p_values,p_values_BH,std_errs,return_dict,OG_sizes,included_OGs,log_possible = get_plot_values(GF_sizes_dict=orthoDB_dict, species_list = species_names, exp_dict=repeat_percentages, sig_list=orthoDB_sig_list, tree_path=tree)
            print(f"{len(included_OGs)} (of {len(orthoDB_sig_list)}) orthogroups included because the LR residuals are normally distributed")

            if len(included_OGs) == 0:
                raise RuntimeError

            repeat_category_filename = repeat_category.replace(" ", "")
            repeat_category_filename = repeat_category_filename.replace(".","")
            TE_inclines = plot_slopes(inclines,intercepts,p_values,p_values_BH,std_errs,return_dict,OG_sizes, x_label = f"{repeat_category} content in percent",  filename = f"{data_dir}correlations/sig_OGs_vs_{repeat_category_filename}_inclines_bh_corrected_PIC.png", sig_list=included_OGs,log_possible=log_possible, svg=svg_bool)


        ## last column of the sig_OGs_[...]_pvalues.tsv lists has one of three:
        #  * x: the orthogorup is not significant according to CAFE
        #  * n: the orthogroup is not significantly correlated after multiple testing correction
        #  * y: the orthogroup is significantly correlated after multiple testing correction

    ## Test stats stuff
    if False:

        repeats_categories_dict = read_repeat_categories(repeat_categories_in_species)
        
        species_names_no_Dmel = species_names
        species_names_no_Dmel.remove("D_melanogaster")

        ## count for all repeat categories in repeats_categories_dict
        for repeat_category in repeats_categories_dict.keys():
            count_all = 0
            count_non_normal_TE = 0
            repeat_percentages = repeats_categories_dict[repeat_category]
            pvalues = []
            for orthogroup, GF_sizes in orthoDB_dict.items():
                count_all += 1
                TE_result,PICs_GF_sizes,x_axis_vec,log_possible = calculate_OG_lin_reg(GF_sizes = GF_sizes, exp_dict= repeat_percentages, tree_path = tree, species_list = species_names_no_Dmel, log10_GF=False, log2_GF=True)
                TE_stat,TE_p_value = test_normality_of_residuals(TE_result,PICs_GF_sizes,x_axis_vec)
                
                if TE_p_value < 0.05:
                    count_non_normal_TE += 1
                pvalues.append(TE_result.pvalue)
            
            reject, p_values_bh, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')
            sig_p_val = len([p_value for p_value in p_values_bh if p_value < 0.05])

            percentage = 100*count_non_normal_TE/count_all
            print(f"repeat category: {repeat_category}\n\t -- {count_non_normal_TE} ({percentage:.2f} %) OGs not normal residuals")
            if sig_p_val>0:
                print(f"\t -- {sig_p_val} OGs significantly correlated after FDR correction")

            """
            repeat category: Retroelements
                    -- 3238 (38.94 %) OGs not normal residuals
            repeat category: DNA transposons
                    -- 3098 (37.26 %) OGs not normal residuals
            repeat category: Rolling-circles
                    -- 2105 (25.32 %) OGs not normal residuals
            repeat category: Unclassified
                    -- 3045 (36.62 %) OGs not normal residuals
            repeat category: Small RNA
                    -- 0 (0.00 %) OGs not normal residuals
            repeat category: Satellites
                    -- 0 (0.00 %) OGs not normal residuals
            repeat category: Simple repeats
                    -- 3294 (39.62 %) OGs not normal residuals
            repeat category: Low complexity
                    -- 3337 (40.13 %) OGs not normal residuals

            --> not a single significant correlation after multiple testing!!

            """            


        if False:
            for orthogroup, GF_sizes in tqdm(orthoDB_dict.items()):
                count_all += 1

                GS_result,PICs_GF_sizes,x_axis_vec,log_possible = calculate_OG_lin_reg(GF_sizes = GF_sizes, exp_dict= genome_sizes_dict, tree_path = tree, species_list = species_names_no_Dmel, log10_GF=False, log2_GF=True)
                GS_stat,GS_p_value = test_normality_of_residuals(GS_result,PICs_GF_sizes,x_axis_vec)
                if GS_p_value < 0.05:
                    count_non_normal_GS += 1

                TE_result,PICs_GF_sizes,x_axis_vec,log_possible = calculate_OG_lin_reg(GF_sizes = GF_sizes, exp_dict= repeat_percentages, tree_path = tree, species_list = species_names_no_Dmel, log10_GF=False, log2_GF=True)
                TE_stat,TE_p_value = test_normality_of_residuals(TE_result,PICs_GF_sizes,x_axis_vec)
                if TE_p_value < 0.05:
                    count_non_normal_TE += 1

            print(f"{count_all} orthogroups, linear models residuals calculated. \n\t -- GS: {count_non_normal_GS} not normal\n\t -- TE: {count_non_normal_TE} not normal")