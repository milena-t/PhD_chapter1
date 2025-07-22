"""
After running CAFE5 20 times, I get between 530 and 570 significant orthogroups in each run. 
Lambda converges nicely to circa 0.4067

I will try to  figure out a good way to select a set of significantly rapidly evolving orthogroups from all runs
"""

import os
import parse_orthogroups as OGs
import matplotlib.pyplot as plt
try:
    import plotting.plot_significant_orthogroups_from_CAFE as plot_OG
except:
    import plot_significant_orthogroups_from_CAFE as plot_OG
import parse_gff as gff


def CAFE_output_paths(in_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence"):
    
    # get results files for all runs
    if in_dir[-1] == "/":
        in_dir = in_dir[:-1]
    runs_list = [] 
    for run in os.listdir(in_dir):
        name = f"{in_dir}/{run}/Base_family_results.txt"
        runs_list.append(name)
    
    return runs_list


def get_union_sig_list(cafe_outputs_list:str):
    """
    get a list of all orthogorups that are significant in at least one run
    """
    all_sig_OGs = []
    for run_path in cafe_outputs_list:
        sig_list, unsig_list = OGs.get_sig_orthogroups(run_path)
        all_sig_OGs.extend(sig_list)
    
    union_sig_OGs = list(set(all_sig_OGs))
    return union_sig_OGs


def count_OG_occurence(cafe_outputs_list:str):
    """
    count in how many runs each orthogroup is significant
    """
    all_sig_OGs = get_union_sig_list(cafe_outputs_list)
    all_sig_OGs = {sig_OG : 0 for sig_OG in all_sig_OGs}

    for run_path in cafe_outputs_list:
        sig_list, unsig_list = OGs.get_sig_orthogroups(run_path)
        for sig_OG in sig_list:
            all_sig_OGs[sig_OG]+=1
    return all_sig_OGs, unsig_list


def get_overlap_OG_sig_list(cafe_outputs_list):
    """
    get a list of orthogroups that are significant in all runs
    """
    if type(cafe_outputs_list) != list:
        # if the output is not a list yet
        cafe_outputs_list = CAFE_output_paths(cafe_outputs_list)

    sig_OG_counts, unsig_list = count_OG_occurence(cafe_outputs_list)
    num_runs = len(cafe_outputs_list)

    overlap_OGs = []
    for sig_OG, sig_runs in sig_OG_counts.items():
        if sig_runs == num_runs:
            overlap_OGs.append(sig_OG)
    return overlap_OGs, unsig_list



def plot_sig_OG_occurence_histogram(sig_counts_list:list[int], filename = ""):

    fig, ax = plt.subplots(1,1, figsize=(15, 12))
    fs = 35
    plt.rcParams.update({'font.size': fs})

    colors = {
        "uniform" : "#F2933A",
        "native" : "#b82946",
        "third" : "#9C4C32",
    }
    ax.set_xlabel('number of runs', fontsize=fs)
    ax.set_ylabel('number of OGs', fontsize=fs)
    ax.tick_params(axis='x', labelsize=fs)
    ax.tick_params(axis='y', labelsize=fs)
    plt.title(f'histogram of the number of runs an OG is significant in', fontsize = fs)
    ax.hist(sig_counts_list, bins = 20, color = colors["uniform"])

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)


CAFE_runs_dir = "/Users/milena/work/PhD_code/PhD_chapter1/data/CAFE_convergence/runs_to_test_convergence"


if __name__ == "__main__":

    out_data = "/Users/milena/work/PhD_code/PhD_chapter1/data/CAFE_convergence"
    cafe_outputs_list = CAFE_output_paths(CAFE_runs_dir)

    
    overlap_sig_OGs, all_OGs_list = get_overlap_OG_sig_list(cafe_outputs_list)
    print(f"{len(overlap_sig_OGs)} orthogroups significant out of {len(all_OGs_list)} in all {len(cafe_outputs_list)} CAFE runs")

    if True:
        ## plot gene family sizes
        orthogroups_orthoDB_filepath = "/Users/milena/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv"
        orthoDB_dict_lists = OGs.parse_orthogroups_dict(orthogroups_orthoDB_filepath)
        orthoDB_dict = OGs.get_GF_sizes(orthoDB_dict_lists)
        tree = "/Users/milena/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_only_species_names.nw"
        species_names = gff.make_species_order_from_tree(tree)
        plot_OG.plot_gene_counts(orthoDB_dict, overlap_sig_OGs, all_OGs_list, species_names, annotation = "orthoDB", 
                                 filename = f"{out_data}/significant_orthogroups_from_CAFE.png", 
                                 transparent_bg=True, 
                                 title = f"Orthogroups that are significant in all {len(cafe_outputs_list)} CAFE5 runs", 
                                 svg = False)

    if False:
        ## plot histogram of the number of runs that OGs are significant in
        sig_OG_counts = count_OG_occurence(cafe_outputs_list)
        sig_counts_list = list(sig_OG_counts.values())
        plot_sig_OG_occurence_histogram(sig_counts_list, f"{out_data}/runs_sig_OGs_hist.png")
    
