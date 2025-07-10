"""
Plot some stuff to evaluate the functional annotation of the gene groups
"""

import plot_basics as my_plotting
import parse_orthogroups as OGs
import plotting.plot_significant_orthogroups_from_CAFE as plot_OG
import matplotlib.pyplot as plt


def filepaths():
    orthogroups_orthoDB_filepath = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv"
    out_dir = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/"
    tree = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_only_species_names.nw"
    return out_dir,orthogroups_orthoDB_filepath,tree

def orthogroups_lists():
    out_dict = {
        "Aobt_expansion" : ["N0.HOG0000035","N0.HOG0000014"],
        "Gene Group 1" : ["N0.HOG0000027","N0.HOG0000059","N0.HOG0000095","N0.HOG0000204","N0.HOG0001077","N0.HOG0000140","N0.HOG0000492","N0.HOG0001030"],
        "Gene Group 3" : ["N0.HOG0000086","N0.HOG0002393","N0.HOG0000085"],
        "Gene Group 17" : ["N0.HOG0000669"],
        "Gene Group 5,8,9" : ["N0.HOG0000541","N0.HOG0000775","N0.HOG0000892","N0.HOG0000401","N0.HOG0009002"]
    }
    return out_dict

def plot_selected_OGs(orthogroups_path:str, OG_IDs:list, tree_path:str, filename:str, out_dir:str = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/", title = ""):
    """
    plot the gene family sizes in a subset of orthogroup_IDs
    """

    # get sorted species names from tree
    species_names_unsorted = my_plotting.plot_tree_manually(tree_path)
    species_coords_sorted = sorted(list(species_names_unsorted.keys()))
    species_list = [species_names_unsorted[species_coord] for species_coord in species_coords_sorted]

    orthoDB_dict_lists = OGs.parse_orthogroups_dict(orthogroups_path)
    orthoDB_dict = OGs.get_GF_sizes(orthoDB_dict_lists)
    OGs_of_interest_dict = {OG_id : orthoDB_dict[OG_id] for OG_id in OG_IDs}
    
    plot_OG.plot_gene_counts(OGs_of_interest_dict, sig_list=OG_IDs, all_cafe_list=OG_IDs, species_names=species_list, annotation = "orthoDB", filename = f"{out_dir}{filename}", transparent_bg=False, title=title)



if __name__ == "__main__":
    
    out_dir,orthogroups_orthoDB_filepath,tree_path = filepaths()
    OG_lists_dict = orthogroups_lists()

    # --> AOBT EXPANSION
    OGs_title = " and ".join(OG_lists_dict["Aobt_expansion"])
    # plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Aobt_expansion"], tree_path=tree_path, filename="Aobt_expansion_GF_sizes.png", title = f"A. obtectus expansion: {OGs_title}")

    # --> DETOXIFICATION
    # plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 1"], tree_path=tree_path, filename="Gene_Group_1_detoxofication_GF_sizes.png", title = "Gene group 1")
    # plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 17"], tree_path=tree_path, filename="Gene_Group_17_detoxofication_GF_sizes.png", title = "Gene group 17")
    # plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=["N0.HOG0000140"], tree_path=tree_path, filename="detoxofication_N0.HOG0000140_GF_sizes.png", title = "Orthogroup N0.HOG0000140")
    # plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 3"], tree_path=tree_path, filename="Gene_Group_3_lipid_metabolism_GF_sizes.png", title = "Gene Group 3")

    # --> REPRODUCTION
    # plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 5,8,9"], tree_path=tree_path, filename="Gene_Group_5_8_9_reproduction_GF_sizes.png", title = "Gene group 5, 8, and 9")
    # plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=["N0.HOG0000401"], tree_path=tree_path, filename="OG_N0.HOG0000401_reproduction_GF_sizes.png", title = "Orthogroup N0.HOG0000401")
    

