import pandas as pd
from collections import Counter
from statistics import mean, stdev
import re

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import seaborn as sns
import os
import subprocess as sp
from tqdm import tqdm

####!!
# on uppmax load biopython/1.80-py3.10.8 to use argparse!
import argparse



def parse_args():
    # Create the parser
    program_description = """
A script to plot the gene number and gene family/prthogroup number from up to three different orthofinder runs on the same set of species. 
the orthogroups input is the Orthogroups output directory that orthofinder makes. Specifically it uses these three files:
Orthogroups.GeneCount.tsv, Orthogroups.txt, and Orthogroups_UnassignedGenes.tsv.

The plot filenames are hardcoded into the end of the script!
"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description=program_description)

    # Add the arguments
    parser.add_argument('--orthogroups1', type=str, required=True, help='Filepath to the first set of hierarchical orthogroups. can also be a directory with the first set of orthogroups (required, default: native annotations)')
    parser.add_argument('--orthogroups2', type=str, help="Filepath to the second set of hierarchical orthogroups. can also be a directory with the second set of orthogroups (not required, default: orthoDB annotations)")
    parser.add_argument('--orthogroups3', type=str, help="Filepath to the third set of hierarchical orthogroups. can also be a directory with the third set of orthogroups (not required, default: species-specific proteinseqs annotations)")

    parser.add_argument('--common_tree', type=str, required=True, help='path to a newick tree containing all species (from an orthofinder run) to determine the plot order of species in the X axis')

    parser.add_argument('--annotations_og1', type=str, help='directory with annotations or proteinfasta files that are the basis for orthogroups set 1. Required for hierarchical orthogroups (file names need to contain species names!)')
    parser.add_argument('--annotations_og2', type=str, help="directory with annotations or proteinfasta files that are the basis for orthogroups set 2. Required for hierarchical orthogroups (file names need to contain species names!)")
    parser.add_argument('--annotations_og3', type=str, help="directory with annotations or proteinfasta files that are the basis for orthogroups set 3. Required for hierarchical orthogroups (file names need to contain species names!)")

    parser.add_argument('--name_og1', type=str, help='legend name for orthogroups set 1')
    parser.add_argument('--name_og2', type=str, help="legend name for orthogroups set 2")
    parser.add_argument('--name_og3', type=str, help="legend name for orthogroups set 3")

    parser.add_argument('--color_og1', type=str, help='line color for orthogroups set 1')
    parser.add_argument('--color_og2', type=str, help="line color for orthogroups set 2")
    parser.add_argument('--color_og3', type=str, help="line color for orthogroups set 3")

    parser.add_argument('--include_unassigned_genes', action=argparse.BooleanOptionalAction, help="include the unassigned genes in t he plot")

    # Parse the arguments
    args = parser.parse_args()

    if os.path.isdir(args.orthogroups1) and args.orthogroups1[-1]!="/": # if the user does not include a trailing "/", add it
        args.orthogroups1=args.orthogroups1+"/"
    if args.orthogroups2 and args.orthogroups2[-1]!="/": # if the user does not include a trailing "/", add it
        if os.path.isdir(args.orthogroups2):
            args.orthogroups2=args.orthogroups2+"/"
    if args.orthogroups3 and args.orthogroups3[-1]!="/": # if the user does not include a trailing "/", add it
        if os.path.isdir(args.orthogroups3):
            args.orthogroups3=args.orthogroups3+"/"

    # if args.species_names:
    #     args.species_names = args.species_names.split(",") # split the csv species names into a list

    if not args.name_og1:
        args.name_og1 = "native" 
    if not args.name_og2:
        args.name_og2 = "orthoDB" 
    if not args.name_og3:
        args.name_og3 = "proteinseqs"

    if not args.color_og1:
        args.color_og1 = "#B82845" # cardinal
    if not args.color_og2:
        args.color_og2 = "#4D7298" # UCLA blue
    if not args.color_og3:
        args.color_og3 = "#ED7D3A" # Pumbkin. Old default: "#413C58" # english violet

    return args



############################################
#### process the output of orthofinder #####
############################################


#### process the Orthogroups.txt output that orthofinder generates into numbers of gene families


def read_orthogroups_file(path_to_orthogroups_file):
    """
    makes a dictionary containing {orthogroup_ID : [list, of, orthogroup, member, gene_IDs]}
    """
    orthogroups_dict = {} 
    with open(path_to_orthogroups_file, "r") as orthogroups_file:
        for orthogroup in orthogroups_file:
            orthogroup = orthogroup.split(": ")
            orthogroups_dict[orthogroup[0]] = orthogroup[1].split(" ")
    return orthogroups_dict

def read_hierarchical_orthogroups_file(path_to_orthogroups_file):
    """
    makes a dictionary containing {orthogroup_ID : [list, of, orthogroup, member, gene_IDs]}
    """
    orthogroups_dict = {} 
    with open(path_to_orthogroups_file, "r") as orthogroups_file:
        next(orthogroups_file) # skip first line because it contains the file headers
        for orthogroup in orthogroups_file:
            og_list = orthogroup.split("\t")
            og_number = og_list[1]
            # hierarchical_og_number = og_list[0]
            # node_number = og_list[2]
            og_members = ", ".join([og_species.strip() for og_species in og_list[3:] if len(og_species.strip())>0])
            orthogroups_dict[og_number] = og_members.split(", ")

    return orthogroups_dict


def get_species_counts(orthogroups_file, species_names, hierarchical = False):
    print(f"\n{orthogroups_file}")
    if hierarchical:
        orthogroups_dict = read_hierarchical_orthogroups_file(orthogroups_file)
    else:    
        orthogroups_dict = read_orthogroups_file(orthogroups_file)
    print(str(len(orthogroups_dict))+" Orthogroups found")
    print("Mean orthogroup size: "+ str(mean([len(og) for og in orthogroups_dict.values()])))

    # count number of gene families per species
    gf_deduplicated = {}
    for OG_id in orthogroups_dict.keys():
        species_presence_list = []
        for species_name in species_names:
            if any(species_name in gene for gene in orthogroups_dict[OG_id]):
                species_presence_list.append(species_name)
        gf_deduplicated[OG_id] = species_presence_list

    gf_merged = [gf for gf in gf_deduplicated.values()] # this makes a list of sublists of all gene families
    gf_merged = [item for sublist in gf_merged for item in sublist] # make one long list without sublists
    species_counts = Counter(gf_merged)
    return(species_counts)

def split_at_second_occurrence(s, char): # split the gene string at the second occurence of "_" to get only the species name
    if s.count(char)<2:
        return s
    else:
        second_occurrence = s.find(char, 2) # start after the first occurence of "_"
        species = s[:second_occurrence]
        return species

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
        species_names = [split_at_second_occurrence(leaf, "_") for leaf in leaf_names]
    return species_names

#####################################################
###### Get Orthogroups per species
#####################################################

### This function is incomplete without the UnassignedGenes file since the GeneCount file does not contain single-gene orthogroups!


### modify the header columns of all the tsv files to exclude the file extension suffixes: 
# sed -i '' 's/_filtered_proteinfasta//g' orthoDB_braker_annotations_orthogroups_GeneCount.tsv
# sed -i '' 's/_filtered_proteinfasta//g' orthoDB_braker_annotations_orthogroups_UnassignedGenes.tsv
# sed -i '' 's/_transcripts//g' proteinref_braker_annotations_orthogroups_GeneCount.tsv
# sed -i '' 's/_transcripts//g' proteinref_braker_annotations_orthogroups_UnassignedGenes.tsv





def make_headers_lookup_table(genecount_table, species_names):
    """
    in case the header of the N0.tsv file or genecount table does not match the species names
    create a dictionary that associates the two
    """

    lookup_dict = {}
    genecount_header = list(genecount_table.columns.values)
    # this assumes that the species name is a substring of the corresponding table header (which is named after the input filenames, so if this script is used with the rest of the pipeline as intended it will work)
    for protein_set_name in genecount_header:
        try:
            species = [species for species in species_names if species in protein_set_name][0]
            lookup_dict[species] = protein_set_name
        except:
            lookup_dict[protein_set_name] = protein_set_name
    return lookup_dict


def get_gf_count_per_species(species_names, genecount_file, unassigned_file):
    # read GeneCounts file
    genecount_table = pd.read_csv(genecount_file, sep="\t", header=0)
    genecount_table.set_index("Orthogroup", inplace=True)

    header_lookup_detect = make_headers_lookup_table(genecount_table, species_names)

    gf_presence_table = genecount_table.map(lambda x: 1 if x > 0 else 0)
    num_orthogroups=gf_presence_table.shape[0]
    # read unassigned file
    with open(unassigned_file, "r") as unassigned:
        unassigned_lines = unassigned.readlines()[1:] # skip header line

    gf_num_dict = {}
    unassigned_count_dict = {}
    for species in species_names:
        species_header = header_lookup_detect[species]
        if species_header in gf_presence_table.columns: 
            gf_num_dict[species]=sum(gf_presence_table[species_header].tolist())
        else:
            gf_num_dict[species]=0
            print(f"{species} not found! ")
        unassigned_count_dict[species] = len([line for line in unassigned_lines if species in line])

    # gf_ordered_list = [gf_num_dict[species]+unassigned_count_dict[species] for species in species_names]
    return gf_num_dict, unassigned_count_dict

# returns the number of individual genes present in each species, 
# in a list ordered like the "species_names" vector

def get_gene_count_per_species(species_names, genecount_file, unassigned_file):
    # read GeneCounts file
    genecount_table = pd.read_csv(genecount_file, sep="\t", header=0)
    genecount_table.set_index("Orthogroup", inplace=True)

    header_lookup_detect = make_headers_lookup_table(genecount_table, species_names)

    num_orthogroups=genecount_table.shape[0]
    # read unassigned file
    with open(unassigned_file, "r") as unassigned:
        unassigned_lines = unassigned.readlines()[1:] # skip header line

    gene_num_dict = {}
    unassigned_count_dict = {}
    for species in species_names:
        species_header = header_lookup_detect[species]
        if species_header in genecount_table.columns:
            gene_num_dict[species]=sum(genecount_table[species_header].tolist())
        else:
            print(f"{species} not found!" )
            gene_num_dict[species]=0

        unassigned_count_dict[species] = len([line for line in unassigned_lines if species in line])
        # print(species)
        # print(unassigned_count_dict[species])

    #genes_ordered_list = [gf_num_dict[species]+unassigned_count_dict[species] for species in species_names]
    return gene_num_dict, unassigned_count_dict



def get_numbers_dicts(species_names, genecount_file, unassigned_file, verbose = False, include_unassigned = True):
    gf_num_dict, unassigned_dict = get_gf_count_per_species(species_names, genecount_file, unassigned_file)
    gene_num_dict, unassigned_dict = get_gene_count_per_species(species_names, genecount_file, unassigned_file)

    if verbose:
        print(f"genes: {gene_num_dict}, \ngene families: {gf_num_dict}, \nunassigned: {unassigned_dict}")
    
    if include_unassigned:
        numbers = {
            "Number of orthogroups" : gf_num_dict,
            "Number of genes" : gene_num_dict, 
            "Unassigned genes" : unassigned_dict
        }
    else:
        numbers = {
            "Number of orthogroups" : gf_num_dict,
            "Number of genes" : gene_num_dict, 
        }

    return(numbers)


def get_gene_conuts_from_annot(species_names, files_dir):
    files_list = os.listdir(files_dir)
    gene_nums_dict = {}
    search_string = "gene\t"
    if "fa" in files_list[0].split(".")[-1] or "fna" in files_list[0].split(".")[-1]:
        search_string = ">"
    for species_name in species_names: 
        filepath = files_dir +"/"+ [file for file in files_list if species_name in file][0]

        # get count per gene with bash commands 
        # in total: grep "search_string" filepath | wc -l

        command = ["grep", search_string, filepath]
        # print(command)
        grep_out = sp.run(command , stdout = sp.PIPE)
        command2 = ["wc", "-l"]
        num_hits = sp.run(command2 , stdout = sp.PIPE, input=grep_out.stdout).stdout.decode("utf-8")

        gene_nums_dict[species_name] = int(num_hits.strip())

    return(gene_nums_dict)


def get_gf_count_per_species_hierarchical(species_names, orthogroups_file):
    # read GeneCounts file
    gene_num_dict = {}

    for species in species_names:
        gene_num_dict[species] = 0


    print(f" ------------> {orthogroups_file}")
    gene_families = pd.read_csv(orthogroups_file, sep="\t")
    
    header_lookup_dict = make_headers_lookup_table(gene_families, species_names)
    # make the dataframe a binary matrix with just presence/absence data
    gene_families = gene_families.map(lambda x: 1 if pd.notna(x) else 0) # .notna().astype(int)
    
    for species in species_names:
        gene_num_dict[species] = gene_families[header_lookup_dict[species]].sum()

    return gene_num_dict


def get_hierarchical_numbers_dicts(species_names, files_dir, orthogroups_file, verbose = False):
    # gf_num_dict, unassigned_dict = get_gf_count_per_species(species_names, genecount_file, unassigned_file)
    gf_num_dict = get_gf_count_per_species_hierarchical(species_names, orthogroups_file)
    gene_num_dict = get_gene_conuts_from_annot(species_names, files_dir)
    # print(gene_num_dict)
    
    if verbose:
        print(f"genes: {gene_num_dict}, \ngene families: {gf_num_dict}")
    
    numbers = {
        "No. gene families" : gf_num_dict,
        "No. genes" : gene_num_dict, 
    }

    return(numbers)


#####################################################
#### compare orthofinder runs #######################
#####################################################

# plot function with lots of stuff to include optionally. The order of the species in the x-axis is determined by the "speciesnames" variable

def plot_general_annotation_comparisons(native = {}, orthoDB = {}, proteinseqs = {},legend_title = "", speciesnames = [], filename = "", genome_size=[], gs_measurement = []):

    fs = 13 # set font size
    # plot each column in the dataframe as a line in the same plot thorugh a for-loop
    fig = plt.figure(figsize=(10,9))
    
    ax = fig.add_subplot(1, 1, 1)
    
    ylab="Number of genes and orthogroups"
    # get a list of lists with [native, orthoDB] number of gene families per species

    legend_labels = []
    ncol_legend = 0

    if len(orthoDB)>0:
        legend_labels = list(orthoDB.keys())
        ncol_legend += 1
        # col = "cornflowerblue"
        col = args.color_og2
        for category in orthoDB:
            values = orthoDB[category]
            values = [values[species] for species in speciesnames]
            if category == legend_labels[0]: # genes
                ax.plot(speciesnames, values, linestyle=':', label = category+f" ({args.name_og2})", color = col)
            elif category == legend_labels[1]: # gene_families
                ax.plot(speciesnames, values, label = category+f" ({args.name_og2})", color = col)
            elif category == legend_labels[2]: # unassigned
                ax.plot(speciesnames, values, linestyle = (0,(5,10)), label = category+f" ({args.name_og2})", color = col)
            plt.yticks(fontsize = fs)
            # add space to see complete species names in the xticks
            plt.subplots_adjust(bottom=0.3)
        ax.set_ylabel(ylab, fontsize = fs)
    
    if len(proteinseqs)>0:
        legend_labels = list(proteinseqs.keys())
        ncol_legend += 1
        # col = "orchid"
        col = args.color_og3
        for category in proteinseqs:
            values = proteinseqs[category]
            values = [values[species] for species in speciesnames]
            if category == legend_labels[0]: # genes
                ax.plot(speciesnames, values, linestyle=':', label = category+f" ({args.name_og3})", color = col)
            elif category == legend_labels[1]: # gene_families
                ax.plot(speciesnames, values, label = category+f" ({args.name_og3})", color = col)
            elif category == legend_labels[2]: # unassigned
                ax.plot(speciesnames, values, linestyle = (0,(5,10)), label = category+f" ({args.name_og3})", color = col)
            plt.yticks(fontsize = fs)
            # add space to see complete species names in the xticks
            plt.subplots_adjust(bottom=0.3)
        ax.set_ylabel(ylab, fontsize = fs)

    if len(native)>0:
        legend_labels = list(native.keys())
        ncol_legend += 1
        col = args.color_og1
        for category in native: # categories are: ["genes", "gene families"]
            values = native[category]
            values = [values[species] for species in speciesnames]
            annotation_method = f" ({args.name_og1})"
            # annotation_method = " (Kaufmann2023)"
            if category == legend_labels[0]: # genes
                ax.plot(speciesnames, values, linestyle=':', label = category+annotation_method, color = col)
            elif category == legend_labels[1]: # gene_families
                ax.plot(speciesnames, values, label = category+annotation_method, color = col)
            elif category == legend_labels[2]: # unassigned
                ax.plot(speciesnames, values, linestyle = (0,(5,10)), label = category+annotation_method, color = col)
            plt.yticks(fontsize = fs)
            # add space to see complete species names in the xticks
            plt.subplots_adjust(bottom=0.3)

        ax.set_ylabel(ylab, fontsize = fs)
        plt.xticks(labels=[species.replace("_", ". ") for species in speciesnames], ticks=speciesnames, rotation = 90, fontsize = fs)

    if len(genome_size)>0 and len(native)>1:
        genome_size = [genome_size[species] for species in species_names]
        col_emp = "mediumseagreen"
        col_ass = "darkcyan"
        ax2 = ax.twinx() 
        # ax2.spines['right'].set_visible(False)
        if len(gs_measurement)>0:
            col_genome_size =[col_emp if col == 1 else col_ass for col in gs_measurement] # TODO this doesn't work yet
        else:
            col_genome_size =col_emp
        ax2.set_ylabel('Genome size in Mb', color = col_genome_size, fontsize = fs)
        ax2.tick_params(axis ='y', labelcolor = col_genome_size, labelsize = fs)  
        plot_2 = ax2.scatter(speciesnames, genome_size, color = col_genome_size)

    # include legend (reversed order)
    # handles, labels = ax.get_legend_handles_labels()
    ax.legend(title=legend_title, loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol = ncol_legend)
    # set grid only for X axis ticks 
    ax.grid(True)
    ax.yaxis.grid(False)

    # save plot
    #ax.set_ylabel(y_label, fontsize = fs)
    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)
    # plt.show()



orthofinder_results_dict = "/Users/miltr339/work/orthofinder/run_with_topological_tree/"

native_orthogroups = orthofinder_results_dict+"native_annotations_orthogroups.txt"
native_genecount = orthofinder_results_dict+"native_annotations_orthogroups_GeneCount.tsv"
native_unassigned = orthofinder_results_dict+"native_annotations_orthogroups_UnassignedGenes.tsv"
orthoDB_orthogroups = orthofinder_results_dict+"orthoDB_braker_annotations_orthogroups.txt"
orthoDB_genecount = orthofinder_results_dict+"orthoDB_braker_annotations_orthogroups_GeneCount.tsv"
orthoDB_unassigned = orthofinder_results_dict+"orthoDB_braker_annotations_orthogroups_UnassignedGenes.tsv"
proteinref_orthogroups = orthofinder_results_dict+"proteinref_braker_annotations_orthogroups.txt"
proteinref_genecount = orthofinder_results_dict+"proteinref_braker_annotations_orthogroups_GeneCount.tsv"
proteinref_unassigned = orthofinder_results_dict+"proteinref_braker_annotations_orthogroups_UnassignedGenes.tsv"

genome_sizes = [180, 842, 471, 399, 250, 204, 258, 461, 589, 223, 949, 375, 701, 971, 1202] # the last three are in the order of C_chinensis(C_analis,C_maculatus) 
#    old species order with (C_chinensis,C_analis)C_maculatus
#    species_names = ["D_melanogaster","I_luminosus","P_pyralis","C_septempunctata","A_verrucosus","T_castaneum","T_molitor","Z_morio","R_ferrugineus","D_ponderosae","A_obtectus","B_siliquastri","C_chinensis","C_analis","C_maculatus"]

if False:

    # old data! left for archiving reasons
    native_annotations = {
        "genes" : {'D_melanogaster': 11218, 'I_luminosus': 22762, 'P_pyralis': 14671, 'C_septempunctata': 14114, 'A_verrucosus': 11027, 'T_castaneum': 11939, 'T_molitor': 17508, 'Z_morio': 24017, 'R_ferrugineus': 15891, 'D_ponderosae': 12181, 'A_obtectus': 30444, 'B_siliquastri': 15426, 'C_chinensis': 29579, 'C_analis': 28288, 'C_maculatus': 33238}, 
        'gene families' : {'D_melanogaster': 7841, 'I_luminosus': 10857, 'P_pyralis': 9346, 'C_septempunctata': 9463, 'A_verrucosus': 7764, 'T_castaneum': 9379, 'T_molitor': 10371, 'Z_morio': 11157, 'R_ferrugineus': 10628, 'D_ponderosae': 8999, 'A_obtectus': 12602, 'B_siliquastri': 10744, 'C_chinensis': 12472, 'C_analis': 12632, 'C_maculatus': 13643}, 
        'unassigned' : {'D_melanogaster': 2308, 'I_luminosus': 1723, 'P_pyralis': 642, 'C_septempunctata': 595, 'A_verrucosus': 713, 'T_castaneum': 166, 'T_molitor': 2258, 'Z_morio': 4312, 'R_ferrugineus': 5590, 'D_ponderosae': 959, 'A_obtectus': 4260, 'B_siliquastri': 2491, 'C_chinensis': 4129, 'C_analis': 3327, 'C_maculatus': 2590}
    }
    orthoDB_annotations = {
        'genes' : {'D_melanogaster': 13821, 'I_luminosus': 48605, 'P_pyralis': 37903, 'C_septempunctata': 36363, 'A_verrucosus': 17554, 'T_castaneum': 17436, 'T_molitor': 19785, 'Z_morio': 29791, 'R_ferrugineus': 18821, 'D_ponderosae': 16622, 'A_obtectus': 55137, 'B_siliquastri': 16616, 'C_chinensis': 29571, 'C_analis': 24450, 'C_maculatus': 39133}, 
        'gene families' : {'D_melanogaster': 8453, 'I_luminosus': 15621, 'P_pyralis': 14333, 'C_septempunctata': 13359, 'A_verrucosus': 11662, 'T_castaneum': 11536, 'T_molitor': 11314, 'Z_morio': 13484, 'R_ferrugineus': 11074, 'D_ponderosae': 10911, 'A_obtectus': 15421, 'B_siliquastri': 11541, 'C_chinensis': 13293, 'C_analis': 12514, 'C_maculatus': 15595}, 
        'unassigned' : {'D_melanogaster': 3572, 'I_luminosus': 7891, 'P_pyralis': 3277, 'C_septempunctata': 4373, 'A_verrucosus': 0, 'T_castaneum': 2003, 'T_molitor': 3173, 'Z_morio': 6829, 'R_ferrugineus': 5620, 'D_ponderosae': 1984, 'A_obtectus': 5377, 'B_siliquastri': 2057, 'C_chinensis': 2145, 'C_analis': 1230, 'C_maculatus': 3186}
    }
    proteinref_annotations = {
        'genes' : {'D_melanogaster': 13894, 'I_luminosus': 62847, 'P_pyralis': 39283, 'C_septempunctata': 36107, 'A_verrucosus': 22989, 'T_castaneum': 17899, 'T_molitor': 20287, 'Z_morio': 35531, 'R_ferrugineus': 21262, 'D_ponderosae': 17928, 'A_obtectus': 60296, 'B_siliquastri': 18224, 'C_chinensis': 35171, 'C_analis': 29625, 'C_maculatus': 47381}, 
        'gene families' : {'D_melanogaster': 8477, 'I_luminosus': 18267, 'P_pyralis': 14864, 'C_septempunctata': 13552, 'A_verrucosus': 13631, 'T_castaneum': 11789, 'T_molitor': 11579, 'Z_morio': 14743, 'R_ferrugineus': 11992, 'D_ponderosae': 11400, 'A_obtectus': 16610, 'B_siliquastri': 12277, 'C_chinensis': 14734, 'C_analis': 14225, 'C_maculatus': 17515}, 
        'unassigned' : {'D_melanogaster': 3682, 'I_luminosus': 10615, 'P_pyralis': 3420, 'C_septempunctata': 4178, 'A_verrucosus': 0, 'T_castaneum': 1891, 'T_molitor': 3021, 'Z_morio': 7497, 'R_ferrugineus': 6949, 'D_ponderosae': 2304, 'A_obtectus': 5929, 'B_siliquastri': 2460, 'C_chinensis': 3402, 'C_analis': 1651, 'C_maculatus': 3857}
    }

    Kaufman_cmac_annotation_comparison = {
        'genes' : {'A_obtectus' : 34704, 'c_maculatus_all_proteinrefs' : 51238, 'C_maculatus__' : 35865, 'c_maculatus_only_orthoDB' : 42465, 'c_maculatus_RNA_combined' : 23428, 'c_maculatus_RNA_simple' : 23226, 'T_castaneum' : 12170},
        'gene families' : {'c_maculatus_all_proteinrefs': 31320, 'c_maculatus_only_orthoDB': 30075, 'C_maculatus__': 24597, 'c_maculatus_RNA_combined': 17595, 'c_maculatus_RNA_simple': 17480, 'A_obtectus': 16299, 'T_castaneum': 11210},
        'unassigned' : {'T_castaneum': 10225, 'C_maculatus__': 2141, 'c_maculatus_only_orthoDB': 2800, 'c_maculatus_all_proteinrefs': 2094, 'c_maculatus_RNA_simple': 138, 'c_maculatus_RNA_combined': 196, 'A_obtectus': 7114}
    }

    # new species order with (C_maculatus,C_analis)C_chinensis
    # species_names = ["D_melanogaster","I_luminosus","P_pyralis","C_septempunctata","A_verrucosus","T_castaneum","T_molitor","Z_morio","R_ferrugineus","D_ponderosae","A_obtectus","B_siliquastri","C_chinensis","C_analis","C_maculatus"]
    # plot_general_annotation_comparisons(native = native_annotations, speciesnames = species_names, filename = "annotation_native_with_reference_tree.png", genome_size=genome_sizes)
    # plot_general_annotation_comparisons(native = native_annotations, orthoDB = orthoDB_annotations, proteinseqs = proteinref_annotations, speciesnames = species_names, filename = "annotation_comparison_with_reference_tree.png", genome_size=genome_sizes)
    # plot_general_annotation_comparisons(native = native_annotations, orthoDB = orthoDB_annotations, proteinseqs = proteinref_annotations, speciesnames = species_names, filename = "annotation_comparison_with_reference_tree.png")

    # species_names = ["T_castaneum", "C_maculatus__", "c_maculatus_only_orthoDB", "c_maculatus_all_proteinrefs", "c_maculatus_RNA_simple", "c_maculatus_RNA_combined", "A_obtectus"]
    # plot_general_annotation_comparisons(native = Kaufman_cmac_annotation_comparison, speciesnames = species_names, filename = "Kaufmann_annotation_comparison.png")
    pass



if __name__ == '__main__':
    
    # example command line input:
    # python3 06b_orthofinder_plotting.py  --orthogroups1 /Users/miltr339/work/orthofinder/native_orthogroups --name_og1 "native" --orthogroups2 /Users/miltr339/work/orthofinder/orthoDB_orthogroups --name_og2 "orthoDB"  --orthogroups3 /Users/miltr339/work/orthofinder/orthoDB_orthogroups_uniform_repeatmasking  --name_og3 "orthoDB repeatmasked" --common_tree /Users/miltr339/Box\ Sync/code/annotation_pipeline/14_species_orthofinder_tree.nw
    # python3 06b_orthofinder_plotting.py  --orthogroups1 /Users/miltr339/work/orthofinder/native_orthogroups/hierarchical/N0.tsv --name_og1 "native" --annotations_og1 /Users/miltr339/work/native_annotations/all_native_annot --orthogroups2 /Users/miltr339/work/orthofinder/orthoDB_orthogroups/hierarchical/N0.tsv --name_og2 "orthoDB" --annotations_og2 /Users/miltr339/work/orthoDB_proteinseqs --orthogroups3 /Users/miltr339/work/orthofinder/orthoDB_orthogroups_uniform_repeatmasking/hierarchical/N0.tsv  --name_og3 "orthoDB repeatmasked" --annotations_og3 /Users/miltr339/work/orthoDB_annotations --common_tree /Users/miltr339/Box\ Sync/code/annotation_pipeline/14_species_orthofinder_tree.nw
    # python3 06b_orthofinder_plotting.py  --orthogroups1 /Users/miltr339/work/orthofinder/native_orthogroups/hierarchical/N0.tsv --name_og1 "native" --annotations_og1 /Users/miltr339/work/native_annotations/all_native_annot --common_tree /Users/miltr339/Box\ Sync/code/annotation_pipeline/14_species_orthofinder_tree.nw
    
    # no blue line for unmasked:
    # python3 06b_orthofinder_plotting.py  --orthogroups1 /Users/miltr339/work/orthofinder/native_orthogroups --name_og1 "native" --color_og2 "#F2933A" --orthogroups2 /Users/miltr339/work/orthofinder/orthoDB_orthogroups_uniform_repeatmasking  --name_og2 "orthoDB repeatmasked" --common_tree /Users/miltr339/Box\ Sync/code/annotation_pipeline/14_species_orthofinder_tree.nw

    args=parse_args()

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


    species_names = make_species_order_from_tree(args.common_tree)
    # print(species_names)

    ## check if the input is a directory 
    if os.path.isdir(args.orthogroups1):
        print(f" --> entered a directory, use non-hierarchical orthogroups.\n     (Hierarchical orthogroups are recommended)\n")

        # set default file names according to orthofinder output
        orthogroups_file = "Orthogroups.txt"
        genecount_file = "Orthogroups.GeneCount.tsv"
        unassigned_file = "Orthogroups_UnassignedGenes.tsv"

        # basic orthogroups statistics
        get_species_counts(args.orthogroups1+orthogroups_file, species_names)
        # first orthogroups set
        print("\nnative orthogroups set:")
        first_numbers = get_numbers_dicts(species_names, args.orthogroups1+genecount_file, args.orthogroups1+unassigned_file, verbose = True, include_unassigned=args.include_unassigned_genes)
        # native tree:
        # (D_melanogaster:0.331686,((P_pyralis:0.206536,I_luminosus:0.190788)0.833694:0.0924935,(C_septempunctata:0.280821,(((D_ponderosae:0.261936,R_ferrugineus:0.196725)0.801546:0.0855842,(A_obtectus:0.0824289,(B_siliquastri:0.0691699,(C_chinensis:0.0758605,C_maculatus:0.078485)0.554869:0.039176)0.363833:0.0292653)0.833076:0.156195)0.372488:0.0366516,(A_verrucosus:0.162389,(T_castaneum:0.0990607,(T_molitor:0.0922383,Z_morio:0.121058)0.253478:0.0259247)0.332612:0.0389913)0.759815:0.0928272)0.186399:0.0231684)0.309737:0.0466942)1:0.331686);

        
        if args.orthogroups2 and len(args.orthogroups2)>0:
            get_species_counts(args.orthogroups2+orthogroups_file, species_names)
            print("\northoDB annotation orthogroups set")
            second_numbers = get_numbers_dicts(species_names, args.orthogroups2+genecount_file, args.orthogroups2+unassigned_file, verbose = True, include_unassigned=args.include_unassigned_genes)
            # orthoDB tree
            # (D_melanogaster:0.348432,(((((B_siliquastri:0.0627412,(C_chinensis:0.0554959,C_maculatus:0.0607227)0.647242:0.0354391)0.45952:0.0373316,A_obtectus:0.10178)0.852536:0.169791,(R_ferrugineus:0.19925,D_ponderosae:0.230231)0.794484:0.0879981)0.367215:0.036785,(C_septempunctata:0.296374,(A_verrucosus:0.144991,((Z_morio:0.124737,T_molitor:0.0890143)0.238657:0.0226841,T_castaneum:0.100365)0.305383:0.0426554)0.751112:0.0958722)0.171263:0.0251832)0.294484:0.0484742,(P_pyralis:0.210399,I_luminosus:0.194079)0.806495:0.0973794)1:0.348432);
            # differently ordered but same topology as the native tree!
        else:
            second_numbers = {}

        if args.orthogroups3 and len(args.orthogroups3)>0:
            get_species_counts(args.orthogroups3+orthogroups_file, species_names)
            print("\northoDB filtered annotation orthogroups set")
            third_numbers = get_numbers_dicts(species_names, args.orthogroups3+genecount_file, args.orthogroups3+unassigned_file, verbose = True, include_unassigned=args.include_unassigned_genes)
        else:
            third_numbers = {}

        plot_general_annotation_comparisons(native = first_numbers, orthoDB = second_numbers, proteinseqs = third_numbers, legend_title = "", speciesnames = species_names, filename = "plot_orthofinder_comparison_uniform_repeats.png")
        plot_general_annotation_comparisons(native = first_numbers, orthoDB = second_numbers, proteinseqs = third_numbers, legend_title = "", speciesnames = species_names, filename = "plot_orthofinder_comparison_with_genome_sizes_uniform_repeats.png", genome_size = genome_sizes_dict)



    else: # input path is not a directory and must be the N0.tsv hierarchical orthogroups file instead
        print(f"\n --> entered a file, use hierarchical orthogroups.\n")


        get_species_counts(args.orthogroups1, species_names, hierarchical=True)
        # first orthogroups set
        print(f"\n {args.name_og1} orthogroups set:")
        first_numbers = get_hierarchical_numbers_dicts(species_names, files_dir= args.annotations_og1, orthogroups_file=args.orthogroups1, verbose = True)
        # native tree:
        # (D_melanogaster:0.331686,((P_pyralis:0.206536,I_luminosus:0.190788)0.833694:0.0924935,(C_septempunctata:0.280821,(((D_ponderosae:0.261936,R_ferrugineus:0.196725)0.801546:0.0855842,(A_obtectus:0.0824289,(B_siliquastri:0.0691699,(C_chinensis:0.0758605,C_maculatus:0.078485)0.554869:0.039176)0.363833:0.0292653)0.833076:0.156195)0.372488:0.0366516,(A_verrucosus:0.162389,(T_castaneum:0.0990607,(T_molitor:0.0922383,Z_morio:0.121058)0.253478:0.0259247)0.332612:0.0389913)0.759815:0.0928272)0.186399:0.0231684)0.309737:0.0466942)1:0.331686);

        
        if args.orthogroups2 and len(args.orthogroups2)>0:
            get_species_counts(args.orthogroups2, species_names, hierarchical=True)
            print(f"\n {args.name_og2} annotation orthogroups set")
            second_numbers = get_hierarchical_numbers_dicts(species_names, files_dir= args.annotations_og2, orthogroups_file=args.orthogroups2, verbose = True)
            # orthoDB tree
            # (D_melanogaster:0.348432,(((((B_siliquastri:0.0627412,(C_chinensis:0.0554959,C_maculatus:0.0607227)0.647242:0.0354391)0.45952:0.0373316,A_obtectus:0.10178)0.852536:0.169791,(R_ferrugineus:0.19925,D_ponderosae:0.230231)0.794484:0.0879981)0.367215:0.036785,(C_septempunctata:0.296374,(A_verrucosus:0.144991,((Z_morio:0.124737,T_molitor:0.0890143)0.238657:0.0226841,T_castaneum:0.100365)0.305383:0.0426554)0.751112:0.0958722)0.171263:0.0251832)0.294484:0.0484742,(P_pyralis:0.210399,I_luminosus:0.194079)0.806495:0.0973794)1:0.348432);
            # differently ordered but same topology as the native tree!
        else:
            second_numbers = {}

        if args.orthogroups3 and len(args.orthogroups3)>0:
            get_species_counts(args.orthogroups3, species_names, hierarchical=True)
            print(f"\n {args.name_og3} annotation orthogroups set")
            third_numbers = get_hierarchical_numbers_dicts(species_names, files_dir= args.annotations_og3, orthogroups_file=args.orthogroups3, verbose = True)
        else:
            third_numbers = {}


        print(f"\n\n ... plotting ...")
        plot_general_annotation_comparisons(native = first_numbers, orthoDB = second_numbers, proteinseqs = third_numbers, legend_title = "", speciesnames = species_names, filename = "plot_orthofinder_comparison.png")
        plot_general_annotation_comparisons(native = first_numbers, orthoDB = second_numbers, proteinseqs = third_numbers, legend_title = "", speciesnames = species_names, filename = "plot_orthofinder_comparison_with_genome_sizes.png", genome_size = genome_sizes_dict)





    
    # plot_general_annotation_comparisons(native = first_numbers, legend_title = "", speciesnames = species_names, filename = "plot_only_native_annotation.png", genome_size = genome_sizes_dict)