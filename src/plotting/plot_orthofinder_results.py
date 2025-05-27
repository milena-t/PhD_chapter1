import pandas as pd
from collections import Counter
from statistics import mean, stdev
import src.parse_gff as gff

import matplotlib.pyplot as plt
import os
import subprocess as sp

####!!
# on uppmax load biopython/1.80-py3.10.8 to use argparse!
import argparse



def parse_args():
    # Create the parser
    program_description = """
A script to plot the gene number and gene family/prthogroup number from up to three different orthofinder runs on the same set of species. 
the orthogroups input is the N0.tsv file from  the phylogenetically hierarchically orthogroups that orthofinder makes. The old Orthogroups directory 
and its contents is deprecated (some of the orthogroup IDs are duplicated for example)

The plot filenames are hardcoded into the end of the script!
"""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,description=program_description)

    # Add the arguments
    parser.add_argument('--orthogroups1', type=str, required=True, help='Filepath to the first set of hierarchical orthogroups')
    parser.add_argument('--orthogroups2', type=str, help="Filepath to the second set of hierarchical orthogroups")
    parser.add_argument('--orthogroups3', type=str, help="Filepath to the third set of hierarchical orthogroups")

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


    # Parse the arguments
    args = parser.parse_args()

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


def get_species_counts(orthogroups_file, species_names):
    print(f"\n{orthogroups_file}")
    
    orthogroups_dict = read_hierarchical_orthogroups_file(orthogroups_file)
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




#####################################################
###### Get Orthogroups per species
#####################################################



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




if __name__ == '__main__':
    

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


    species_names = gff.make_species_order_from_tree(args.common_tree)
    # print(species_names)


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