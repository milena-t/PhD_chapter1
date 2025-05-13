# plot basic annotation stats to compare the two annotation methods
#  * gene numbers phylogenetically sorted
#  * length distribution histogram
#  * number of introns histogram
#  * number of repeat motifs in intron regions
#  * number/proportion of single exon genes (and how many of them are single-exon in both annotations? check the 1-to-1 overlap)

import subprocess as sp
import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import parse_gff as gff 
# import src.parse_gff as gff 
from Bio import SeqIO, Phylo


def get_gene_conuts_from_annot(species_names, files_dir):
    files_list = os.listdir(files_dir)
    gene_nums_dict = {}
    gene_nums_dict_from_annot = {}
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

        ### if you want to compare the counting with grep to parse_gff:
        # if search_string == "gene\t":
        #     annot_dir = gff.parse_gff3_general(filepath=filepath, only_genes=True, verbose=False)
        #     gene_nums_dict_from_annot[species_name] = len(annot_dir)
        # if gene_nums_dict[species_name] != gene_nums_dict_from_annot[species_name]:
        #     print(f"{files_dir}/{species_name} -->  grep: {gene_nums_dict[species_name]} ,  parse_gff: {gene_nums_dict_from_annot[species_name]}")

    return(gene_nums_dict)


def get_coords(tree):
    """
    Traverse tree and assign x/y positions for each clade.
    The tree orientation is intended with leaves pointing upwards in a normal x/y coordinate system
    """
    coords = {}
    x = 0
    max_depth = 0 # get position of the deepest leaf so that you can align all the leaf labels later
    leaf_names = []
    
    # recursion through every node of the tree
    def assign(clade, depth):
        nonlocal x
        nonlocal max_depth
        
        if depth>max_depth:
            max_depth = depth
        if clade.is_terminal():
            coords[clade] = (x, depth) # depth is the distance from the root
            x += 1 # orientation of the node along the x-axis (left/right orientation)
        else:
            for child in clade.clades:
                assign(child, depth + clade.branch_length if clade.branch_length else depth)
            # if the node has multiple children, determine the x coordinates of all the children and take the mean to center the line in the middle above the children
            child_coords = [coords[c][0] for c in clade.clades]
            coords[clade] = (sum(child_coords) / len(child_coords), depth)

    assign(tree.root, 0)

    return coords, max_depth



def plot_tree_manually(species_tree, ax_tree, add_leaf_label=False):
    """
    plot a phylogenetic tree (newick format in a file in species_tree)
    manually so that the leaves point upwards. 
    ax_tree is the plot axis defined in fig, ax_tree = plt.subplots()
    """
    tree = Phylo.read(species_tree, "newick")
    tree.ladderize()
    coords, max_depth = get_coords(tree)
    leaf_names = {}
    # Draw manually
    for clade in tree.find_clades(order="level"):
        x, y = coords[clade]

        # Draw a horizontal line from this node to its children
        for child in clade.clades:
            x2, y2 = coords[child]

            # Vertical line
            ax_tree.plot([x2, x2], [y, y2], color="black")
            # Horizontal line
            ax_tree.plot([x, x2], [y, y], color="black")

        # Draw labels on tips
        if clade.is_terminal():
            ax_tree.plot([x, x], [y, y + max_depth-y + 0.2], color="black")
            if add_leaf_label:
                ax_tree.text(x, y + max_depth-y + 0.25, clade.name, ha="center", va="bottom", rotation=90)
            leaf_names[x] = clade.name

    # Adjust and flip y-axis so tree grows upwards
    ax_tree.set_ylim(ax_tree.get_ylim()[::-1])  # Invert Y axis
    ax_tree.invert_yaxis()
    ax_tree.axis("off")
    return leaf_names


def plot_gene_counts(native_annot_dir, species_tree, orthoDB_annot_dir="", orthoDB_filtered_annot_dir="", filename = "only_genome_sizes_14_species.png"):
    """
    plot gene counts from annotations (or proteinfasta, but preferably annotation), with a species tree on the x-axis
    """
    
    species_names = gff.make_species_order_from_tree(species_tree)
    
    
    fs = 15 # set font size
    # plot each column in the dataframe as a line in the same plot thorugh a for-loop
    # fig = plt.figure(figsize=(10,8))
    # ax = fig.add_subplot(1, 1, 1)

    fig, (ax_data, ax_tree) = plt.subplots(2, 1, figsize=(10, 15), gridspec_kw={'height_ratios': [1, 2]})
    
    
    species_names_unsorted = plot_tree_manually(species_tree, ax_tree)
    print(species_names_unsorted)
    species_coords_sorted = sorted(list(species_names_unsorted.keys()))
    species_names = [species_names_unsorted[species_coord] for species_coord in species_coords_sorted]

    ylab="Number of annotated genes"
    # get a list of lists with [native, orthoDB] number of gene families per species

    legend_labels = []
    ncol_legend = 0

    native_gene_nos = get_gene_conuts_from_annot(species_names, native_annot_dir)

    native_gene_list = [native_gene_nos[species] for species in species_names]
    ax_data.plot(species_names, native_gene_list, label = "native annotation", color = "#b82946") # red

    if len(orthoDB_filtered_annot_dir)>0:
        orthoDB_filtered_gene_nos = get_gene_conuts_from_annot(species_names, orthoDB_filtered_annot_dir)
        orthoDB_filtered_gene_list = [orthoDB_filtered_gene_nos[species] for species in species_names]
        ax_data.plot(species_names, orthoDB_filtered_gene_list, label = "orthoDB TE-filtered", color = "#F2933A") # orange
        ymax = max(native_gene_list+orthoDB_filtered_gene_list)*1.1
    
    if len(orthoDB_annot_dir)>0:
        orthoDB_gene_nos = get_gene_conuts_from_annot(species_names, orthoDB_annot_dir)
        orthoDB_gene_list = [orthoDB_gene_nos[species] for species in species_names]
        ax_data.plot(species_names, orthoDB_gene_list, label = "orthoDB uniform annotation", color = "#4d7298") # blue
        ymax = max(native_gene_list+orthoDB_gene_list)*1.1

    if len(orthoDB_annot_dir)>0 and len(orthoDB_filtered_annot_dir)>0:
        ymax = max(native_gene_list+orthoDB_gene_list+orthoDB_filtered_gene_list)*1.1

    ax_data.set_ylabel(ylab, fontsize = fs)
    # plt.xticks(ticks=range(len(species_names)), labels=[species.replace("_", ". ") for species in species_names], rotation = 90, fontsize = fs)
    ax_data.set_xticklabels([species.replace("_", ". ") for species in species_names], rotation=90, fontsize=fs)
    
    legend = ax_data.legend(fontsize = fs)
    # rotate legend text by 90 degrees (but it looks like shit)
    # legend = ax.legend(fontsize = fs, ncol=2)
    # for text in legend.get_texts():
    #     text.set_rotation(90)

    # set grid only for X axis ticks 
    ax_data.grid(True)
    ax_data.yaxis.grid(False)
    
    ax_data.set_ylim(5e3,ymax)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)


def plot_Kaufmann_annotation_comparison(Kaufman_cmac_annotation_comparison, Kaufmann_labels, filename = "Kaufmann_annotation_comparison.png"):
    fs = 15 # set font size
    # plot each column in the dataframe as a line in the same plot thorugh a for-loop
    fig = plt.figure(figsize=(10,8))
    
    ax = fig.add_subplot(1, 1, 1)
    
    lines = ["genes", "gene families"]

    ylab="Number of genes or gene families"
    # get a list of lists with [native, orthoDB] number of gene families per species

    legend_labels = []
    ncol_legend = 0

    categories = [key for key in Kaufmann_labels.keys()]
    x_ticks = [Kaufmann_labels[key] for key in categories]
    genes_numbers = [Kaufman_cmac_annotation_comparison[lines[0]][key] for key in categories]
    gene_family_numbers = [Kaufman_cmac_annotation_comparison[lines[1]][key] for key in categories]

    ax.plot(x_ticks, genes_numbers, label = "number of genes", color = "#b82946", linestyle = ":") # red dotted
    ax.plot(x_ticks, gene_family_numbers, label = "number of gene families", color = "#b82946") # red

    ax.set_ylabel(ylab, fontsize = fs)
    plt.xticks(labels=x_ticks, ticks=x_ticks, rotation = 90, fontsize = fs)
    
    ax.legend(fontsize = fs)
    # set grid only for X axis ticks 
    ax.grid(True)
    ax.yaxis.grid(False)

    ymax = max(genes_numbers)*1.1
    ymin = min(gene_family_numbers)*0.9
    ax.set_ylim(ymin,ymax)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)




def get_lengths_list(fasta_filepath):
    seq_lengths = []
    for record in SeqIO.parse(fasta_filepath, "fasta"):
        seq_lengths.append(len(record.seq))
    return seq_lengths


def plot_histogram_protein_lengths(native_path:str, orthoDB_path:str, species_name:str, no_bins = 20, max_length = 1000, filename = "protein_lengths_histogram.png"):
    """
    plot a histogram of the length distribution of the proteins in the input fasta files
    """
    native_lengths = [length for length in get_lengths_list(native_path) if length < max_length]
    orthoDB_lengths = [length for length in get_lengths_list(orthoDB_path) if length < max_length]

    fig, ax = plt.subplots(1,1, figsize=(15, 12))
    fs = 20
    plt.rcParams.update({'font.size': fs})

    colors = {
        "orthoDB" : "#F2933A",
        "native" : "#b82946"
    }
    plt.hist([native_lengths, orthoDB_lengths], bins=no_bins, histtype="bar", color = [colors["native"], colors["orthoDB"]], label=["native", "orthoDB"])
    
    ax.tick_params(axis='both', labelsize=fs)
    plt.xlabel(f"protein length (Aminoacids, up to {max_length})", fontsize = fs)
    plt.ylabel("", fontsize = fs)
    plt.legend(fontsize = fs)
    species_title = species_name.replace("_", ". ")
    plt.title(f"{species_title} protein length distribution")
    
    plt.savefig(filename, dpi = 300, transparent = True)
    print(f"plot saved in the current working directory as: {filename}")


def plot_all_species_protein_length_distribution(native_files:dict, orthoDB_files:dict, columns = 3, no_bins = 20, max_length = 1000, filename = "protein_lengths_histogram.png", legend_in_last = True):
    """
    plot a grid of histograms for all species. input is dictionaries with species names as keys and filepaths to aminoacid fasta files as values
    """
    cols = columns
    rows = int(len(native_files)/cols)  +1
    fig, axes = plt.subplots(rows, cols, figsize=(12, 15))
    fs = 15

    colors = {
        "orthoDB" : "#F2933A",
        "native" : "#b82946"
    }

    for idx, species in enumerate(native_files.keys()):
        
        native_lengths = [length for length in get_lengths_list(native_files[species]) if length < max_length]
        orthoDB_lengths = [length for length in get_lengths_list(orthoDB_files[species]) if length < max_length]
        
        # Calculate row and column indices for the current subplot
        row = idx // cols
        col = idx % cols
        species_name = species.replace("_", ". ")
        print(f"\tin position {row+1},{col+1}: \t{species_name}")

        # Plot histogram on the corresponding subplot axis
        axes[row, col].hist([native_lengths, orthoDB_lengths], bins=no_bins, histtype="bar", color = [colors["native"], colors["orthoDB"]])
        axes[row, col].set_title(f'{species_name}')
        axes[row, col].set_xlabel('')
        axes[row, col].set_ylabel('')
    
    # add legend to last plot square
    if legend_in_last == True:
        idx_max = len(native_files.keys())
        row = idx_max // cols
        col = idx_max % cols
        axes[row, col].axis('off')
        axes[row, col].set_title(f'{species_name}')
        handles = []
        labels = [] 
        handles.append(mpatches.Patch(color=colors["native"]))
        labels.append("native")
        handles.append(mpatches.Patch(color=colors["orthoDB"]))
        labels.append("orthoDB")
        axes[row, col].legend(handles, labels, fontsize = fs, loc='center', title_fontsize = fs)
    
    # Set a single x-axis label for all subplots
    x_label = f"protein length (Aminoacids, up to {max_length})"
    fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=12)
    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    plt.savefig(filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {filename}")



if __name__ == "__main__":

    ## plot basic gene counts
    if True:

        try:
            tree = "/Users/miltr339/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw"
            species_names = gff.make_species_order_from_tree(tree)
        except:
            tree = "/Users/milena/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw"
            species_names = gff.make_species_order_from_tree(tree)

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


        Kaufman_cmac_annotation_comparison = {
            'genes' : {'A_obtectus' : 34704, 'c_maculatus_all_proteinrefs' : 51238, 'C_maculatus__' : 35865, 'c_maculatus_only_orthoDB' : 42465, 'c_maculatus_RNA_combined' : 23428, 'c_maculatus_RNA_simple' : 23226, 'T_castaneum' : 12170},
            'gene families' : {'c_maculatus_all_proteinrefs': 31320, 'c_maculatus_only_orthoDB': 30075, 'C_maculatus__': 24597, 'c_maculatus_RNA_combined': 17595, 'c_maculatus_RNA_simple': 17480, 'A_obtectus': 16299, 'T_castaneum': 11210},
            'unassigned' : {'T_castaneum': 10225, 'C_maculatus__': 2141, 'c_maculatus_only_orthoDB': 2800, 'c_maculatus_all_proteinrefs': 2094, 'c_maculatus_RNA_simple': 138, 'c_maculatus_RNA_combined': 196, 'A_obtectus': 7114}
        }
        Kaufmann_labels = {'C_maculatus__' : "RNA same population", 'c_maculatus_only_orthoDB' : "no RNA", 'c_maculatus_RNA_combined' : "DE RNA different population", 'c_maculatus_RNA_simple' : "RNA different population"}
        # plot_Kaufmann_annotation_comparison(Kaufman_cmac_annotation_comparison, Kaufmann_labels, filename = "Kaufmann_annotation_comparison.png")


        ## plot just the gene numbers from the two (three) annotation methods
        orthoDB_annot = "/Users/miltr339/work/orthoDB_annotations/"
        orthoDB_TE_filtered = "/Users/miltr339/work/orthoDB_proteinseqs_TE_filtered/" # "/proj/naiss2023-6-65/Milena/gene_family_analysis/orthofinder_only_orthoDB_annotations/protein_sequences_TE_filtered/"
        native_annot = "/Users/miltr339/work/native_annotations/all_native_annot/"
        try:
            tree = "/Users/milena/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw"
            plot_gene_counts(native_annot_dir=native_annot, species_tree=tree, orthoDB_filtered_annot_dir=orthoDB_TE_filtered, filename="only_genome_size_14_species.png")
        except:
            tree = "/Users/miltr339/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw"
            plot_gene_counts(native_annot_dir=native_annot, species_tree=tree, orthoDB_filtered_annot_dir=orthoDB_TE_filtered, filename="only_genome_size_14_species.png")
        # plot_gene_counts(native_annot_dir=native_annot, orthoDB_annot_dir=orthoDB_annot, species_names=species_names, orthoDB_filtered_annot_dir=orthoDB_TE_filtered, filename="only_genome_size_14_species_with_TE_filtering.png")

    
    ### plot protein length histograms
    ## filepaths
    if True:
        native_dir = "/Users/miltr339/work/native_proteinseqs"
        native_files = {
            "A_obtectus" : f"{native_dir}/A_obtectus.faa",
            "A_verrucosus" : f"{native_dir}/A_verrucosus.faa",
            "B_siliquastri" : f"{native_dir}/B_siliquastri.faa",
            "C_chinensis" : f"{native_dir}/C_chinensis.faa",
            "C_maculatus" : f"{native_dir}/C_maculatus.faa",
            "C_septempunctata" : f"{native_dir}/C_septempunctata.faa",
            "D_melanogaster" : f"{native_dir}/D_melanogaster.faa",
            "D_ponderosae" : f"{native_dir}/D_ponderosae.faa",
            "I_luminosus" : f"{native_dir}/I_luminosus.faa",
            "P_pyralis" : f"{native_dir}/P_pyralis.faa",
            "R_ferrugineus" : f"{native_dir}/R_ferrugineus.faa",
            "T_castaneum" : f"{native_dir}/T_castaneum.faa",
            "T_molitor" : f"{native_dir}/T_molitor.faa",
            "Z_morio" : f"{native_dir}/Z_morio.faa",
        }

        orthoDB_dir = "/Users/miltr339/work/orthoDB_proteinseqs_TE_filtered"
        orthoDB_files = {
            "A_obtectus" : f"{orthoDB_dir}/A_obtectus_filtered_proteinfasta_TE_filtered.fa",
            "A_verrucosus" : f"{orthoDB_dir}/A_verrucosus_filtered_proteinfasta_TE_filtered.fa",
            "B_siliquastri" : f"{orthoDB_dir}/B_siliquastri_filtered_proteinfasta_TE_filtered.fa",
            "C_analis" : f"{orthoDB_dir}/C_analis_filtered_proteinfasta_TE_filtered.fa",
            "C_chinensis" : f"{orthoDB_dir}/C_chinensis_filtered_proteinfasta_TE_filtered.fa",
            "C_maculatus" : f"{orthoDB_dir}/C_maculatus_filtered_proteinfasta_TE_filtered.fa",
            "C_septempunctata" : f"{orthoDB_dir}/C_septempunctata_filtered_proteinfasta_TE_filtered.fa",
            "D_melanogaster" : f"{orthoDB_dir}/D_melanogaster_filtered_proteinfasta_TE_filtered.fa",
            "D_ponderosae" : f"{orthoDB_dir}/D_ponderosae_filtered_proteinfasta_TE_filtered.fa",
            "I_luminosus" : f"{orthoDB_dir}/I_luminosus_filtered_proteinfasta_TE_filtered.fa",
            "P_pyralis" : f"{orthoDB_dir}/P_pyralis_filtered_proteinfasta_TE_filtered.fa",
            "R_ferrugineus" : f"{orthoDB_dir}/R_ferrugineus_filtered_proteinfasta_TE_filtered.fa",
            "T_castaneum" : f"{orthoDB_dir}/T_castaneum_filtered_proteinfasta_TE_filtered.fa",
            "T_molitor" : f"{orthoDB_dir}/T_molitor_filtered_proteinfasta_TE_filtered.fa",
            "Z_morio" : f"{orthoDB_dir}/Z_morio_filtered_proteinfasta_TE_filtered.fa",
        }

    # get individual plots for all species
    # for species in native_files.keys():
    #     plot_histogram_protein_lengths(native_files[species], orthoDB_files[species], species_name=species, filename = f"protein_lengths_histogram_{species}.png")
    
    # plot_all_species_protein_length_distribution(native_files, orthoDB_files)