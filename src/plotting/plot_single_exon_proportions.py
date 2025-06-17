# make barplot that shows the proportion of single-exon genes in species-specific genes compared to the proportion in the genome overall

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import re
import numpy as np
from pathlib import Path
import parse_gff as gff
from collections import Counter
from Bio import SeqIO



def split_at_second_occurrence(s, char = "_"): # split the gene string at the second occurence of "_" to get only the species name
    if s.count(char)<2:
        return s
    else:
        first_occurrence = s.find(char)
        second_occurrence = s.find(char, first_occurrence + 1)
        species = s[:second_occurrence]
        return species

def make_species_order_from_tree(newick_tree_path):
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


def get_single_exon_genes(dir_path, species_list, write_to_file = False, outfile_name = "", include_total_gene_num = False):
    # get a dictionary with {species_name : [list, of, single, exon, transcripts]} from files generated in 04c_make_single_exon_gene_lists.sh

    directory = Path(dir_path)
    # List all files in the directory
    all_files = [f for f in directory.iterdir() if f.is_file()]

    species_dict = {}
    num_single_exon_genes = {}
    num_transcripts = {}

    all_break = False
    for infile_path in all_files:

        print_statements = False
        infile_path = str(infile_path)
        species = infile_path.split("/")[-1]
        species = species.replace("native_", "")
        species = species.replace("orthoDB_", "")
        species = split_at_second_occurrence(species, "_")
        if len(species.split("_")[0])>1:
            species = species.split("_")[0][0].upper()+"_"+species.split("_")[-1]
        # replace with exact names from the tree if the species exists in the tree 
        species_mod = next((species_name for species_name in species_list if species in species_name), None) # replace the above parsed species name with the corresponding name from species_list
        if species_mod!=None:
            species = species_mod
        
        if len(species.split("_")[0])>1 and len(species.split("_"))>1:
            # if the species name is tribolium_castaneum instead of T_castaneum for example
            parts = species.split("_")
            species = f"{parts[0][0].upper()}_{parts[1]}"

        infile_dict = {}
        list_key = ""
        num_key = ""
        transcripts_num_key = ""
        with open(infile_path, "r") as infile:
            dict_elements = infile.readlines()
            for species_line in dict_elements:
                try:
                    key = species_line.split(":")[0]
                    value = ":".join(species_line.split(":")[1:]) # some of the transcript names have ":" in them to spite me personally
                except:
                    print_statements = True
                    print(species_line[0:150])
                    continue
                value_list = [gene.split("|")[-1] for gene in value.strip().split(",")] # TODO # native annotations have some stuff going on before "|" and to actually make it work you need to get rid of that
                infile_dict[key]=value_list
                if "list" in key:
                    list_key = key
                if "single exon transcripts" in key:
                    num_key = key
                #if "number of transcripts" in key:
                if "gene features" in key:    
                    transcripts_num_key = key
                
 
            if all_break:
                break
            if print_statements:
                print(species)
                print(list_key)
                print(infile_dict.keys())
                print_statements = False
        # species_dict[species] = infile_dict["list of single-exon transcript IDs"]
        try:
            species_dict[species] = infile_dict[list_key]
        except:
            print(f"{species} single-exon list did not work with {list_key}")
        
        try:
            num_single_exon_genes[species] = int(infile_dict[num_key][0])
        except:
            print(f"{species} single-exon number did not work with {num_key}")

        if include_total_gene_num:
            try:
                num_transcripts[species] = int(infile_dict[transcripts_num_key][0])
            except:
                keys = list(infile_dict.keys())
                print(f"{species} transcript number did not work with {transcripts_num_key}, available keys are: {keys}")


    if write_to_file:
        gff.write_dict_to_file(species_dict, f"{dir_path}/../{outfile_name}")

    # num_single_exon_genes = {species : len(IDs) for species, IDs in species_dict.items()}
    if not include_total_gene_num:
        return(species_dict, num_single_exon_genes)
    elif include_total_gene_num:
        return(species_dict, num_single_exon_genes, num_transcripts)


def get_species_specific_genes(orthogroups_filepath, species_list, write_to_file = False, outfile_name = ""):
    # get a dictionary with {species_name : [list, of, species-specific, transcripts]}    

    list_spec_genes = {} # {species : [list, of, species-specific, transcripts]}
    num_spec_genes = {} # dictionary of species-specific genes {species : number_of_species-specific orthogroups}

    with open(orthogroups_filepath+"/Orthogroups.txt", "r") as orthogroups_file:
        orthogroups_lines = orthogroups_file.readlines()
        
        complete_gene_count_dict = {}
        complete_gene_count = 0

        # set up dictionary
        for species in species_list:
            list_spec_genes[species] = []
            complete_gene_count_dict[species] = 0
            
        for orthogroup in orthogroups_lines:
            orthogroup_members = orthogroup.split(": ")[1].strip().split(" ")
            complete_gene_count = complete_gene_count+len(orthogroup_members)

            # check if it's a single-species orthogroup
            # og_species_list = list(set([split_at_second_occurrence(gene, "_") for gene in orthogroup_members]))
            # og_species_list = [gene.rsplit("_",1)[0] for gene in orthogroup_members]
            og_species_list = [split_at_second_occurrence(gene, "_") for gene in orthogroup_members]
            og_species_list = [next((species_name for species_name in species_list if species_name in gene_species), gene_species) for gene_species in og_species_list]
            
            # count number of species in all genes
            all_species_counter = Counter(og_species_list)
            for species in species_list:
                complete_gene_count_dict[species] = complete_gene_count_dict[species]+all_species_counter[species]
            
            og_species_list = list(set(og_species_list)) 
            if len(og_species_list)<2 and og_species_list[0] in species_list:

                species = og_species_list[0]
                species = next((species_name for species_name in species_list if species in species_name), species) # replace the above parsed species name with the corresponding name from species_list
                # species_genes = [gene.replace("__", "_")[:-2].split(species+"_")[-1] for gene in orthogroup_members if species in gene] # get gene IDs for each species, replace the instances of "__" with "_" remove the "_1" at the end that orthofinder for some reason adds
                species_genes = []
                for gene in orthogroup_members:
                    # if species in gene:
                    species_genes.append(gene[:-2].replace(f"{species}_", "")) # remove the "_1" that orthofinder adds after each ID, and the species name
                list_spec_genes[species].extend(species_genes)

        # fill the number of species specific genes into a dictionary        
        for species in species_list:
            num_spec_genes[species] = len(list_spec_genes[species])

        if write_to_file:
            gff.write_dict_to_file(list_spec_genes, f"{orthogroups_filepath}/{outfile_name}")

    return(list_spec_genes, num_spec_genes, complete_gene_count_dict)


def check_proteinfasta_single_exon(proteinfasta_filepath, annotation_filepath = "", species = ""):
    fasta_headers_list = []
    transcript_IDs_single_exon = []
    for record in SeqIO.parse(proteinfasta_filepath, "fasta"):
        transcript_id = record.id[:-2]
        transcript_id = transcript_id.split(f"{species}_")[-1] # remove the species name and "_1" suffix from the orthfinder input files
        fasta_headers_list.append(transcript_id)
    # print(fasta_headers_list[:10])
    # print(f"{len(fasta_headers_list)} transcripts in {proteinfasta_filepath}")
    gff_parsed = gff.parse_gff3_general(annotation_filepath, verbose=False)
    # gff_parsed[fasta_headers_list[0]].show
    # print(gff_parsed[fasta_headers_list[0]].features[gff.FeatureCategory.Exon])
    for fasta_header in fasta_headers_list:
        if gff_parsed[fasta_header].category == gff.FeatureCategory.Gene:
            if len(gff_parsed[fasta_header].child_ids_list) == 1:
                transcript_IDs_single_exon.append(fasta_header)

    return(transcript_IDs_single_exon)

    


def get_transcript_number(proteinfasta_filepath):
    fasta_headers_list = []
    for record in SeqIO.parse(proteinfasta_filepath, "fasta"):
        fasta_headers_list.append(record.id)
    return(len(fasta_headers_list))


def read_gene_dicts(infile_path):
    species_dict = {}
    number_of_genes_dict = {}
    with open(infile_path, "r") as infile:
        dict_elements = infile.readlines()
        for species_line in dict_elements:
            species, genes = species_line.split(":")
            gene_id_list = [gene.split(".")[0] for gene in genes.split(",")] # remove the transcript ID at the end if it exists
            species_dict[species]=gene_id_list
            number_of_genes_dict[species]=len(gene_id_list)
    return(species_dict, number_of_genes_dict)


# get a dictionary with species specific genes that are also single-exon
def get_gene_overlap(single_exon_dict, species_specific_dict, verbose = True):

    overlapping_species = [single_exon_species for single_exon_species in list(single_exon_dict.keys()) if single_exon_species in list(species_specific_dict.keys())]
    single_exon_exclusive = [single_exon_species for single_exon_species in list(single_exon_dict.keys()) if single_exon_species not in list(species_specific_dict.keys())]
    species_specific_exclusive = [species_specifc_species for species_specifc_species in list(species_specific_dict.keys()) if species_specifc_species not in list(species_specific_dict.keys())]
    if verbose:
        print(f"overlapping species: {len(overlapping_species)}")
        print(f"single-exon genes species: {len(single_exon_dict.keys())}\t,\t species-specific genes species: {len(species_specific_dict.keys())}")
        if len(single_exon_exclusive)>0:
            print(f"\tspecies only present in the single-exon input: {single_exon_exclusive}")
        if len(species_specific_exclusive):
            print(f"\tspecies only present in the species specifc input: {species_specific_exclusive}")

    # initialize output dictionary that has values that are gene IDs that are species specific and only have one exon
    overlap_dict = {}
    number_of_genes_dict = {} 

    for species in overlapping_species:
        species = next((species_name for species_name in species_names if species in species_name), None) # replace the above parsed species name with the corresponding name from species_names
        overlap = [gene for gene in species_specific_dict[species] if gene in single_exon_dict[species]]
        overlap_dict[species] = overlap
        number_of_genes_dict[species] = len(overlap)
        # break
    
    return(overlap_dict, number_of_genes_dict)




### formatting of gene number y ticks
def format_func(value, tick_number):
    if value >= 1e6:
        return f'{value/1e6:.0f}M'
    elif value >= 1e3:
        return f'{value/1e3:.0f}k'
    else:
        return str(value)



def plot_single_exon_no_species_specific_three_annot(native_numbers, orthodb_numbers = {}, orthoDB_unmasked_numbers = {}, species_names = [], filename = "", L50_values = {}):

    print(f" plotting for these {len(species_names)} species: \n{species_names}")

    # X coordinates for the groups
    x = np.arange(len(species_names))

    # figure proportions according to the data included (longer or shorter)

    # fontsize scales with the dpi somehow which i have to do extra because i change the aspect ratio manually below
    fs = 35 # 37 originally
    
    if len(orthodb_numbers) == 0 and len(orthoDB_unmasked_numbers) == 0:
        width = 0.35 # (this is a fraction of the standardized 1 unit of space between axis ticks)
        aspect_ratio = 17 / 12 # nice for a presentation
    elif len(orthodb_numbers) > 0 and len(orthoDB_unmasked_numbers) > 0:
        width = 0.22
        fs = fs*1.15 
        aspect_ratio = 24 / 12
    else:
        width = 0.175
        fs = fs*0.9 # a smaller fontsize is nicer for the longer figure since it needs to be shown larger overall anyways
        aspect_ratio = 27 / 12

    height_pixels = 2000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig = plt.figure(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    ax = fig.add_subplot(111)


    colors = {
        "native" : "#b82946", 
        "orthodb" : "#F2933A", 
        "orthoDB_unmasked" : "#4d7298",
    }
    hatch_color = '#ffffff' # '#E2D4CA' #kind of eggshell white
    plt.rcParams['hatch.color'] = hatch_color


    #### plot native annotation ####

    num_single_exon_genes = native_numbers["num_single_exon_genes"]
    num_total_genes = native_numbers["num_total_genes"]

    single_exon_genes = [num_single_exon_genes[species] for species in species_names]
    multi_exon_genes = [num_total_genes[species]-num_single_exon_genes[species] for species in species_names]

    color = [colors["native"], colors["native"]]
    hatching = ["//", "" ]


    if len(orthodb_numbers)>0 and len(orthoDB_unmasked_numbers) == 0:
        category = " (native)"
        x_subtr = width/2
    if len(orthodb_numbers)>0 and len(orthoDB_unmasked_numbers) > 0:
        category = " (native)"
        x_subtr = width # 2*width
    else:
        category = ""
        x_subtr = 0
    
    print(f"bar width = {width}")

    # total number of genes (with single-exons hatched)
    native_rects1_base = ax.bar(x - x_subtr, single_exon_genes, width, label='proportion of which are single-exon', color= color[0], hatch=hatching[0])
    native_rects1_top = ax.bar(x - x_subtr, multi_exon_genes, width, bottom=single_exon_genes, label='all genes'+category, color= color[1], hatch=hatching[1])


    #### plot orthoDB uniform masking annotation ####

    if len(orthodb_numbers) > 0:
        num_single_exon_genes = orthodb_numbers["num_single_exon_genes"]
        num_total_genes = orthodb_numbers["num_total_genes"]

        single_exon_genes = [num_single_exon_genes[species] for species in species_names]
        multi_exon_genes = [num_total_genes[species]-num_single_exon_genes[species] for species in species_names]

        color = [colors["orthodb"], colors["orthodb"]]
        hatching = ["//", ""]

        if len(orthoDB_unmasked_numbers) == 0:
            # total number of genes (with single-exons hatched)            
            orthodb_rects1_base = ax.bar(x + width/2, single_exon_genes, width, label='proportion of which are single-exon', color= color[0], hatch=hatching[0])            
            orthodb_rects1_top = ax.bar(x + width/2, multi_exon_genes, width, bottom=single_exon_genes, label='all genes (uniform)', color= color[1], hatch=hatching[1])  
        elif len(orthoDB_unmasked_numbers) > 0:
            # total number of genes (with single-exons hatched)            
            orthodb_rects1_base = ax.bar(x, single_exon_genes, width, label='proportion of which are single-exon', color= color[0], hatch=hatching[0])            
            orthodb_rects1_top = ax.bar(x, multi_exon_genes, width, bottom=single_exon_genes, label='all genes (uniform)', color= color[1], hatch=hatching[1])  

        plt.rcParams.update({'hatch.color': hatch_color})
        ymax = max(num_total_genes.values())*1.25
    
    
    #### plot orthoDB filtered annotation ####

    if len(orthoDB_unmasked_numbers) > 0:
        num_single_exon_genes = orthoDB_unmasked_numbers["num_single_exon_genes"]
        num_total_genes = orthoDB_unmasked_numbers["num_total_genes"]

        single_exon_genes = [num_single_exon_genes[species] for species in species_names]
        multi_exon_genes = [num_total_genes[species]-num_single_exon_genes[species] for species in species_names]

        color = [colors["orthoDB_unmasked"], colors["orthoDB_unmasked"]]
        hatching = ["//", ""]
           
        orthodb_rects1_base = ax.bar(x + x_subtr, single_exon_genes, width, label='proportion of which are single-exon', color= color[0], hatch=hatching[0])            
        orthodb_rects1_top = ax.bar(x + x_subtr, multi_exon_genes, width, bottom=single_exon_genes, label='all genes (uniform no remasking)', color= color[1], hatch=hatching[1])  
        plt.rcParams.update({'hatch.color': hatch_color})
        ymax = max(num_total_genes.values())*1.25


    #### set up labels and stuff ####
    
    ax.set_ylabel('Number of genes', fontsize=fs+4)
    ax.set_title('Proportion of single-exon genes', fontsize=fs+4)
    ax.set_xticks(x)
    if len(L50_values) > 0:
        ax.set_xlabel('Species (L50 value of the assembly)', fontsize=fs+4)
        xtick_labels = [species.replace("_", ". ")+f" ({L50_values[species]})" for species in species_names]
    else:
        ax.set_xlabel('', fontsize=fs+4)
        xtick_labels = [species.replace("_", " ") for species in species_names]
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=fs)
    ax.set_yticklabels([f'{int(tick)/1e3:.0f}k' for tick in ax.get_yticks()], fontsize=fs)

    # make custom legend patch for the dashed bars
    plt.rcParams.update({'hatch.color': "#3f3832ff"})
    dashed_handle = mpatches.Patch(hatch = "//", alpha = 0.0)
    dashed_label = "proportion of genes that are single-exon"

    # Legend with custom order
    handles, labels = ax.get_legend_handles_labels()

    new_order = [1,3,5]
    handles = [handles[idx] for idx in new_order]
    labels = [labels[idx] for idx in new_order]
    handles.append(dashed_handle)
    labels.append(dashed_label)
    new_order = [0,2,1,3]
    handles = [handles[idx] for idx in new_order]
    labels = [labels[idx] for idx in new_order]

    ax.legend(handles, labels, fontsize=fs, ncol=2, loc='upper center')

    # add space at the top of the plot for the legend
    ax.set_ylim(0, int(ymax))
    ax.set_xlim(-0.5, len(xtick_labels)-0.5)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)



if __name__ == '__main__':

    ## calculate all the numbers from input data

    species_names = make_species_order_from_tree("/Users/miltr339/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw")

    # from quast results for all assemblies
    L50_values = {
        "A_obtectus" : 5,
        "A_verrucosus" : 8944,
        "B_siliquastri" : 5,
        "C_analis" : 1051,
        "C_chinensis" : 231,
        "C_maculatus" : 6, # the non-superscaffolded value is 38,
        "C_septempunctata" : 4,
        "D_melanogaster" : 3,
        "D_ponderosae" : 87,
        "I_luminosus" : 1966,
        "P_pyralis" : 5,
        "R_ferrugineus" : 268,
        "T_castaneum" : 5,
        "T_molitor" : 6,
        "Z_morio" : 4
    }

    ### original plotting with species specific and single-exon genes and their intersection
    if False:

        single_exon_dict_orthoDB, num_single_exon_dict_orthoDB = get_single_exon_genes("/Users/milena/work/single_exon_genes/orthoDB", species_names, write_to_file=False, outfile_name="orthoDB_single_exon_transcripts_list_14_species.txt")
        single_exon_dict_orthoDB_unmasked, num_single_exon_dict_orthoDB_unmasked = get_single_exon_genes("/Users/milena/work/single_exon_genes/orthoDB_old", species_names, write_to_file=False, outfile_name="orthoDB_single_exon_transcripts_list_14_species.txt")
        single_exon_dict_native, num_single_exon_dict_native = get_single_exon_genes("/Users/milena/work/single_exon_genes/native", species_names, write_to_file=False, outfile_name="native_single_exon_transcripts_list_14_species.txt")
        
        ## TODO: get for filtered annotations
        ## filter the orthoDB dictionaries to include only the overlap-filtered transcripts
        orthoDB_filtered_proteinseqs_dir = "/Users/milena/work/orthoDB_proteinseqs/overlap_filtered_proteinseqs/"
        filtered_proteinfasta = {
            "A_obtectus" : f"{orthoDB_filtered_proteinseqs_dir}A_obtectus_filtered_proteinfasta_overlap_filtered.fa",
            "A_verrucosus" : f"{orthoDB_filtered_proteinseqs_dir}A_verrucosus_filtered_proteinfasta_overlap_filtered.fa",
            "B_siliquastri" : f"{orthoDB_filtered_proteinseqs_dir}B_siliquastri_filtered_proteinfasta_overlap_filtered.fa",
            "C_chinensis" : f"{orthoDB_filtered_proteinseqs_dir}C_chinensis_filtered_proteinfasta_overlap_filtered.fa",
            "C_maculatus" : f"{orthoDB_filtered_proteinseqs_dir}C_maculatus_filtered_proteinfasta_overlap_filtered.fa",
            "C_septempunctata" : f"{orthoDB_filtered_proteinseqs_dir}C_septempunctata_filtered_proteinfasta_overlap_filtered.fa",
            "D_melanogaster" : f"{orthoDB_filtered_proteinseqs_dir}D_melanogaster_filtered_proteinfasta_overlap_filtered.fa",
            "D_ponderosae" : f"{orthoDB_filtered_proteinseqs_dir}D_ponderosae_filtered_proteinfasta_overlap_filtered.fa",
            "I_luminosus" : f"{orthoDB_filtered_proteinseqs_dir}I_luminosus_filtered_proteinfasta_overlap_filtered.fa",
            "P_pyralis" : f"{orthoDB_filtered_proteinseqs_dir}P_pyralis_filtered_proteinfasta_overlap_filtered.fa",
            "R_ferrugineus" : f"{orthoDB_filtered_proteinseqs_dir}R_ferrugineus_filtered_proteinfasta_overlap_filtered.fa",
            "T_castaneum" : f"{orthoDB_filtered_proteinseqs_dir}T_castaneum_filtered_proteinfasta_overlap_filtered.fa",
            "T_molitor" : f"{orthoDB_filtered_proteinseqs_dir}T_molitor_filtered_proteinfasta_overlap_filtered.fa",
            "Z_morio" : f"{orthoDB_filtered_proteinseqs_dir}Z_morio_filtered_proteinfasta_overlap_filtered.fa"
        }
        orthoDB_annotations = {
            "A_obtectus" : "/Users/milena/work/orthoDB_annotations/A_obtectus_isoform_filtered.gff",
            "A_verrucosus" : "/Users/milena/work/orthoDB_annotations/A_verrucousus_isoform_filtered.gff",
            "B_siliquastri" : "/Users/milena/work/orthoDB_annotations/B_siliquastri_isoform_filtered.gff",
            "C_analis" : "/Users/milena/work/orthoDB_annotations/C_analis_isoform_filtered.gff",
            "C_chinensis" : "/Users/milena/work/orthoDB_annotations/C_chinensis_isoform_filtered.gff",
            "C_maculatus" : "/Users/milena/work/orthoDB_annotations/C_maculatus_isoform_filtered.gff",
            "C_septempunctata" : "/Users/milena/work/orthoDB_annotations/C_septempunctata_isoform_filtered.gff",
            "D_melanogaster" : "/Users/milena/work/orthoDB_annotations/D_melanogaster_isoform_filtered.gff",
            "D_ponderosae" : "/Users/milena/work/orthoDB_annotations/D_ponderosae_isoform_filtered.gff",
            "I_luminosus" : "/Users/milena/work/orthoDB_annotations/I_luminosus_isoform_filtered.gff",
            "P_pyralis" : "/Users/milena/work/orthoDB_annotations/P_pyralis_isoform_filtered.gff",
            "R_ferrugineus" : "/Users/milena/work/orthoDB_annotations/R_ferrugineus_isoform_filtered.gff",
            "T_castaneum" : "/Users/milena/work/orthoDB_annotations/T_castaneum_isoform_filtered.gff",
            "T_molitor" : "/Users/milena/work/orthoDB_annotations/T_molitor_isoform_filtered.gff",
            "Z_morio" : "/Users/milena/work/orthoDB_annotations/Z_morio_isoform_filtered.gff"
    }
        
        single_exon_dict_filtered_orthoDB = {}
        num_single_exon_dict_filtered_orthoDB = {}

        for species in filtered_proteinfasta.keys():
            print(f"{species}:")
            single_exon_dict_filtered_orthoDB[species] = check_proteinfasta_single_exon(filtered_proteinfasta[species], annotation_filepath=orthoDB_annotations[species], species=species)
            num_single_exon_dict_filtered_orthoDB[species] = len(single_exon_dict_filtered_orthoDB[species])
            print(f"\tnumber of single-exon transcripts unfiltered: {len(single_exon_dict_orthoDB[species])} \n\tnumber of single-exon transcripts filtered: {len(single_exon_dict_filtered_orthoDB[species])}")
        



        species_specific_dict_orthoDB, num_species_specific_genes_orthoDB, species_gene_numbers_orthoDB = get_species_specific_genes("/Users/milena/work/orthofinder/orthoDB_orthogroups", species_names, write_to_file=False, outfile_name="species_specific_transcripts_lists_orthoDB.txt")
        species_specific_dict_native, num_spceies_specific_genes_native, species_gene_numbers_native = get_species_specific_genes("/Users/milena/work/orthofinder/native_orthogroups", species_names, write_to_file=False, outfile_name="species_specific_transcripts_lists_native.txt")
        if False:
            species_specific_dict_filtered_orthoDB = {}
            num_species_specific_dict_filtered_orthoDB = {}
            for species in filtered_proteinfasta.keys():
                species_specific_dict_filtered_orthoDB[species] = filter_transcript_list_by_proteinfasta(species_specific_dict_orthoDB[species], filtered_proteinfasta[species], species=species)
                print(f"{species}: \n\tnumber of species specific transcripts unfiltered: {len(species_specific_dict_orthoDB[species])} \n\tnumber of species specific transcripts filtered: {len(species_specific_dict_filtered_orthoDB[species])}")
                num_species_specific_dict_filtered_orthoDB[species] = len(species_specific_dict_filtered_orthoDB[species])
        if True:
            species_gene_numbers_filtered_orthoDB = {}
            for species in filtered_proteinfasta.keys():
                species_gene_numbers_filtered_orthoDB[species] = get_transcript_number(filtered_proteinfasta[species])
                print(f"{species} has {species_gene_numbers_filtered_orthoDB[species]} transcripts in the orthoDB annotation after filtering")


        # overlap_dict_orthoDB, number_of_overlap_genes_orthoDB = get_gene_overlap(single_exon_dict=single_exon_dict_orthoDB, species_specific_dict=species_specific_dict_orthoDB, verbose = False)
        # overlap_dict_native, number_of_overlap_genes_native = get_gene_overlap(single_exon_dict=single_exon_dict_native, species_specific_dict=species_specific_dict_native, verbose = False)



        ## put all the numbers together to sort the plotting out easier
        if False:
            orthodb_numbers = {
                "num_single_exon_genes" : num_single_exon_dict_orthoDB,
                "num_species_spec_genes" : num_species_specific_genes_orthoDB,
                "num_total_genes" : species_gene_numbers_orthoDB,
                "num_species_spec_genes_with_one_exon" : number_of_overlap_genes_orthoDB
            }
            native_numbers = {
                "num_single_exon_genes" : num_single_exon_dict_native,
                "num_species_spec_genes" : num_spceies_specific_genes_native,
                "num_total_genes" : species_gene_numbers_native,
                "num_species_spec_genes_with_one_exon" : number_of_overlap_genes_native
            }
            # plot_single_exon_species_spec_genes_both_annot(native_numbers, orthodb_numbers, species_names, filename = "single-exon_genes_barplot_14_species.png", L50_values = L50_values)
            plot_single_exon_species_spec_genes_both_annot(native_numbers, orthodb_numbers, species_names, filename = "single-exon_genes_barplot_14_species_no_species_specific.png", L50_values = L50_values, include_species_specific=False)

        else:
            # TODO double check single exon genes with gff parser

            orthodb_numbers = {
                "num_single_exon_genes" : num_single_exon_dict_orthoDB,
                "num_total_genes" : species_gene_numbers_orthoDB,
                "num_species_spec_genes" : L50_values, # just a "nonsense" placeholder dictionary so that the iteration works, will not be plotted if "include_species_specific = False"
                "num_species_spec_genes_with_one_exon" : L50_values
            }
            orthodb_filtered_numbers = {
                "num_single_exon_genes" : num_single_exon_dict_filtered_orthoDB,
                "num_total_genes" : species_gene_numbers_filtered_orthoDB,
                "num_species_spec_genes" : L50_values,
                "num_species_spec_genes_with_one_exon" : L50_values
            }
            native_numbers = {
                "num_single_exon_genes" : num_single_exon_dict_native,
                "num_total_genes" : species_gene_numbers_native,
                "num_species_spec_genes" : L50_values,
                "num_species_spec_genes_with_one_exon" : L50_values
            }
        
        plot_single_exon_species_spec_genes_both_annot(native_numbers, orthodb_numbers=orthodb_numbers, orthoDB_filtered_numbers=orthodb_filtered_numbers, species_names=species_names, filename = "single-exon_genes_barplot_14_species_three_annot.png", L50_values = L50_values, include_species_specific=False)

    if True:

        single_exon_dict_orthoDB, num_single_exon_dict_orthoDB, num_transcripts_dict_orthoDB = get_single_exon_genes("/Users/miltr339/work/single_exon_genes/orthoDB", species_names, write_to_file=False, outfile_name="orthoDB_single_exon_transcripts_list_14_species.txt", include_total_gene_num = True)
        single_exon_dict_orthoDB_unmasked, num_single_exon_dict_orthoDB_unmasked, num_transcripts_dict_orthoDB_unmasked = get_single_exon_genes("/Users/miltr339/work/single_exon_genes/orthoDB_old", species_names, write_to_file=False, outfile_name="orthoDB_single_exon_transcripts_list_14_species.txt", include_total_gene_num = True)
        single_exon_dict_native, num_single_exon_dict_native, num_transcripts_dict_native = get_single_exon_genes("/Users/miltr339/work/single_exon_genes/native", species_names, write_to_file=False, outfile_name="native_single_exon_transcripts_list_14_species.txt", include_total_gene_num = True)

        orthodb_numbers = {
            "num_single_exon_genes" : num_single_exon_dict_orthoDB,
            "num_total_genes" : num_transcripts_dict_orthoDB
        }
        print(f"\t --> orthodb_numbers = {orthodb_numbers}")
        orthodb_unmasked_numbers = {
            "num_single_exon_genes" : num_single_exon_dict_orthoDB_unmasked,
            "num_total_genes" : num_transcripts_dict_orthoDB_unmasked
        }
        print(f"\t --> orthodb_unmasked_numbers = {orthodb_unmasked_numbers}")
        native_numbers = {
            "num_single_exon_genes" : num_single_exon_dict_native,
            "num_total_genes" : num_transcripts_dict_native
        }
        print(f"\t --> native_numbers = {native_numbers}")
        print()
        

        data = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/"
        plot_single_exon_no_species_specific_three_annot(native_numbers, orthodb_numbers = orthodb_numbers, orthoDB_unmasked_numbers = orthodb_unmasked_numbers, species_names = species_names, filename = f"{data}single_exon_Genes_14_species_3_annotations.png", L50_values = L50_values)