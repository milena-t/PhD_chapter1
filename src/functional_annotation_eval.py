"""
Plot some stuff to evaluate the functional annotation of the gene groups
"""

import plot_basics as my_plotting
import parse_orthogroups as OGs
import plotting.plot_significant_orthogroups_from_CAFE as plot_OG
import make_orthogroup_flybaseID_table as parse_DAVID
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt
import numpy as np
import base64
import re
from bs4 import BeautifulSoup

from collections import Counter


def filepaths():
    orthogroups_orthoDB_filepath = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv"
    out_dir = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/"
    tree = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_only_species_names.nw"
    DAVID_path = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/functional_table_per_OG.tsv"
    return out_dir,orthogroups_orthoDB_filepath,tree, DAVID_path

def filepaths_work():
    orthogroups_orthoDB_filepath = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv"
    out_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/"
    tree = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_only_species_names.nw"
    DAVID_path = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/functional_table_per_OG.tsv"
    return out_dir,orthogroups_orthoDB_filepath,tree, DAVID_path

def orthogroups_lists():
    out_dict_old_single_CAFE_run = {
        "Aobt_expansion" : ["N0.HOG0000035","N0.HOG0000014"],
        "Gene Group 1" : ["N0.HOG0000027","N0.HOG0000059","N0.HOG0000095","N0.HOG0000140","N0.HOG0000204","N0.HOG0000492","N0.HOG0001030","N0.HOG0001077"],
        "Gene Group 3" : ["N0.HOG0000086","N0.HOG0002393","N0.HOG0000085"],
        "Gene Group 18" : ["N0.HOG0000669"],
        "Gene Group 5,8" : ["N0.HOG0000541"],
        "Gene Group 8" : ["N0.HOG0000775","N0.HOG0000892"],
        "Gene Group 9" : ["N0.HOG0000401","N0.HOG0009002"],
        "Gene Group 20" : ["N0.HOG0000030","N0.HOG0000450","N0.HOG0000526","N0.HOG0000756","N0.HOG0000761","N0.HOG0010885"],
        "Gene Group 30": ["N0.HOG0000037","N0.HOG0000177","N0.HOG0001445","N0.HOG0000038","N0.HOG0000194","N0.HOG0000345","N0.HOG0000467"],
        "Gene Group 7" : ["N0.HOG0000056","N0.HOG0000454","N0.HOG0000436","N0.HOG0000480","N0.HOG0009039"],
        "Gene Group 15": ["N0.HOG0000307","N0.HOG0001194","N0.HOG0003035"],
        "Gene Group 24": ["N0.HOG0000108","N0.HOG0000039","N0.HOG0000044","N0.HOG0001108"],
        "Gene Group 26" : ["N0.HOG0000078","N0.HOG0000112","N0.HOG0000113","N0.HOG0000173","N0.HOG0000174","N0.HOG0000215","N0.HOG0000276","N0.HOG0000305","N0.HOG0000385","N0.HOG0000424","N0.HOG0000445","N0.HOG0000462","N0.HOG0000506","N0.HOG0000555","N0.HOG0000560","N0.HOG0000562","N0.HOG0000582","N0.HOG0000590","N0.HOG0000624","N0.HOG0000653","N0.HOG0000662","N0.HOG0000668","N0.HOG0000709","N0.HOG0000715","N0.HOG0000730","N0.HOG0000798","N0.HOG0000836","N0.HOG0000854","N0.HOG0000859","N0.HOG0000885","N0.HOG0000888","N0.HOG0000890","N0.HOG0000902","N0.HOG0000915","N0.HOG0000916","N0.HOG0000960","N0.HOG0001028","N0.HOG0001040","N0.HOG0001065","N0.HOG0001094","N0.HOG0001097","N0.HOG0001100","N0.HOG0001311","N0.HOG0001487","N0.HOG0001719","N0.HOG0001774","N0.HOG0001781","N0.HOG0001850","N0.HOG0001925","N0.HOG0001928","N0.HOG0002240","N0.HOG0002511","N0.HOG0002693","N0.HOG0003778","N0.HOG0004409","N0.HOG0005613","N0.HOG0008866","N0.HOG0010108","N0.HOG0010919","N0.HOG0012168","N0.HOG0013047","N0.HOG0014474","N0.HOG0017274","N0.HOG0000415","N0.HOG0001112","N0.HOG0001856"],
        "Gene Group 13" : ["N0.HOG0000278","N0.HOG0001036"],
        "Gene Group 4": ["N0.HOG0000120","N0.HOG0000141","N0.HOG0000365","N0.HOG0000378","N0.HOG0007183"],
        "Acyl_CoA_synthesis" : ["N0.HOG0000284","N0.HOG0000397","N0.HOG0000613"],
        "Gene Group 11" : ["N0.HOG0000525"],
    }
    out_dict = {
        "Aobt_expansion" : ["N0.HOG0000035","N0.HOG0000014"],
        "Acyl_CoA_synthesis" : ["N0.HOG0000284","N0.HOG0000397","N0.HOG0000613"],
    }
    return out_dict_old_single_CAFE_run

def plot_selected_OGs(orthogroups_path:str, OG_IDs:list[list[str]], colors:list, labels:list, tree_path:str, filename:str, out_dir:str = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/", title = "", transparent_bg=True, svg = False):
    """
    plot the gene family sizes in a subset of orthogroup_IDs
    takes a list of lists of OG IDs, grouped by which should have the same color. 
    The colors list is then these colors, in the same order as the lists in OG_IDs
    """
    fs = 22 # set font size
    # plot each column in the dataframe as a line in the same plot thorugh a for-loop
    fig = plt.figure(figsize=(15,10))
    ax = fig.add_subplot(1, 1, 1)

    # get sorted species names from tree
    species_names_unsorted = my_plotting.plot_tree_manually(tree_path)
    species_coords_sorted = sorted(list(species_names_unsorted.keys()))
    species_names = [species_names_unsorted[species_coord] for species_coord in species_coords_sorted]

    orthoDB_dict_lists = OGs.parse_orthogroups_dict(orthogroups_path)
    orthoDB_dict = OGs.get_GF_sizes(orthoDB_dict_lists)
    
    plot_name = f"{out_dir}{filename}"

    if type(OG_IDs[0]) != list:
        OG_IDs = [OG_IDs]
    if len(OG_IDs) != len(colors) and len(colors) == len(labels):
        raise RuntimeError(f"\t --> list with orthogroup IDs lists and list with colors do not have the same length!")
    elif len(OG_IDs) != len(labels) and len(colors) == len(labels):
        raise RuntimeError(f"\t --> list with orthogroup IDs lists and list with labels do not have the same length!")
    elif len(colors) != len(labels):
        raise RuntimeError(f"\t --> list with colors and list with labels do not have the same length!")

    ymax = 0
    for i,OG_IDs_list in enumerate(OG_IDs):
        
        ax.plot(0, 0, color = colors[i], alpha = 1, linewidth =2, label = labels[i]) 
        OGs_of_interest_dict = {OG_id : orthoDB_dict[OG_id] for OG_id in OG_IDs_list}

        for orthogroup in OG_IDs_list:
            gene_family_members = []
            for species in species_names:
                try:
                    gene_family_members.append(OGs_of_interest_dict[orthogroup][species])
                except:
                    gene_family_members.append(0)
        
            if max(gene_family_members) > ymax:
                ymax = max(gene_family_members)

            ax.plot(species_names, gene_family_members, color = colors[i], alpha = 1, linewidth =2) # originally 0.8


    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x < 0  else f'{int(x)}'))

    ylab="number of gene family members"
    ax.set_ylabel(ylab, fontsize = fs)
    plt.xticks(labels=[species.replace("_", ". ") for species in species_names], ticks=species_names, rotation = 90, fontsize = fs)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, fontsize = fs, loc='upper center')

    plt.title(title, fontsize = fs)

    #ax.legend(fontsize = fs)
    # set grid only for X axis ticks 
    ax.grid(True)
    ax.yaxis.grid(False)

    # ymax = ymax*1.25
    ymax = ymax*1.5

    if ymax>15:
        ax.set_ylim(-4.5,ymax)
    else:
        ax.set_ylim(-1.5,ymax)
    ax.tick_params(axis='y', labelsize=fs)

    plt.tight_layout()

    if svg:
        plot_name = plot_name.replace(".png", ".svg")
        plt.savefig(plot_name, transparent = transparent_bg)
    else:
        plt.savefig(plot_name, dpi = 300, transparent = transparent_bg)
    print("Figure saved in the current working directory directory as: "+plot_name)
    
    return plot_name


def inline_svgs_in_html(html_path, output_path = ""):
    if output_path == "":
        output_path = html_path.replace(".html", "_svg_embed.html")
    with open(html_path, 'r', encoding='utf-8') as f:
        html_content = f.read()

    # Find all <img src="...svg">
    pattern = r'<img\s+[^>]*src="([^"]+\.svg)"[^>]*>'
    
    def replace_img_with_svg(match):
        svg_path = match.group(1)
        img_tag = match.group(0)
        src = match.group(1)
        try:
            try:
                with open(svg_path, 'r', encoding='utf-8') as svg_file:
                    svg_content = svg_file.read()
            except:
                svg_path = svg_path.replace("file:///", "")
                with open(svg_path, 'r', encoding='utf-8') as svg_file:
                    svg_content = svg_file.read()
            img_soup = BeautifulSoup(img_tag, 'html.parser').img
            svg_soup = BeautifulSoup(svg_content, 'html.parser').svg

            # Transfer relevant attributes (like width, height, style, class)
            for attr in ['width', 'height', 'style', 'class']:
                if img_soup.has_attr(attr):
                    svg_soup[attr] = img_soup[attr]

            if 'height' in svg_soup.attrs:
                del svg_soup['height']

            return str(svg_soup)
            # return svg_content
        except FileNotFoundError:
            print(f"⚠️ Warning: SVG file not found: {svg_path}")
            return match.group(0)  # keep original <img> if file not found

    # Replace all matched <img> tags with the actual SVG content
    updated_html = re.sub(pattern, replace_img_with_svg, html_content)

    # Write to output
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(updated_html)

    print(f"--> Inlined SVGs written to: {output_path}")




def get_base64_for_html_embed(image_path):
    with open(image_path, "rb") as image:
        encoded = base64.b64encode(image.read()).decode('utf-8')
    return encoded


def get_OG_david_annot(table_path:str):
    """
    parse the table created from the flybase ID functional annotation into a dictionary
    orthogroups = { orthogroup_ID : [ Gene_Group_A , Gene_Group_B , ... ] }
    gene_groups = { gene_group : "function" }
    """
    orthogroups = {}
    gene_groups = {}
    with open(table_path, "r") as table:
        table_lines = table.readlines()
        for line in table_lines[1:]:
            Orthogroup_ID,CAFE_p_value,Gene_Group,Group_function,repeat_correlation_slope,repeat_correlation_p_value,GS_correlation_slope,GS_correlation_p_value,Gene_Name,D_melanogaster,I_luminosus,P_pyralis,C_septempunctata,A_verrucosus,T_castaneum,Z_morio,T_molitor,D_ponderosae,R_ferrugineus,A_obtectus,B_siliquastri,C_chinensis,C_maculatus,max_delta_GF,transcript_ID_native,Flybase,Flybase_summary = line.strip().split("\t")
            
            if Orthogroup_ID in orthogroups:
                orthogroups[Orthogroup_ID].append(Gene_Group)
            else:
                orthogroups[Orthogroup_ID] = [Gene_Group]
            
            if Gene_Group not in gene_groups:
                gene_groups[Gene_Group] = Group_function
    
    return orthogroups, gene_groups



def investigate_large_gene_families(tree_path:str, DAVID_table_path:str, orthogroups_path:str, percentile:int = 95, verbose = False, min_GF_size = 2):
    """
    go through the orthogroups by species, and get each that is the largest nth percentile of gene family size in that species
    Then check the DAVID gene groups, and count how these orthogroups are represented there to see if some 
    gene groups are more often expanded in one species than others.
    """

    # get sorted species names from tree
    species_names_unsorted = my_plotting.plot_tree_manually(tree_path)
    species_coords_sorted = sorted(list(species_names_unsorted.keys()))
    species_list = [species_names_unsorted[species_coord] for species_coord in species_coords_sorted]
    orthogroups_GG_dict, gene_groups_function = get_OG_david_annot(DAVID_table_path)
    
    GF_function_counts = {}
    for species in species_list:
        species_dict_list = OGs.parse_orthogroups_dict(orthogroups_path, species=species)
        OGs_species = [OG_id for OG_id,transcript_list in species_dict_list.items() if transcript_list != ['']]
        GF_sizes_list = [ len(transcript_IDs) if transcript_IDs != [''] else 0 for transcript_IDs in species_dict_list.values()]
        sizes = np.array(GF_sizes_list)
        percentile_size = np.percentile(sizes, q = percentile)

        if percentile_size < min_GF_size:
            if verbose:
                print(f"\n   {species} percentile size adjusted from {int(percentile_size)} to {min_GF_size}")
            percentile_size = min_GF_size
        expanded_OGs = { orthogroup : members_list  for orthogroup, members_list in species_dict_list.items() if len(members_list) > percentile_size}
        if verbose:
            print(f" * {species} : number of gene families with more than {int(percentile_size)} members (upper {100-percentile}th percentile) = {len(expanded_OGs)} (of {len(OGs_species)} orthogroups)")

        gene_groups_list = []
        functions = []
        for orthogroup in expanded_OGs.keys():
            try:
                gene_groups = orthogroups_GG_dict[orthogroup]
            except:
                gene_groups = ["unsignificant"] # if the orthogroup is not in the table then it is not significantly rapidly evolving according to CAFE
            gene_groups_list = gene_groups_list + gene_groups
            
            functions_list = [gene_groups_function[gg] if gg != "unsignificant" else "unsignificant" for gg in gene_groups ]
            functions = functions + functions_list
        
        gene_groups_counts = Counter(gene_groups_list)
        functions_conunts = Counter(functions)
                
        if verbose:
            print(f"\t{len(gene_groups_counts)-2} unique Gene Groups with functional annotations, these ones appear more than once:")
            for function, count in functions_conunts.items():
                if count>1 and function not in ["unsignificant", "None"]:
                    print(f"\t  - {count} GFs annotated as : {function}")
        
        GF_function_counts[species] = functions_conunts
        if verbose:
            print(f"")
    return GF_function_counts

        





if __name__ == "__main__":
    
    # out_dir,orthogroups_orthoDB_filepath,tree_path,DAVID_path = filepaths()
    out_dir,orthogroups_orthoDB_filepath,tree_path,DAVID_path = filepaths_work()
    OG_lists_dict = orthogroups_lists()

    ##  IMPORT SVG TO HTML --> makes it so the html file has no external figure dependencies and can be sent by email
    if True:
        # html_path = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/my_thoughts.html"
        html_path = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/my_thoughts.html"
        inline_svgs_in_html(html_path=html_path)

    # --> GENERAL ""ENRICHMENT"" OF GENE GROUP FUNCTION IN RAPIDLY EXPANDING ORTHOGROUPS
    if False:
        investigate_large_gene_families(tree_path=tree_path, DAVID_table_path=DAVID_path, orthogroups_path=orthogroups_orthoDB_filepath, verbose = False)

    # --> AOBT EXPANSION
    OGs_title = " and ".join(OG_lists_dict["Aobt_expansion"])
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Aobt_expansion"], tree_path=tree_path, filename="Aobt_expansion_GF_sizes.png", title = f"A. obtectus expansion: {OGs_title}")

    # --> DETOXIFICATION
    if False:
        cols_list = [
            "#a9c5e2",
            "#434b4c",
            "#DE6449",
            ] # first light blue: "#a9c5e2"
        labels_list = [
            "Cluster 1: Cytochrome P450", 
            "Cluster 3: lipid metabolic process",
            "Cluster 18: aldehyde oxidase", 
            ]
        IDs_lists = [
            OG_lists_dict["Gene Group 1"],
            OG_lists_dict["Gene Group 3"],
            OG_lists_dict["Gene Group 18"],
        ]
        image_path = plot_selected_OGs(
            orthogroups_path=orthogroups_orthoDB_filepath, 
            OG_IDs=IDs_lists, colors=cols_list, labels=labels_list, 
            tree_path=tree_path, filename="detoxificatoin_clusters.png", 
            out_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", 
            title = "Detoxification-related clusters", 
            transparent_bg=True, svg = True)

    # --> REPRODUCTION
    if False:
        cols_list = [
            "#5C3348",
            "#957186",
            "#D9B8C4",
            ] # first light blue: "#a9c5e2"
        labels_list = [
            "Cluster 8: immunity and reproduction", 
            "Cluster 9: sexual reproduction",
            "Cluster 5 and 8: protease inhibitor (immunity and reproduction)", 
            ]
        IDs_lists = [
            OG_lists_dict["Gene Group 8"],
            OG_lists_dict["Gene Group 9"],
            OG_lists_dict["Gene Group 5,8"],
        ]
        image_path = plot_selected_OGs(
            orthogroups_path=orthogroups_orthoDB_filepath, 
            OG_IDs=IDs_lists, colors=cols_list, labels=labels_list, 
            tree_path=tree_path, filename="sexual_reproduction.png", 
            out_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", 
            title = "reproduction and immunity clusters", 
            transparent_bg=True, svg = False)

    # -->  ODORANT BINDING AND PHEROMONE SENSING
    if False:
        cols_list = [
            "#a9c5e2",
            "#91584B",
            "#899D58",
            ] # first light blue: "#a9c5e2"
        labels_list = [
            "Cluster 7: odorant binding", 
            "Cluster 30: pheromone sensing", 
            "Cluster 20: transmembrane transport (in antennae)",
            ]
        IDs_lists = [
            OG_lists_dict["Gene Group 7"],
            OG_lists_dict["Gene Group 30"],
            OG_lists_dict["Gene Group 20"],
        ]
        image_path = plot_selected_OGs(
            orthogroups_path=orthogroups_orthoDB_filepath, 
            OG_IDs=IDs_lists, colors=cols_list, labels=labels_list, 
            tree_path=tree_path, filename="pheromone_sensing_clusters.png", 
            out_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", 
            title = "Detoxification-related clusters", 
            transparent_bg=True, svg = True)

    # --> CHITIN AND CUTICULAR PROTEIN
    if False:
        cols_list = [
            # "#A9C4D9",
            "#331E36",
            "#5C48AD",
            "#7F98C7",
            ] # first light blue: "#a9c5e2"
        labels_list = [
            # "Cluster 26: glycolysis and early development",
            "Cluster 11: Adenosine deaminase-related growth factor",
            "Cluster 15: chitin-related", 
            "Cluster 24: Cuticular protein", 
            ]
        IDs_lists = [
            # OG_lists_dict["Gene Group 26"],
            OG_lists_dict["Gene Group 11"],
            OG_lists_dict["Gene Group 15"],
            OG_lists_dict["Gene Group 24"],
        ]
        image_path = plot_selected_OGs(
            orthogroups_path=orthogroups_orthoDB_filepath, 
            OG_IDs=IDs_lists, colors=cols_list, labels=labels_list, 
            tree_path=tree_path, filename="early_development.png", 
            out_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", 
            title = "chitin formation and Adenosine deaminase-related growth factor", 
            transparent_bg=True, svg = True)

        ## only gene group 26
        cols_list = [
            "#719EC1",
            ] # first light blue: "#a9c5e2"
        labels_list = [
            "Cluster 26: glycolysis and early development",
            ]
        IDs_lists = [
            OG_lists_dict["Gene Group 26"],
        ]
        image_path = plot_selected_OGs(
            orthogroups_path=orthogroups_orthoDB_filepath, 
            OG_IDs=IDs_lists, colors=cols_list, labels=labels_list, 
            tree_path=tree_path, filename="early_development_GF_cluster_26.png", 
            out_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", 
            title = "glycolysis and early development", 
            transparent_bg=True, svg = True)

    # --> ESTERASE AND MATING BEHAVIOUR
    # This group is expanding in elateriformia
    if False:
        cols_list = [
            "#885E5E",
            "#62A87C",
            "#4990C7",
            ] # first light blue: "#a9c5e2"
        labels_list = [
            "Cluster 4: Esterase and mating behavior", 
            "Acyl-CoA synthesis related",
            "Acyl-CoA synthetase", 
            ]
        IDs_lists = [
            OG_lists_dict["Gene Group 4"],
            ["N0.HOG0000284","N0.HOG0000397"],
            ["N0.HOG0000613"]
        ]
        image_path = plot_selected_OGs(
            orthogroups_path=orthogroups_orthoDB_filepath, 
            OG_IDs=IDs_lists, colors=cols_list, labels=labels_list, 
            tree_path=tree_path, filename="fluorescence_elateriformia.png", 
            out_dir = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", 
            title = "Expansions in Elateriformia", 
            transparent_bg=True, svg = False)

