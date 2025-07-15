"""
Plot some stuff to evaluate the functional annotation of the gene groups
"""

import plot_basics as my_plotting
import parse_orthogroups as OGs
import plotting.plot_significant_orthogroups_from_CAFE as plot_OG
import make_orthogroup_flybaseID_table as parse_DAVID
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
    out_dict = {
        "Aobt_expansion" : ["N0.HOG0000035","N0.HOG0000014"],
        "Gene Group 1" : ["N0.HOG0000027","N0.HOG0000059","N0.HOG0000095","N0.HOG0000204","N0.HOG0001077","N0.HOG0000140","N0.HOG0000492","N0.HOG0001030"],
        "Gene Group 3" : ["N0.HOG0000086","N0.HOG0002393","N0.HOG0000085"],
        "Gene Group 17" : ["N0.HOG0000669"],
        "Gene Group 5,8,9" : ["N0.HOG0000541","N0.HOG0000775","N0.HOG0000892","N0.HOG0000401","N0.HOG0009002"],
        "Gene Group 30": ["N0.HOG0000037","N0.HOG0000177","N0.HOG0001445","N0.HOG0000038","N0.HOG0000194","N0.HOG0000345","N0.HOG0000467"],
        "Gene Group 7" : ["N0.HOG0000056","N0.HOG0000454","N0.HOG0000436","N0.HOG0000480","N0.HOG0009039"],
        "Gene Group 16": ["N0.HOG0000307","N0.HOG0001194","N0.HOG0003035"],
        "Gene Group 23": ["N0.HOG0000108","N0.HOG0000039","N0.HOG0000044","N0.HOG0001108","N0.HOG0000067"],
        "Gene Group 4": ["N0.HOG0000120","N0.HOG0000141","N0.HOG0000365","N0.HOG0000378","N0.HOG0007183"],
        "Acyl_CoA_synthesis" : ["N0.HOG0000284","N0.HOG0000397","N0.HOG0000613"],
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
    
    plot_name = f"{out_dir}{filename}"
    plot_OG.plot_gene_counts(OGs_of_interest_dict, sig_list=OG_IDs, all_cafe_list=OG_IDs, 
                             species_names=species_list, annotation = "orthoDB", 
                             filename = plot_name, transparent_bg=False, 
                             title=title, svg=False)
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

    # --> GENERAL ""ENRICHMENT"" OF GENE GROUP FUNCTION IN RAPIDLY EXPANDING ORTHOGROUPS
    if False:
        investigate_large_gene_families(tree_path=tree_path, DAVID_table_path=DAVID_path, orthogroups_path=orthogroups_orthoDB_filepath, verbose = False)

    # --> AOBT EXPANSION
    OGs_title = " and ".join(OG_lists_dict["Aobt_expansion"])
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Aobt_expansion"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Aobt_expansion_GF_sizes.png", title = f"A. obtectus expansion: {OGs_title}")

    # --> DETOXIFICATION
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 1"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_1_detoxofication_GF_sizes.png", title = "Gene group 1")
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 17"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_17_detoxofication_GF_sizes.png", title = "Gene group 17")
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=["N0.HOG0000140"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="detoxofication_N0.HOG0000140_GF_sizes.png", title = "Orthogroup N0.HOG0000140")
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 3"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_3_lipid_metabolism_GF_sizes.png", title = "Gene Group 3")

    # --> REPRODUCTION
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 5,8,9"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_5_8_9_reproduction_GF_sizes.png", title = "Gene group 5, 8, and 9")
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=["N0.HOG0000401"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="OG_N0.HOG0000401_reproduction_GF_sizes.png", title = "Orthogroup N0.HOG0000401")
    
    # -->  ODORANT BINDING AND PHEROMONE SENSING
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 30"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_30_pheromone_sensing_GF_sizes.png", title = "Gene group 30")
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 7"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_7_odorant_binding_GF_sizes.png", title = "Gene group 7")

    # --> CHITIN AND CUTICULAR PROTEIN
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 16"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_16_chitin_related_GF_sizes.png", title = "Gene group 16")
    image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 23"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_23_cuticular_protein_GF_sizes.png", title = "Gene group 23")    

    # --> ESTERASE AND MATING BEHAVIOUR
    # This group is expanding in elateriformia
    if False:
        # latop path
        image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 4"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_4_esterase_GF_sizes.png", title = "Gene group 4")    
    if True:
        # work computer path
        image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 4"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Gene_Group_4_esterase_GF_sizes.png", title = "Gene group 4")    
        image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=["N0.HOG0000613"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="N0_HOG0000613_GF_sizes.png", title = "N0.HOG0000613 (Acyl-CoA synthetase family member 2)")    
        image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Acyl_CoA_synthesis"], tree_path=tree_path, out_dir="/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/", filename="Acyl_CoA_synthesis.png", title = ", ".join(OG_lists_dict["Acyl_CoA_synthesis"]))    

        pass

    ##  IMPORT SVG TO HTML
    if False:
        # html_path = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/my_thoughts.html"
        html_path = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/my_thoughts.html"
        inline_svgs_in_html(html_path=html_path)
