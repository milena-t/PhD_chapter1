"""
Plot some stuff to evaluate the functional annotation of the gene groups
"""

import plot_basics as my_plotting
import parse_orthogroups as OGs
import plotting.plot_significant_orthogroups_from_CAFE as plot_OG
import matplotlib.pyplot as plt
import base64
import re
from bs4 import BeautifulSoup


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
        "Gene Group 5,8,9" : ["N0.HOG0000541","N0.HOG0000775","N0.HOG0000892","N0.HOG0000401","N0.HOG0009002"],
        "Gene Group 30": ["N0.HOG0000037","N0.HOG0000177","N0.HOG0001445","N0.HOG0000038","N0.HOG0000194","N0.HOG0000345","N0.HOG0000467"],
        "Gene Group 7" : ["N0.HOG0000056","N0.HOG0000454","N0.HOG0000436","N0.HOG0000480","N0.HOG0009039"],
        "Gene Group 16": ["N0.HOG0000307","N0.HOG0001194","N0.HOG0003035"],
        "Gene Group 23": ["N0.HOG0000108","N0.HOG0000039","N0.HOG0000044","N0.HOG0001108","N0.HOG0000067"]
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
                             title=title, svg=True)
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


if __name__ == "__main__":
    
    out_dir,orthogroups_orthoDB_filepath,tree_path = filepaths()
    OG_lists_dict = orthogroups_lists()

    # --> AOBT EXPANSION
    OGs_title = " and ".join(OG_lists_dict["Aobt_expansion"])
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Aobt_expansion"], tree_path=tree_path, filename="Aobt_expansion_GF_sizes.png", title = f"A. obtectus expansion: {OGs_title}")

    # --> DETOXIFICATION
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 1"], tree_path=tree_path, filename="Gene_Group_1_detoxofication_GF_sizes.png", title = "Gene group 1")
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 17"], tree_path=tree_path, filename="Gene_Group_17_detoxofication_GF_sizes.png", title = "Gene group 17")
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=["N0.HOG0000140"], tree_path=tree_path, filename="detoxofication_N0.HOG0000140_GF_sizes.png", title = "Orthogroup N0.HOG0000140")
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 3"], tree_path=tree_path, filename="Gene_Group_3_lipid_metabolism_GF_sizes.png", title = "Gene Group 3")

    # --> REPRODUCTION
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 5,8,9"], tree_path=tree_path, filename="Gene_Group_5_8_9_reproduction_GF_sizes.png", title = "Gene group 5, 8, and 9")
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=["N0.HOG0000401"], tree_path=tree_path, filename="OG_N0.HOG0000401_reproduction_GF_sizes.png", title = "Orthogroup N0.HOG0000401")
    
    # -->  ODORANT BINDING AND PHEROMONE SENSING
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 30"], tree_path=tree_path, filename="Gene_Group_30_pheromone_sensing_GF_sizes.png", title = "Gene group 30")
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 7"], tree_path=tree_path, filename="Gene_Group_7_odorant_binding_GF_sizes.png", title = "Gene group 7")

    # --> CHITIN AND CUTICULAR PROTEIN
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 16"], tree_path=tree_path, filename="Gene_Group_16_chitin_related_GF_sizes.png", title = "Gene group 16")
    # image_path = plot_selected_OGs(orthogroups_path=orthogroups_orthoDB_filepath, OG_IDs=OG_lists_dict["Gene Group 23"], tree_path=tree_path, filename="Gene_Group_23_cuticular_protein_GF_sizes.png", title = "Gene group 23")    



    ##  IMPORT SVG TO HTML
    html_path = "/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/my_thoughts.html"
    inline_svgs_in_html(html_path=html_path)
