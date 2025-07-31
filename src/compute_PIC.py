""" 
compute phylogenetically independent contrasts from:
 * species tree (ultrametric!)
 * dictionary with { species : trait_value } where the species names correspond to the leaf names in the tree

How PICs work is nicely explained here: https://bio.libretexts.org/Bookshelves/Evolutionary_Developmental_Biology/Phylogenetic_Comparative_Methods_(Harmon)/04%3A_Fitting_Brownian_Motion/4.02%3A_Estimating_Rates_using_Independent_Contrasts

There will always be one less PIC than there are values in the input dict, because the contrasts show the evolutionary change between 
the leaves computed from the trait values and the branch lengths of the tree.
Basically, we compare the traits between two sister species with the rate of evolution between them. PIC generalizes this approach
over the whole tree with the help of a "pruning algorithm" (Felsenstein 1985 and 2004)

"""

from ete3 import Tree
import numpy as np


def match_trait_dict_to_tree_leaves(tree, trait_values):
    """
    check if all trait value keys match all tree leaves
    """
    discarded_species = []
    trait_vaues_matched = {}
    for species in trait_values.keys():
        
        # test that species name appears at all anywhere
        try:
            node = tree.search_nodes(name=species)
        except:
            discarded_species.append(species)

        # test that species name is unique to one leaf
        if len(node) == 1:
            node = node[0]
            if node.is_leaf():
                trait_vaues_matched[species] = trait_values[species]
            else:
                discarded_species.append(species)
        else:
            discarded_species.append(species)
    
    return trait_vaues_matched, discarded_species
        

def calculate_PIC(tree_path:str, trait_values:dict[str,float]):
    """
    calculate phylogenetically independent contrasts of trait values in respect to the tree
    the trait values are in a dictionary where the keys correspond to the leaf names in the tree!
    """
    tree_struct = Tree(tree_path)

    trait_values, discarded_species_list = match_trait_dict_to_tree_leaves(tree_struct, trait_values)

    
"""
TODO continue here:

Chatgpt suggests this kind of recursive traversion but
 1. it's not completely correct
 2. not sure I like it

def compute_pic(node):
    if node.is_leaf():
        return node.value, node.dist

    # Get children
    left, right = node.get_children()

    # Recurse
    x1, v1 = compute_pic(left)
    x2, v2 = compute_pic(right)

    # Calculate independent contrast
    contrast = (x1 - x2) / np.sqrt(v1 + v2)
    node.add_feature("contrast", contrast)

    # Ancestral value (weighted average)
    v_combined = 1 / (1/v1 + 1/v2)
    x_combined = (x1/v1 + x2/v2) * v_combined

    return x_combined, v_combined + node.dist


"""


if __name__ == "__main__":
    
    ultrametric_tree = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_ultrametric.nw"
    test_OG = {
        "D_melanogaster" : 16,
        "I_luminosus" : 73,
        "P_pyralis" : 28,
        "C_septempunctata" : 10,
        "A_verrucosus" : 29,
        "T_castaneum" : 35,
        "Z_morio" : 2,
        "T_molitor" : 52,
        "D_ponderosae" : 8,
        "R_ferrugineus" : 18,
        "A_obtectus" : 9,
        "B_siliquastri" : 14,
        "C_chinensis" : 10,
        "C_maculatus" : 12,
    }
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
    	
    newick_str = "((Human:0.2,Chimp:0.2):0.1,Gorilla:0.3);"
    traits = {"Human": 5.1,"Chimp": 4.8,"Gorilla": 6.3}

    calculate_PIC(tree_path=ultrametric_tree, trait_values=test_OG)
    	
    	
    	
    	
    	
    	
    	
    	
    	
    