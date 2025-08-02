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
                print(f"  !!! {species} ({trait_values[species]}) is not present in the tree.")
        else:
            discarded_species.append(species)
            print(f"  !!! {species} ({trait_values[species]}) is not present in the tree.")
    
    return trait_vaues_matched, discarded_species


def test_all_leaves_have_traits(tree, trait_values):
    """
    test if all leaves in the tree have a trait value associated
    Error if not!
    """
    trait_names = list(trait_values.keys())

    for node in tree.traverse("postorder"):
        if node.is_leaf() == False:
            continue
        elif node.is_leaf() and node.name in trait_names:
            # node.add_feature(value = trait_values[node.name])
            node.add_feature("value", trait_values[node.name])
        else: # node is leaf and node.name not in trait_names
            node.delete()
            print(f"  !!! {node.name} is a tree leaf with no assigned trait --> delete from tree")
    return tree

        

def calculate_PIC(tree_path:str, trait_values:dict[str,float], verbose = False, max_iterations = 0):
    """
    calculate phylogenetically independent contrasts of trait values in respect to the tree
    the trait values are in a dictionary where the keys correspond to the leaf names in the tree!
    """
    if max_iterations == 0:
        max_iterations = len(trait_values)**2

    tree = Tree(tree_path)

    ##  prepare data: remove trait values if they're not in the tree and remove tree leaves if they don't have trait values
    trait_values, discarded_species_list = match_trait_dict_to_tree_leaves(tree, trait_values)
    tree = test_all_leaves_have_traits(tree, trait_values)

    PICs = []
    ## Traverse through tree from the leaves up
    node_number = 1

    while len(PICs) < len(trait_values)-1:
        for node in tree.traverse("postorder"):

            if len(node.children)== 2:
                children = node.children
                child1 = children[0]
                child2 = children[1]

                if child1.is_leaf() and child2.is_leaf():

                    # calculate standardized contrast between two leaves
                    c = child1.value - child2.value
                    sum_dist = np.sqrt(child1.dist + child2.dist)
                    s = c/sum_dist

                    PICs.append(s)

                    # prune leaves and update MRCA to become a new leaf
                    node_dist = 1 / (1/child1.dist + 1/child2.dist)
                    node_value = (child1.value/child1.dist + child2.value/child2.dist) * node_dist
                    # node_dist = node.dist + child1.dist*child2.dist/(child1.dist+child2.dist)
                    # node_value = (child1.value/child1.dist + child2.value/child2.dist) / (1/child1.dist + 1/child2.dist)

                    # node.add_feature(value = node_value)
                    node.add_feature("value", node_value)
                    node.dist = node_dist
                    node.name = f"anc({child1.name},{child2.name})"

                    if verbose:
                        print(f"s: {s:.4}\t for {child1.name} and {child2.name}")

                    node.remove_child(child1)
                    node.remove_child(child2)
                    # print(f" after deleting children is node {node} a leaf: {node.is_leaf()}")
                    # print(tree)

        node_number += 1
        if node_number == max_iterations:
            raise RuntimeError(f"time out while tree traversing, {node_number} loops")
    
    return PICs
    
    
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

    ## check against R calculations:
    # /Box Sync/code/annotation_pipeline/repeatmasking_eval/gene_numbers_stats.Rmd
    
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
    
    # coleoptera_PICs = calculate_PIC(tree_path=ultrametric_tree, trait_values=genome_sizes_dict, verbose=True)
    # coleoptera_PICs.sort()
    # print(coleoptera_PICs)
    # 	
    # newick_str = "((Human:0.2,Chimp:0.2):0.1,Gorilla:0.3);"
    # traits = {"Human": 5.1,"Chimp": 4.8,"Gorilla": 6.3}
    # R results: -1.9091883  0.4743416

    primates_tree2 = "((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);"
    primates_traits_a2 = {"Homo":4.09434, "Pongo":3.61092,"Macaca":2.37024, "Ateles":2.02815, "Galago":-1.46968}
    # R results: 0.7459333 1.1929263 1.5847416 3.3583189 
    primates_traits_b2 = {"Homo":4.74493,"Pongo":3.33220,"Macaca":3.36730, "Ateles":2.89037, "Galago":2.30259}
    # R results: 0.7176125 0.8678969 0.8970604 2.1798897 

    primate_PICs = calculate_PIC(primates_tree2, primates_traits_a2)
    primate_PICs.sort()
    print(primate_PICs)
    
    	
    	
    	
    	
    	
    	
    	
    	
    	
    