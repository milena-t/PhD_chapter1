import unittest
import src.compute_PIC as PIC
from ete3 import Tree

class TEST_PIC(unittest.TestCase):

    def test_documentation(self):
        self.assertTrue(PIC.calculate_PIC.__doc__)
        self.assertTrue(PIC.match_trait_dict_to_tree_leaves.__doc__)
        self.assertTrue(PIC.test_all_leaves_have_traits.__doc__)

    def test_associate_data(self):
        """
        test if the code correctly associates trait values 
        """
        ultrametric_tree= "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_ultrametric.nw"
        species_tree = Tree(ultrametric_tree)
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
        trait_values, discarded_species = PIC.match_trait_dict_to_tree_leaves(species_tree, genome_sizes_dict)
        assert len(discarded_species) == 1 # C_analis is in the genome sizes but not in the species tree so it is discarded

        # mock trees to manually calculate
        newick_str = "((Human:0.2,Chimp:0.2):0.1,Gorilla:0.3);"
        primate_tree = Tree(newick_str)
        traits = {"Human": 5.1,"Chimp": 4.8,"Gorilla": 6.3}
        trait_values, discarded_species = PIC.match_trait_dict_to_tree_leaves(primate_tree, traits)
        assert len(discarded_species) == 0

        

        