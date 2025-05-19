"""Tests parsing an annotation file with code in code in src.parse_gff"""
import unittest

import src.parse_gff as gff


class TestFeatureClass(unittest.TestCase):

    def test_Feature(self):
        """The class 'Feature' returns an object of type feature."""
        # self.assertTrue(Feature.__doc__)
        self.assertEqual(str(type(gff.Feature("1", "contig1", gff.FeatureCategory.Gene, 1, 100, "+", "."))), "<class 'src.parse_gff.Feature'>")
        self.assertIsInstance(gff.Feature("1", "contig1", gff.FeatureCategory.Gene, 1, 100, "+", "."), gff.Feature)

        
    def test_read_annotation(self):
        """All annotation files are read correctly into the default structure"""
        native_annot_dir = "/Users/milena/work/native_annot_gff/"
        native_annotations = {
            "A_obtectus" : f"{native_annot_dir}acanthoscelides_obtectus.gff.isoform_filtered",
            "A_verrucosus" : f"{native_annot_dir}asbolus_verrucosus.gff.isoform_filtered",
            "B_siliquastri" : f"{native_annot_dir}bruchidius_siliquastri.gff.isoform_filtered",
            "C_analis" : f"{native_annot_dir}callosobruchus_analis.gff.isoform_filtered",
            "C_chinensis" : f"{native_annot_dir}callosobruchus_chinensis.gff.isoform_filtered",
            "C_maculatus" : f"{native_annot_dir}callosobruchus_maculatus.gff.isoform_filtered",
            "C_tempunctata" : f"{native_annot_dir}coccinella_septempunctata.gff.isoform_filtered",
            "D_ponderosae" : f"{native_annot_dir}dendroctonus_ponderosae.gff.isoform_filtered",
            "D_melanogaster" : f"{native_annot_dir}drosophila_melanogaster.gff.isoform_filtered",
            "I_luminosus" : f"{native_annot_dir}ignelater_luminosus.gff.isoform_filtered",
            "I_luminosus_mod" : f"{native_annot_dir}ignelater_luminosus.gff.isoform_filtered_modified",
            "P_pyralis" : f"{native_annot_dir}photinus_pyralis.gff.isoform_filtered",
            "rferrugineus" : f"{native_annot_dir}rhynchophorus_ferrugineus.gff.isoform_filtered",
            "T_molitor" : f"{native_annot_dir}tenebrio_molitor.gff.isoform_filtered",
            "T_castaneum" : f"{native_annot_dir}tribolium_castaneum.gff.isoform_filtered",
            "Z_morio" : f"{native_annot_dir}zophobas_morio.gff.isoform_filtered",
        }
        for species in native_annotations.keys():
            self.assertTrue(gff.parse_gff3_general(native_annotations[species], verbose = False))

    def test_read_annotation_by_contig(self):    
        """All annotation files can also be read by-contig"""
        native_annot_dir = "/Users/milena/work/native_annot_gff/"
        native_annotations = {
            "A_obtectus" : f"{native_annot_dir}acanthoscelides_obtectus.gff.isoform_filtered",
            "A_verrucosus" : f"{native_annot_dir}asbolus_verrucosus.gff.isoform_filtered",
            "B_siliquastri" : f"{native_annot_dir}bruchidius_siliquastri.gff.isoform_filtered",
            "C_analis" : f"{native_annot_dir}callosobruchus_analis.gff.isoform_filtered",
            "C_chinensis" : f"{native_annot_dir}callosobruchus_chinensis.gff.isoform_filtered",
            "C_maculatus" : f"{native_annot_dir}callosobruchus_maculatus.gff.isoform_filtered",
            "C_tempunctata" : f"{native_annot_dir}coccinella_septempunctata.gff.isoform_filtered",
            "D_ponderosae" : f"{native_annot_dir}dendroctonus_ponderosae.gff.isoform_filtered",
            "D_melanogaster" : f"{native_annot_dir}drosophila_melanogaster.gff.isoform_filtered",
            "I_luminosus" : f"{native_annot_dir}ignelater_luminosus.gff.isoform_filtered",
            "I_luminosus_mod" : f"{native_annot_dir}ignelater_luminosus.gff.isoform_filtered_modified",
            "P_pyralis" : f"{native_annot_dir}photinus_pyralis.gff.isoform_filtered",
            "rferrugineus" : f"{native_annot_dir}rhynchophorus_ferrugineus.gff.isoform_filtered",
            "T_molitor" : f"{native_annot_dir}tenebrio_molitor.gff.isoform_filtered",
            "T_castaneum" : f"{native_annot_dir}tribolium_castaneum.gff.isoform_filtered",
            "Z_morio" : f"{native_annot_dir}zophobas_morio.gff.isoform_filtered",
        }

        for species in native_annotations.keys():
            self.assertTrue(gff.parse_gff3_by_contig(native_annotations[species], verbose = False))


    # def test_is_zero_responds_correctly_to_ints(self):
    #     """The function 'is_zero' responds correctly to integers."""
    #     self.assertTrue(is_zero(0))
    #     self.assertFalse(is_zero(1))

    # def test_is_zero_raises_an_exception_upon_non_ints(self):
    #     """The function 'is_zero' raises an exception upon non-ints."""
    #     self.assertRaises(TypeError, is_zero, {1, 2})
    #     self.assertRaises(TypeError, is_zero, "I am a string")
        