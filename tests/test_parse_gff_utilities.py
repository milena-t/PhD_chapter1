"""Tests utility functions that rely on parsed gff files in in src.parse_gff"""
import unittest

import src.parse_gff as gff


class TestFeatureClass(unittest.TestCase):

    def test_general_utilities(self):
        self.assertEqual(gff.split_at_second_occurrence("Species_name_file_description.extension"), "Species_name")

    def test_single_exon_transcripts(self):
        """
        Test all the code used to get single-exon transcripts
        """
        
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

        for species in native_annotations.values():
            single_exon, multi_exon = gff.get_single_exon_transcripts(gff.parse_gff3_general(species))
            self.assertTrue(single_exon)
            self.assertTrue(multi_exon)


    def test_get_transcript_lengths(self):
        """
        Test all the code used to get transcript lengths
        """
        
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

        for species in native_annotations.values():
            self.assertTrue(gff.get_transcript_lengths(gff.parse_gff3_general(species, verbose=False)))