import unittest
import src.plotting.plot_gene_position_comparison as comp

class TestGeneComparison(unittest.TestCase):
    def test_read_overlaps(self):
        """
        Data for this is not on my laptop! only at work
        """
        overlap_transcript_path = "/Users/miltr339/work/gene_position_comparison_native_vs_orhtoDB/"
        overlap_files_transcript = {
            "A_obtectus" : f"{overlap_transcript_path}A_obtectus_transcript_overlap_stats_complete.txt",
            "A_verrucosus" : f"{overlap_transcript_path}A_verrucosus_transcript_overlap_stats_complete.txt",
            "B_siliquastri" : f"{overlap_transcript_path}B_siliquastri_transcript_overlap_stats_complete.txt",
            "C_analis" : f"{overlap_transcript_path}C_analis_transcript_overlap_stats_complete.txt",
            "C_chinensis" : f"{overlap_transcript_path}C_chinensis_transcript_overlap_stats_complete.txt",
            "C_maculatus" : f"{overlap_transcript_path}C_maculatus_superscaffolded_transcript_overlap_stats_complete.txt",
            "C_septempunctata" : f"{overlap_transcript_path}C_septempunctata_transcript_overlap_stats_complete.txt",
            "D_melanogaster" : f"{overlap_transcript_path}D_melanogaster_transcript_overlap_stats_complete.txt",
            "D_ponderosae" : f"{overlap_transcript_path}D_ponderosae_transcript_overlap_stats_complete.txt",
            "I_luminosus" : f"{overlap_transcript_path}I_luminosus_transcript_overlap_stats_complete.txt",
            "P_pyralis" : f"{overlap_transcript_path}P_pyralis_transcript_overlap_stats_complete.txt",
            "R_ferrugineus" : f"{overlap_transcript_path}R_ferrugineus_transcript_overlap_stats_complete.txt",
            "T_castaneum" : f"{overlap_transcript_path}T_castaneum_transcript_overlap_stats_complete.txt",
            "T_molitor" : f"{overlap_transcript_path}T_molitor_transcript_overlap_stats_complete.txt",
            "Z_morio" : f"{overlap_transcript_path}Z_morio_transcript_overlap_stats_complete.txt",
        }

        for species in overlap_files_transcript.keys():
            self.assertTrue(comp.read_overlaps(overlaps_filepath=overlap_files_transcript[species]))
