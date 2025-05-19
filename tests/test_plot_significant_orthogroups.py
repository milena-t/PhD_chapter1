import unittest
import src.plotting.plot_significant_orthogroups_correlations as cor_plot
import src.plotting.plot_significant_orthogroups_from_CAFE as sig_plot

class TestCAFE_sig_orthogroups(unittest.TestCase):

    def test_CAFE_sig_orthogroups(self):
        """
        Test all functions related to reading the orthofinder output
        """
        orthogroups_native = "/Users/milena/Box Sync/code/CAFE/CAFE_input_native_from_N0.tsv"
        orthogroups_orthoDB = "/Users/milena/Box Sync/code/CAFE/CAFE_input_orthoDB_TE_filtered.tsv"

        orthogroups_input_native = sig_plot.read_orthogroups_input(orthogroups_native)
        orthogroups_input_orthoDB = sig_plot.read_orthogroups_input(orthogroups_orthoDB)
        self.assertTrue(orthogroups_input_native)
        self.assertTrue(orthogroups_input_orthoDB)

        sig_native = "/Users/milena/Box Sync/code/CAFE/native_from_N0_Base_family_results.txt"
        sig_orthoDB = "/Users/milena/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"

        sig_list_n, all_list_n = sig_plot.get_sig_orthogroups(sig_native)
        self.assertTrue(sig_list_n)
        self.assertTrue(all_list_n)
        sig_list, all_list = sig_plot.get_sig_orthogroups(sig_orthoDB)
        self.assertTrue(sig_list)
        self.assertTrue(all_list)

        sig_plot.get_means(orthogroups_input_native, sig_list, all_list)
