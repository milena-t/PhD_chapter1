# plot the results from compare_gene_positions.sh

import csv

import parse_gff
import matplotlib.pyplot as plt
import numpy as np

# to correctly compute the overlap, first correct the contig names of the orthoDB annotation to match the native annotations with correct_native_contig_names_for_bedtools_intersect.py
# files generated on uppmax are here: /proj/naiss2023-6-65/Milena/annotation_pipeline/annotation_evaluation/gene_position_comparison_native_vs_orhtoDB

def split_at_second_occurrence(s, char = "_"): # split the gene string at the second occurence of "_" to get only the species name
    if s.count(char)<2:
        return s
    else:
        second_occurrence = s.find(char, 2) # start after the first occurence of "_"
        species = s[:second_occurrence]
        return species

# plot the overlap of a single file with the filepath
def plot_gene_overlap_hist(outfile_overlap, filename, reference_annotation, x_label = 'orthDB genes overlapping with one native gene', color = "#4d7298"):
    with open(outfile_overlap, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        no_overlap = [int(row[2]) for row in reader]  # Extract the second column (index 1)

    fig, ax = plt.subplots(1,1, figsize=(12, 12))

    fs = 20
    # customize x axis
    x_max = 15
    x_min = -0.5
    # Create custom ticks: every integer until 10, then every 2nd integer
    ticks_up_to_10 = np.arange(int(x_min+0.5), 11, 1)  # Integers from x_min to 10
    ticks_after_10 = np.arange(10, x_max + 1, 2)  # Every 2nd integer after 10
    custom_ticks = np.unique(np.concatenate((ticks_up_to_10, ticks_after_10)))  # Combine and remove duplicates
    ax.set_xlim(x_min, x_max)
    no_genes = len(no_overlap)
    ax.set_ylim(0,no_genes)
    ax.set_xticks(custom_ticks)
    ax.set_xticklabels(ax.get_xticks(), fontsize=fs)
    ax.set_yticklabels(ax.get_yticks(), fontsize=fs)   

    plt.hist(no_overlap, bins=range(min(no_overlap), max(no_overlap) + 2), color=color, align='left')

    # Add labels and title
    plt.xlabel(x_label, fontsize = fs)
    plt.ylabel('Frequency', fontsize = fs)
    header_name = outfile_overlap.split("/")[-1].split("_gene_overlap_")[0]
    
    plt.title(f'{header_name} gene overlap {reference_annotation} ({no_genes} in ref.)', fontsize = fs)

    # filename = header_name+"_gene_overlap.png"
    plt.savefig(filename, dpi = 300, transparent = True)
    print(f"plot saved in the current working directory as: {filename}")


def plot_all_species_gene_overlap(filenames_list, outfile_overlap_path, plot_filename = "gene_overlap_14_species_transcripts.png", x_label = "how many orthoDB genes overlap with the position of a single native gene (no. all genes)", cols = 3, color = "#4d7298"):
    # plot all species at once in a grid
    cols = cols
    rows = int(len(filenames_list)/cols) # +1
    fig, axes = plt.subplots(rows, cols, figsize=(12, 15))

    print(f"plot in {rows} x {cols} grid")

    # customize x axis
    x_max = 11 # originally 15
    x_min = -0.5
    # Create custom ticks: every integer until 10, then every 2nd integer
    ticks_up_to_10 = np.arange(int(x_min+0.5), 11, 1)  # Integers from x_min to 10
    ticks_after_10 = np.arange(10, x_max + 1, 2)  # Every 2nd integer after 10
    custom_ticks = np.unique(np.concatenate((ticks_up_to_10, ticks_after_10)))  # Combine and remove duplicates


    # Loop over each file path and corresponding subplot axis
    for idx, file_path in enumerate(filenames_list):
        file_path = outfile_overlap_path+file_path
        

        # Calculate row and column indices for the current subplot
        row = idx // cols
        col = idx % cols
        species_name = file_path.split("/")[-1]
        print(f"\tin position {row+1},{col+1}: \t{species_name}")

        with open(file_path, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            no_overlap = [int(row[2]) for row in reader]  # Extract the second column (index 1)
            no_genes = len(no_overlap)

        # Plot histogram on the corresponding subplot axis
        axes[row, col].hist(no_overlap, bins=range(min(no_overlap), max(no_overlap) + 2), color = color, align='left')
        # header_name = file_path.split("/")[-1].split("_gene_overlap_")[0]
        header_name = file_path.split("/")[-1]
        header_name = split_at_second_occurrence(header_name)
        axes[row, col].set_title(f'{header_name} ({no_genes} in ref.)')
        axes[row, col].set_xlabel('')
        axes[row, col].set_ylabel('')

        # unify x axis
        axes[row, col].set_ylim(0, no_genes)
        axes[row, col].set_xlim(x_min, x_max)
        axes[row, col].set_xticks(custom_ticks)


    # Set a single x-axis label for all subplots
    fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=12)
    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # Show the plot
    # plt.show()
    plt.savefig(plot_filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {plot_filename}")

# old files
# outfile_overlap_path = "/Users/miltr339/work/annotation_transcript_overlap_stuff/"

outfile_overlap_path = "/Users/miltr339/work/gene_position_comparison_native_vs_orhtoDB/"
filenames_list_orthoDB_query = [
    "A_obtectus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "A_verrucosus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "B_siliquastri_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "C_analis_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "C_chinensis_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    #"C_maculatus_Lu2024_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "C_maculatus_superscaffolded_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "C_septempunctata_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "D_melanogaster_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "D_ponderosae_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "I_luminosus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "P_pyralis_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "R_ferrugineus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "T_castaneum_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "T_molitor_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "Z_morio_transcript_overlap_stats_numbers_only_orthodb_query.txt"
]

filenames_list_native_query = [
    "A_obtectus_transcript_overlap_stats_numbers_only_native_query.txt",
    "A_verrucosus_transcript_overlap_stats_numbers_only_native_query.txt",
    "B_siliquastri_transcript_overlap_stats_numbers_only_native_query.txt",
    "C_analis_transcript_overlap_stats_numbers_only_native_query.txt",
    "C_chinensis_transcript_overlap_stats_numbers_only_native_query.txt",
    #"C_maculatus_Lu2024_transcript_overlap_stats_numbers_only_native_query.txt",
    "C_maculatus_superscaffolded_transcript_overlap_stats_numbers_only_native_query.txt",
    "C_septempunctata_transcript_overlap_stats_numbers_only_native_query.txt",
    "D_melanogaster_transcript_overlap_stats_numbers_only_native_query.txt",
    "D_ponderosae_transcript_overlap_stats_numbers_only_native_query.txt",
    "I_luminosus_transcript_overlap_stats_numbers_only_native_query.txt",
    "P_pyralis_transcript_overlap_stats_numbers_only_native_query.txt",
    "R_ferrugineus_transcript_overlap_stats_numbers_only_native_query.txt",
    "T_castaneum_transcript_overlap_stats_numbers_only_native_query.txt",
    "T_molitor_transcript_overlap_stats_numbers_only_native_query.txt",
    "Z_morio_transcript_overlap_stats_numbers_only_native_query.txt"
]



plot_all_species_gene_overlap(filenames_list_orthoDB_query, outfile_overlap_path, plot_filename = "transcript_overlap_14_species_orthoDB_query.png", x_label = "number of native transcripts overlaping with one orthoDB transcript (no. transcripts in reference)", color = "#b82946") # red color
plot_all_species_gene_overlap(filenames_list_native_query, outfile_overlap_path, plot_filename = "transcript_overlap_14_species_native_query.png", x_label = "number of orthoDB transcripts overlaping with one native transcript (no. transcripts in reference)", color = "#F2933A") #yellow color

exon_overlaps_dir = "/Users/miltr339/work/gene_position_comparison_native_vs_orhtoDB/exon_overlaps/"

# it's the exons even though it says transcript in the filename
filenames_exon_native_query = [
    "A_obtectus_transcript_overlap_stats_numbers_only_native_query.txt",
    "A_verrucosus_transcript_overlap_stats_numbers_only_native_query.txt",
    "B_siliquastri_transcript_overlap_stats_numbers_only_native_query.txt",
    "C_analis_transcript_overlap_stats_numbers_only_native_query.txt",
    "C_chinensis_transcript_overlap_stats_numbers_only_native_query.txt",
    #"C_maculatus_Lu2024_transcript_overlap_stats_numbers_only_native_query.txt",
    "C_maculatus_superscaffolded_transcript_overlap_stats_numbers_only_native_query.txt",
    #"C_maculatus_transcript_overlap_stats_numbers_only_native_query.txt",
    "C_septempunctata_transcript_overlap_stats_numbers_only_native_query.txt",
    "D_melanogaster_transcript_overlap_stats_numbers_only_native_query.txt",
    "D_ponderosae_transcript_overlap_stats_numbers_only_native_query.txt",
    "I_luminosus_transcript_overlap_stats_numbers_only_native_query.txt",
    "P_pyralis_transcript_overlap_stats_numbers_only_native_query.txt",
    "R_ferrugineus_transcript_overlap_stats_numbers_only_native_query.txt",
    "T_castaneum_transcript_overlap_stats_numbers_only_native_query.txt",
    "T_molitor_transcript_overlap_stats_numbers_only_native_query.txt",
    "Z_morio_transcript_overlap_stats_numbers_only_native_query.txt"
]

filenames_exon_orthoDB_query = [
    "A_obtectus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "A_verrucosus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "B_siliquastri_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "C_analis_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "C_chinensis_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    #"C_maculatus_Lu2024_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "C_maculatus_superscaffolded_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    #"C_maculatus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "C_septempunctata_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "D_melanogaster_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "D_ponderosae_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "I_luminosus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "P_pyralis_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "R_ferrugineus_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "T_castaneum_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "T_molitor_transcript_overlap_stats_numbers_only_orthodb_query.txt",
    "Z_morio_transcript_overlap_stats_numbers_only_orthodb_query.txt"
]

plot_all_species_gene_overlap(filenames_exon_orthoDB_query, exon_overlaps_dir, plot_filename = "exon_overlap_14_species_orthoDB_query.png", x_label = "number of native exons overlaping with one orthoDB exon (no. exons in reference)", color = "#b82946") # red color
plot_all_species_gene_overlap(filenames_exon_native_query, exon_overlaps_dir, plot_filename = "exon_overlap_14_species_native_query.png", x_label = "number of orthoDB exons overlaping with one native exon (no. exons in reference)", color = "#F2933A") #yellow color





outfile_overlap_path_cmac = "/Users/miltr339/work/c_maculatus/annotation_comparison/gene_overlap/"

braker_ref_files = [
    "gene_pos_RNA_combined_vs_native_numbers_only_braker_ref.txt",
    "gene_pos_RNA_simple_vs_native_numbers_only_braker_ref.txt",
    "gene_pos_all_proteinrefs_vs_native_numbers_only_braker_ref.txt",
    "gene_pos_orthoDB_vs_native_numbers_only_braker_ref.txt"
]

native_ref_files = [
    "gene_pos_RNA_combined_vs_native_numbers_only_native_ref.txt",
    "gene_pos_RNA_simple_vs_native_numbers_only_native_ref.txt",
    "gene_pos_all_proteinrefs_vs_native_numbers_only_native_ref.txt",
    "gene_pos_orthoDB_vs_native_numbers_only_native_ref.txt"
]

# plot_all_species_gene_overlap(native_ref_files, outfile_overlap_path_cmac, plot_filename = "gene_overlap_cmac_native_query.png", x_label = "number of orthoDB genes overlaping with one native gene", cols = 2)
# plot_all_species_gene_overlap(braker_ref_files, outfile_overlap_path_cmac, plot_filename = "gene_overlap_cmac_orthodb_query.png", x_label = "number of native genes overlaping with one orthoDB gene", cols = 2)

# plot_gene_overlap_hist(outfile_overlap_path+"C_maculatus_Lu2024_gene_overlap_stats_numbers_only.txt", "C_maculatus_Lu2024_gene_overlap_native_ref.png", "native query", x_label = 'orthDB genes overlapping with one native gene')
# plot_gene_overlap_hist(outfile_overlap_path+"C_maculatus_Lu2024_gene_overlap_stats_numbers_only_orthodb_ref.txt", "C_maculatus_Lu2024_gene_overlap_orthoDB_ref.png", "orthodb query", x_label = 'native genes overlapping with one orthoDB gene')
# plot_gene_overlap_hist(outfile_overlap_path+"T_castaneum_gene_overlap_stats_numbers_only.txt", "T_castaneum_gene_overlap_stats_native_query.png", "native query", x_label = 'orthDB genes overlapping with one native gene', color = "#cd5b5b")
# plot_gene_overlap_hist(outfile_overlap_path+"T_castaneum_gene_overlap_stats_numbers_only_orthodb_ref.txt", "T_castaneum_gene_overlap_stats_orthoDB_query.png", "orthodb query", x_label = 'native genes overlapping with one orthoDB gene')
# plot_gene_overlap_hist(outfile_overlap_path+"T_molitor_gene_overlap_stats_numbers_only.txt", "T_molitor_gene_overlap_stats_native_query.png", "native query", x_label = 'orthDB genes overlapping with one native gene')
# plot_gene_overlap_hist(outfile_overlap_path+"T_molitor_gene_overlap_stats_numbers_only_orthodb_ref.txt", "T_molitor_gene_overlap_stats_orthoDB_query.png", "orthodb query", x_label = 'native genes overlapping with one orthoDB gene')

if False:
    aobt_native = "/Users/milena/work/a_obtectus/a_obtectus_native_isoform_filtered.gff"
    aobt_orthodb = "/Users/milena/work/a_obtectus/a_obtectus_orthodb_isoform_filtered.gff"

    gff_native = parse_gff.parse_gff3_general(aobt_native)
    aobt_single_exon, aobt_multi_exon = parse_gff.get_single_exon_transcripts(gff_native)
    plt.hist(aobt_multi_exon.values())
    plt.hist(aobt_single_exon.values())
    plt.show()
