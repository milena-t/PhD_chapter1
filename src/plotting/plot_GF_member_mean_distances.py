"""
Analyse the relative positions of all transcripts in an orthogroup to try and infer the transcript mechanism.
Assumed that all paralogs that are in close proximity might be tandem duplicated, and dispersed paralogs are
duplicated through different mechanisms
"""

import parse_gff as gff
import parse_orthogroups as OGs
from statistics import mean
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


def filepaths_native():
    native_annot_dir = "/Users/miltr339/work/native_annotations/all_native_annot/"
    native_annotations = {
        "A_obtectus" : f"{native_annot_dir}A_obtectus_annotation_isoform_filtered.gff",
        "A_verrucosus" : f"{native_annot_dir}A_verrucosus_annotation_isoform_filtered.gff",
        "B_siliquastri" : f"{native_annot_dir}B_siliquastri_annotation_isoform_filtered.gff",
        # "C_analis" : f"{native_annot_dir}C_analis_annotation_isoform_filtered.gff",
        "C_chinensis" : f"{native_annot_dir}C_chinensis_annotation_isoform_filtered.gff",
        "C_maculatus" : f"{native_annot_dir}C_maculatus_superscaffolded_liftover_annotation.gff",
        "C_septempunctata" : f"{native_annot_dir}C_septempunctata_annotation_isoform_filtered.gff",
        "D_melanogaster" : f"{native_annot_dir}D_melanogaster_annotation_isoform_filtered.gff",
        "D_ponderosae" : f"{native_annot_dir}D_ponderosae_annotation_isoform_filtered.gff",
        "I_luminosus" : f"{native_annot_dir}I_luminosus_annotation_isoform_filtered.gff",
        "P_pyralis" : f"{native_annot_dir}P_pyralis_annotation_isoform_filtered.gff",
        "R_ferrugineus" : f"{native_annot_dir}R_ferrugineus_annotation_isoform_filtered.gff",
        "T_castaneum" : f"{native_annot_dir}T_castaneum_annotation_isoform_filtered.gff",
        "T_molitor" : f"{native_annot_dir}T_molitor_annotation_isoform_filtered.gff",
        "Z_morio" : f"{native_annot_dir}Z_morio_annotation_isoform_filtered.gff",
    }
    
    native_proteinseqs_dir = "/Users/miltr339/work/native_proteinseqs/"
    native_proteinseqs={
        "A_obtectus" : f"{native_proteinseqs_dir}A_obtectus.faa",
        "A_verrucosus" : f"{native_proteinseqs_dir}A_verrucosus.faa",
        "B_siliquastri" : f"{native_proteinseqs_dir}B_siliquastri.faa",
        "C_chinensis" : f"{native_proteinseqs_dir}C_chinensis.faa",
        "C_maculatus" : f"{native_proteinseqs_dir}C_maculatus.faa",
        "C_septempunctata" : f"{native_proteinseqs_dir}C_septempunctata.faa",
        "D_melanogaster" : f"{native_proteinseqs_dir}D_melanogaster.faa",
        "D_ponderosae" : f"{native_proteinseqs_dir}D_ponderosae.faa",
        "I_luminosus" : f"{native_proteinseqs_dir}I_luminosus.faa",
        "P_pyralis" : f"{native_proteinseqs_dir}P_pyralis.faa",
        "R_ferrugineus" : f"{native_proteinseqs_dir}R_ferrugineus.faa",
        "T_castaneum" : f"{native_proteinseqs_dir}T_castaneum.faa",
        "T_molitor" : f"{native_proteinseqs_dir}T_molitor.faa",
        "Z_morio" : f"{native_proteinseqs_dir}Z_morio.faa",
    }

    # orthogroups_native = "/Users/miltr339/work/orthofinder/native_orthogroups/N0.tsv"
    # orthogroups_native = "/Users/milena/work/orthofinder/native_orthogroups/N0.tsv"
    # sig_native = "/Users/miltr339/Box Sync/code/CAFE/native_from_N0_Base_family_results.txt"
    orthogroups_native = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/N0.tsv"
    sig_native = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_native_Base_Family_results.txt"

    return native_annotations,orthogroups_native,sig_native,native_proteinseqs

def filepaths_orthoDB():
    orthoDB_annot_dir = "/Users/miltr339/work/orthoDB_annotations/"
    orthoDB_annotations = {
        "A_obtectus" : f"{orthoDB_annot_dir}A_obtectus_braker_isoform_filtered.gff",
        "A_verrucosus" : f"{orthoDB_annot_dir}A_verrucosus_braker_isoform_filtered.gff",
        "B_siliquastri" : f"{orthoDB_annot_dir}B_siliquastri_braker_isoform_filtered.gff",
        # "C_analis" : f"{orthoDB_annot_dir}C_analis_braker_isoform_filtered.gff",
        "C_chinensis" : f"{orthoDB_annot_dir}C_chinensis_braker_isoform_filtered.gff",
        "C_maculatus" : f"{orthoDB_annot_dir}C_maculatus_superscaffolded_annotation_isoform_filtered.gff",
        "C_septempunctata" : f"{orthoDB_annot_dir}C_septempunctata_braker_isoform_filtered.gff",
        "D_melanogaster" : f"{orthoDB_annot_dir}D_melanogaster_braker_isoform_filtered.gff",
        "D_ponderosae" : f"{orthoDB_annot_dir}D_ponderosae_braker_isoform_filtered.gff",
        "I_luminosus" : f"{orthoDB_annot_dir}I_luminosus_braker_isoform_filtered.gff",
        "P_pyralis" : f"{orthoDB_annot_dir}P_pyralis_braker_isoform_filtered.gff",
        "R_ferrugineus" : f"{orthoDB_annot_dir}R_ferrugineus_braker_isoform_filtered.gff",
        "T_castaneum" : f"{orthoDB_annot_dir}T_castaneum_braker_isoform_filtered.gff",
        "T_molitor" : f"{orthoDB_annot_dir}T_molitor_braker_isoform_filtered.gff",
        "Z_morio" : f"{orthoDB_annot_dir}Z_morio_braker_isoform_filtered.gff",
    }

    orthoDB_proteinseqs_dir = "/Users/miltr339/work/orthoDB_proteinseqs_TE_filtered/"
    orthoDB_proteinseqs = {
        "A_obtectus" : f"{orthoDB_proteinseqs_dir}A_obtectus_filtered_proteinfasta_TE_filtered.fa",
        "A_verrucosus" : f"{orthoDB_proteinseqs_dir}A_verrucosus_filtered_proteinfasta_TE_filtered.fa",
        "B_siliquastri" : f"{orthoDB_proteinseqs_dir}B_siliquastri_filtered_proteinfasta_TE_filtered.fa",
        # "C_analis" : f"{orthoDB_proteinseqs_dir}C_analis_filtered_proteinfasta_TE_filtered.fa",
        "C_chinensis" : f"{orthoDB_proteinseqs_dir}C_chinensis_filtered_proteinfasta_TE_filtered.fa",
        "C_maculatus" : f"{orthoDB_proteinseqs_dir}C_maculatus_filtered_proteinfasta_TE_filtered.fa",
        "C_septempunctata" : f"{orthoDB_proteinseqs_dir}C_septempunctata_filtered_proteinfasta_TE_filtered.fa",
        "D_melanogaster" : f"{orthoDB_proteinseqs_dir}D_melanogaster_filtered_proteinfasta_TE_filtered.fa",
        "D_ponderosae" : f"{orthoDB_proteinseqs_dir}D_ponderosae_filtered_proteinfasta_TE_filtered.fa",
        "I_luminosus" : f"{orthoDB_proteinseqs_dir}I_luminosus_filtered_proteinfasta_TE_filtered.fa",
        "P_pyralis" : f"{orthoDB_proteinseqs_dir}P_pyralis_filtered_proteinfasta_TE_filtered.fa",
        "R_ferrugineus" : f"{orthoDB_proteinseqs_dir}R_ferrugineus_filtered_proteinfasta_TE_filtered.fa",
        "T_castaneum" : f"{orthoDB_proteinseqs_dir}T_castaneum_filtered_proteinfasta_TE_filtered.fa",
        "T_molitor" : f"{orthoDB_proteinseqs_dir}T_molitor_filtered_proteinfasta_TE_filtered.fa",
        "Z_morio" : f"{orthoDB_proteinseqs_dir}Z_morio_filtered_proteinfasta_TE_filtered.fa",
    }

    # orthogroups_orthoDB = "/Users/miltr339/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # orthogroups_orthoDB = "/Users/milena/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # sig_orthoDB = "/Users/miltr339/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"
    
    orthogroups_orthoDB = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv"
    sig_orthoDB = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_uniform_Base_Family_results.txt"

    return orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs


def analyze_transcript_positions(transcripts_list:list, annotation_dict:dict):
    """
    The transcript_list contains transcript IDs of all members of a gene family in one species,
    the annotation_dict is the parsed annotation of that species.

    This function then analyzes the transcript positions:
        * How many orthogroups are confined to one contig vs. distributed over multiple contigs? 
          (keep in mind assembly contiguity! sometimes contig == chromosome, and sometimes not)
        * The mean distance of transcripts in an orthogroup if they are confined to the same transcript
          (the distances are calculated middle-to-middle. Since the transcripts are different lengths, I
          calculate mean transcript length and subtract that from the mean distance to get an estimation 
          of mean inter-transcript distance)
    """
    same_contig = True
    middle_positions = [] # middle positions of each transcript.
    transcript_lengths = [] # transcript lengths of all transcripts in the orthogroup
    not_found_transcripts = []
    
    if len(transcripts_list) <= 1:
        mean_distance = 0
        return same_contig, mean_distance

    try:
        first_contig = annotation_dict[transcripts_list[0]].contig
    except:
        ## if the _1 suffix from the orthofinder transcripts makes problems remove it here
        transcripts_list = [transcript[:-2] for transcript in transcripts_list]
        try:
            first_contig = annotation_dict[transcripts_list[0]].contig
        except:
            same_contig = False
            mean_distance = transcripts_list[0]
            return same_contig, mean_distance

    for transcript in transcripts_list:
        try:
            transcript = annotation_dict[transcript]
        except:
            not_found_transcripts.append(transcript)
            same_contig = False
            mean_distance = transcript
            return same_contig, mean_distance

        if first_contig != transcript.contig:
            same_contig = False
            mean_distance = 0
            return same_contig, mean_distance

        tr_start = transcript.start
        if transcript.end < transcript.start:
            tr_start = transcript.end
        
        transcript_length = abs(transcript.start - transcript.end)
        transcript_lengths.append(transcript_length)
        half_length = transcript_length*0.5
        middle_positions.append(tr_start+half_length)

    mean_middle_distance = int(mean(middle_positions))
    mean_length = int(mean(transcript_lengths))
    mean_distance = mean_middle_distance - mean_length

    if mean_distance<0:
        raise ValueError(f"the mean distance of transcripts is {mean_distance}, which cannot be < 0")

    return same_contig, mean_distance



def analyze_orthogroup_position_species(OGs_dict, species_annotation):
    """
    Calculate species-wide orthogroups statistics:
        * percent of orthogroups where not all gene family members are on the same contig
        * dict:
            {
                orthogroup_ID : [ num_members, mean_distance ] , 
                orthogroup_ID : [ num_members, mean_distance ] , 
                ...
            }
    """
    all_OGs = 0
    same_contig_OGs = 0
    out_dict = {}
    orthogroups_with_not_found_transcripts = 0
    
    for orthogroup, transcripts_list in OGs_dict.items():
        # if len(transcripts_list)>2:
        #     print(transcripts_list)
        #     raise RuntimeError("testing")
        same_contig, mean_distance = analyze_transcript_positions(transcripts_list, species_annotation)
        all_OGs += 1
        if same_contig:
            same_contig_OGs += 1
        if type(mean_distance) == str:
            orthogroups_with_not_found_transcripts += 1
        if len(transcripts_list)>1 and same_contig:
                out_dict[orthogroup] = [len(transcripts_list), mean_distance]
    
    same_contig_proportion = float(same_contig_OGs) / float(all_OGs)

    if orthogroups_with_not_found_transcripts >0:
        print(f"{orthogroups_with_not_found_transcripts} (of {all_OGs}) orthogroups have transcripts that were not found in the annotation")

    return same_contig_proportion, out_dict


def plot_transcript_distance(same_contig_proportion, GF_positions_dict, species):
    """
    plot the output from analyze_orthogroup_position_species
    """

    filename = f"mean_transcript_distance_in_gene_families_{species}.png"

    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    fs = 15

    nun_members_vec = []
    mean_distance_vec = []
    
    for values in GF_positions_dict.values():
        nun_members_vec.append(values[0])
        mean_distance_vec.append(values[1])
    
    if "_" in species:
        species = species.replace("_", ". ")
    
    percent = int(same_contig_proportion*100)
    plt.title(f"{species} : {percent}% of gene families have all members on the same contig", fontsize = fs)
    ax.scatter(nun_members_vec, mean_distance_vec, color = "#8E8E8E")

    ylab = f"mean distance between transcripts in the gene family (bp)"
    ax.set_ylabel(ylab, fontsize = fs)
    xlab = f"number of gene family members"
    ax.set_xlabel(xlab, fontsize = fs)

    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x < 0 else f'{x / 1e6:.0f} Mb'))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x < 0 else f"{int(x)}"))
    ax.tick_params(axis ='y', labelsize = fs)  
    ax.tick_params(axis ='x', labelsize = fs)  
    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = False)
    print(f"plot saved in current working directory as {filename}")


def plot_all_OGs_transcript_distances(same_contig_proportion_all_species, GF_positions_dict_all_species, columns = 3, filename = "mean_transcript_distance_in_gene_families.png", L50_values = {}, N50_values = {}):
    species_list = list(same_contig_proportions.keys())
    cols = columns
    rows = int(len(species_list)/cols)  +1
    fig, axes = plt.subplots(rows, cols, figsize=(12, 17))
    fs = 18

    for idx, species in enumerate(species_list):
        same_contig_proportion = same_contig_proportion_all_species[species]
        percent = int(same_contig_proportion*100)
        GF_positions_dict = GF_positions_dict_all_species[species]

        nun_members_vec = []
        mean_distance_vec = []
        for values in GF_positions_dict.values():
            nun_members_vec.append(values[0])
            mean_distance_vec.append(values[1])
        
        # Calculate row and column indices for the current subplot
        row = idx // cols
        col = idx % cols
        species_name = species.replace("_", ". ")

        # Plot histogram on the corresponding subplot axis
        axes[row, col].scatter(nun_members_vec, mean_distance_vec, color = "#8E8E8E")
        try:
            N50_value = N50_values[species]
            if N50_value>10e6:
                N50_value = f'{N50_value / 1e6:.0f} Mb'
            elif N50_value>1e6:
                N50_value = f'{N50_value / 1e6:.1f} Mb'
            else:
                N50_value = f'{N50_value / 1e3:.0f} kb'
            axes[row, col].set_title(f'{species_name} \n({percent}% of gene families) \nL50: {L50_values[species]}, N50: {N50_value}', fontsize = fs)
        except:
            axes[row, col].set_title(f'{species_name} ({percent}% of gene families)', fontsize = fs)
        axes[row, col].set_xlabel('')
        axes[row, col].set_ylabel('')
        
        y_lim = max(mean_distance_vec)
        y_scale = ""
        if y_lim>10e6:
            y_scale= "1.0Mb"
            y_function_formatter = lambda x, pos: '' if x <= 0 else f'{x / 1e6:.0f} Mb'
        elif y_lim>1e6:
            y_scale= "0.1Mb"
            y_function_formatter = lambda x, pos: '' if x <= 0 else f'{x / 1e6:.1f} Mb'
        else:
            y_scale= "1.0kb"
            y_function_formatter = lambda x, pos: '' if x <= 0 else f'{x / 1e3:.0f} kb'

        axes[row, col].yaxis.set_major_formatter(FuncFormatter(y_function_formatter))
        # make integers and skip ticks that are in between integers
        axes[row, col].xaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if str(x).split(".")[-1][0]!= "0" else f"{int(x)}"))
        axes[row, col].tick_params(axis='y', labelsize=fs)
        axes[row, col].tick_params(axis='x', labelsize=fs)
        print(f"\tin position {row+1},{col+1} ;  y scale: {y_scale}  --> {species_name}")
    
    # make last plot empty
    idx_max = len(species_list)
    row = idx_max // cols
    col = idx_max % cols
    axes[row, col].axis('off')
    axes[row, col].set_title(f'')


            # Set a single x-axis label for all subplots
    x_label = f"number of gene family members"
    fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=fs)
    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    plt.savefig(filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {filename}")



if __name__ == "__main__":
    
    
    orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs = filepaths_orthoDB()
    native_annotations, orthogroups_native, sig_native, native_proteinseqs = filepaths_native()
    
    # from quast results for all assemblies
    L50_values = {
        "A_obtectus" : 5,
        "A_verrucosus" : 8944,
        "B_siliquastri" : 5,
        "C_analis" : 1051,
        "C_chinensis" : 231,
        "C_maculatus" : 6, # the non-superscaffolded value is 38,
        "C_septempunctata" : 4,
        "D_melanogaster" : 3,
        "D_ponderosae" : 87,
        "I_luminosus" : 1966,
        "P_pyralis" : 5,
        "R_ferrugineus" : 268,
        "T_castaneum" : 5,
        "T_molitor" : 6,
        "Z_morio" : 4
    }
    N50_values = {
        "A_obtectus" : 108704056,
        "A_verrucosus" : 6513,
        "B_siliquastri" : 39857073, ##!
        "C_analis" : 240251,
        "C_chinensis" : 800721,
        "C_maculatus" : 101463975,
        "C_septempunctata" : 41442133, ##!
        "D_melanogaster" : 25286936, ##!
        "D_ponderosae" : 16386110,
        "I_luminosus" : 95792,
        "P_pyralis" : 47017841,
        "R_ferrugineus" : 471583,
        "T_castaneum" : 20636470, ##!
        "T_molitor" : 20827499,
        "Z_morio" : 48007226
    }

    # species = "B_siliquastri"
    same_contig_proportions:dict = {}
    gene_family_values:dict = {}

    for species in orthoDB_annotations.keys():

        print(f"\n\t---> {species}")
        sig_list, all_list = OGs.get_sig_orthogroups(sig_orthoDB)
        species_OGs_dict = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_list=sig_list, species=species)
        # for orthogroup, transcripts in species_OGs_dict.items():
        #     print(f"{orthogroup} : {transcripts}")

        species_annotation = gff.parse_gff3_general(orthoDB_annotations[species], verbose=False, keep_feature_category=gff.FeatureCategory.Transcript)

        same_contig_proportion, out_dict = analyze_orthogroup_position_species(species_OGs_dict, species_annotation)

        percent = int(same_contig_proportion*100)
        print(f"{species} : {percent}% of gene families have all members on the same contig")

        same_contig_proportions[species] = same_contig_proportion
        gene_family_values[species] = out_dict

        # for OG_id, out_values in out_dict.items():
        #     print(f"{OG_id} :  num. members {out_values[0]}, mean distance = {out_values[1]}")
        
        # plot_transcript_distance(same_contig_proportion, out_dict, species)

    plot_all_OGs_transcript_distances(same_contig_proportion_all_species=same_contig_proportions, GF_positions_dict_all_species=gene_family_values, filename="/Users/miltr339/work/PhD_code/PhD_chapter1/data/mean_transcript_distance_in_gene_families.png", L50_values = L50_values, N50_values = N50_values)