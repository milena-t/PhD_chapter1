"""
Associate the surroundings of genes that are part of rapidly evolving gene families
for TE presence via a stacked barplot. Each bar is a base, and the coverage is how many TEs of each class are
annotated on that base
"""

import numpy as np
import parse_gff as gff
import parse_repeats as repeats
import parse_orthogroups as OGs
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from math import ceil

from tqdm import tqdm

def get_all_repeat_categories(repeats_dict):
    all_contigs = []
    for contig in repeats_dict.keys():
        all_contigs.append(list(set([repeat_instance.repeat_category for repeat_instance in repeats_dict[contig]])))
    all_types = [rep_type for contig_list in all_contigs for rep_type in contig_list]
    all_types = list(set(all_types))
    return(all_types)


def get_sig_transcripts(orthoDB_orthogroups):
    """
    get a list of all transcripts in the orthoDB dictionary (presumeably significant transcripts)
    """
    all_transcript_IDs = []
    for transcripts_list in orthoDB_orthogroups.values():
        for transcript_id in transcripts_list:
            transcript_id = transcript_id[:-2] # remove the "_1" suffix
            if len(transcript_id)>0:
                all_transcript_IDs.append(transcript_id)
    return all_transcript_IDs


def filter_sig_OGs_by_size(orthoDB_orthogroups:dict, species:str, q:int, verbose=False):
    """
    return a list of orthogroup IDs, where the GF size in the species is above the q'th percentile
    """
    GF_sizes_species = {}
    sizes = []
    for OG_id, transcripts_list in orthoDB_orthogroups.items():
        GF_sizes_species[OG_id] = len(transcripts_list)
        sizes.append(len(transcripts_list))

    ## calculate size threshold
    OGs_filtered = []
    sizes = np.array(sizes)
    percentile_size = np.percentile(sizes, q = q)
    if verbose:
        print(f" ---> {species}: \n\tmax gene family: {max(sizes)} \n\tmin gene family: {min(sizes)}\n--> {q}th percentile : {percentile_size}")

    for OG_id, size in GF_sizes_species.items():
        if size > percentile_size:
            OGs_filtered.append(OG_id)
    if verbose:
        print(f"before filtering: {len(orthoDB_orthogroups)} --> after filtering: {len(OGs_filtered)}")

    return OGs_filtered


def make_cumulative_TE_table(orthogroups_path:str, n:int, species:str, repeats_annot_path:str, genome_annot_path:str, sig_orthogroups = [], count_transcripts = False):
    """
    This function only analyzes one species at a time!
    make a table for the surrounding n bases upstream and downstream of each gene where each base is a row 
    and each column is the sum of how often this base is annotated as the TE-category across all transcripts

    optionally: provide a sig_orthogroups list to not include all transcripts that are in orthogroups_path
    
    if count_transcripts=True it only returns a list that includes all the transcripts that were included in the computtion
    """
    orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_path, sig_orthogroups, species=species)
    all_transcript_IDs = get_sig_transcripts(orthoDB_orthogroups)
    print(f"{species} --> {len(orthoDB_orthogroups)} orthogroups with {len(all_transcript_IDs)} transcripts (includes orthogroups with no members in {species}")

    print(f"\t*  parse gene-annotation {genome_annot_path}")
    genes_dict = gff.parse_gff3_general(genome_annot_path, verbose=False, keep_feature_category=gff.FeatureCategory.Transcript)
    print(f"\t*  parse repeat-annotation {repeats_annot_path}")
    repeats_dict = repeats.parse_repeats_repeatmasker_outfile(repeats_annot_path, verbose=False)

    repeats_categories = get_all_repeat_categories(repeats_dict=repeats_dict)
    repeats_categories.sort()
    print(f"\t*  all repeat categories: {repeats_categories}")

    before_transcript = { cat : [0]*n for cat in repeats_categories} ## dict with lists for each category from 0 to n, each sub-list covering the category at base n
    after_transcript = { cat : [0]*n for cat in repeats_categories} 
    """
    number of transcript where the base i positions before transcript start is covered by each TE category. 
    Order given by the list in repeats_categories
    {                 0                 n
        category1 : [ 0 , 5, 43, 0, 0, ...],
                      n                2*n
        category2 : [ 0 , 7, 37, 4, 0, ...], 
        ...
    }
    """
    missing_in_annot_transcripts = []
    contigs_with_no_repeats = [] 

    all_transcripts_list = []
    print(f"\t*  calculate transcript surroundings")
    for transcript_id in tqdm(all_transcript_IDs):
        try:
            transcript = genes_dict[transcript_id]
            if count_transcripts:
                all_transcripts_list.append(transcript_id)
                continue
        except:
            # raise RuntimeError(f"{transcript_id} can not be found in {genome_annot_path}")
            missing_in_annot_transcripts.append(transcript_id)
            print(f"\ttranscript: {transcript_id}")
            continue
        # start and end of the interval surrounding this transcript
        int_start = transcript.start - n
        int_stop = transcript.end + n
        
        # print(f"\t\t -  {transcript_id} interval from {int_start} to {int_stop}")
        
        ################
        ### start filling first half before the coding region
        
            # collect all repeats that are in the pre-transcript interval
            # this fails if there's no repeats on the transcript.contig, which can happen in very fragmented assemblies
        try:
            repeat_before_transcript = []
            for repeat in repeats_dict[transcript.contig]:
                rep_stop_in_interval = repeat.end < transcript.start and repeat.end > int_start
                rep_start_in_interval = repeat.start < transcript.start and repeat.start > int_start
                rep_longer_than_interval = repeat.start < int_start and repeat.end >int_stop

                if rep_stop_in_interval or rep_start_in_interval or rep_longer_than_interval:
                    repeat_before_transcript.append(repeat)
        except:
            contigs_with_no_repeats.append(transcript.contig)
            # print(f"{transcript.contig} has no repeats on it in {repeats_annot_path}")
            continue
        
        # iterate through every base in the interval
        for index, base in enumerate(range(int_start, transcript.start)):
            for repeat in repeat_before_transcript:
                if base >= repeat.start and base <= repeat.end:
                    # print(f"{repeat.start} < {base} [{index}] < {repeat.end}, {repeat}")
                    before_transcript[repeat.repeat_category][index] += 1
        
        ################
        ### fill out the second half after the coding region, same as above but post-transcript interval
        repeat_after_transcript = [repeat for repeat in repeats_dict[transcript.contig] if (repeat.start < int_stop and repeat.start>transcript.end) or (repeat.end > transcript.end and repeat.end<int_stop) or (repeat.start < transcript.end and repeat.end > int_stop)]
        for index, base in enumerate(range(transcript.end, int_stop)):
            for repeat in repeat_after_transcript:
                if base >= repeat.start and base <= repeat.end:
                    after_transcript[repeat.repeat_category][index] += 1
    
    if len(missing_in_annot_transcripts)>0:
        print(f"{len(missing_in_annot_transcripts)} transcripts (of {len(all_transcript_IDs)} total) not found in annotation and were skipped (in C. maculatus this might be due to the liftover?)")
    if len(contigs_with_no_repeats)>0:
        print(f"{len(contigs_with_no_repeats)} contigs with significant genes but no repeats on them")
    
    if not count_transcripts:
        return before_transcript, after_transcript
    elif count_transcripts:
        return all_transcripts_list




def plot_TE_abundance(before_filepath:str, after_filepath:str, sig_transcripts:int, all_before_filepath:str = "", all_after_filepath:str = "", all_transcripts:int = 0, filename = "cumulative_repeat_presence_around_transcripts.png", legend = True):
    """
    plot the cumulative repeat presence per base before and after a transcript (before and after infile paths)
    infiles generated from make_cumulative_TE_table and saved to text file
    """
    
    before_dict = gff.read_dict_from_file(before_filepath)
    before_dict = { key : [int(v)/sig_transcripts*100 for v in value] for key, value in before_dict.items()}
    after_dict = gff.read_dict_from_file(after_filepath)
    after_dict = { key : [int(v)/sig_transcripts*100 for v in value] for key, value in after_dict.items()}

    all_before_dict = {}
    all_after_dict = {}
    if all_before_filepath !="" and all_after_filepath !="" :
        
        all_before_dict = gff.read_dict_from_file(all_before_filepath)
        all_before_dict = { key : [int(v)/all_transcripts*100 for v in value] for key, value in all_before_dict.items()}
        all_after_dict = gff.read_dict_from_file(all_after_filepath)
        all_after_dict = { key : [int(v)/all_transcripts*100 for v in value] for key, value in all_after_dict.items()}
        
    colors = {
        'Unknown' : "#C1C1C1" , # light grey
        # orange
        'DNA' : "#FF9000" , # Princeton orange
        # green
        'LTR' : "#6E8448" , # reseda green
        'RC' : "#8EA861" , # asparagus 
        # red
        'tRNA' : "#C14953" , # bittersweet shimmer
        'rRNA' : "#D0767E" , # old rose
        'snRNA' : "#7A2A30" , # wine
        # blue 
        'LINE' : "#3476AD" , # UCLA blue
        'SINE': "#72A8D5" , # ruddy blue
        'SINE?': "#72A8D5" , # ruddy blue
        # '' : "#2A618D" , #lapis lazuli
        # dark red-brown
        'Low_complexity' : "#3A3335" , # Jet 
        'Satellite' : "#564D4F" , #Wenge 
        'Simple_repeat' : "#827376" , #Taupe gray
    }

    fs = 25 # set font size
    if legend:
        fig, ax = plt.subplots(1, 1, figsize=(20, 10))
    else:
        fig, ax = plt.subplots(1, 1, figsize=(20, 12))
    rep_classes = list(before_dict.keys())
    num_bp = len(before_dict[rep_classes[0]])
    x_before = range(-num_bp, 0)
    x_after = range(num_bp)

    max_percentage = 0
    for rep_class in rep_classes:

        max_before = max(before_dict[rep_class])
        max_after = max(after_dict[rep_class])
        if max_before>max_percentage:
            max_percentage=max_before
        if max_after>max_percentage:
            max_percentage=max_after

        ax.plot(x_before, before_dict[rep_class], label = rep_class, color = colors[rep_class])
        ax.plot(x_after, after_dict[rep_class], color = colors[rep_class])
        
        if all_before_dict !={} and all_after_dict !={} and all_transcripts!=0:
            max_before = max(all_before_dict[rep_class])
            max_after = max(all_after_dict[rep_class])
            if max_before>max_percentage:
                max_percentage=max_before
            if max_after>max_percentage:
                max_percentage=max_after

            ax.plot(x_before, all_before_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 10)))                    
            ax.plot(x_after, all_after_dict[rep_class], color = colors[rep_class], linestyle = (0, (1, 10)))                
    
    max_percentage=int(max_percentage*1.4)
    if max_percentage == 0 or max_percentage>100:
        max_percentage= 100

    plt.vlines(x= 0, ymin=0, ymax=max_percentage, colors="#000000", linestyles="dashed", label="transcript border")
    plt.xticks(range(-num_bp, num_bp+1, int(num_bp/5)), fontsize = fs)
    plt.yticks(range(0, max_percentage+1, 10), fontsize = fs)

    species = gff.split_at_second_occurrence(before_filepath.split("/")[-1])
    species = species.replace("_", ". ")

    if legend:
        ax.set_xlim([-num_bp, num_bp*1.35])
        plt.legend(loc = "upper right", fontsize = fs)
        plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream \n({num_sig_transcripts} significant transcripts of {all_transcripts} in CAFE analysis)", fontsize = fs*1.25)
    else:
        
        solid = Line2D([0], [0], color='black', linestyle='-', linewidth=2)
        dotted = Line2D([0], [0], color='black', linestyle=':', linewidth=2)
        handles = [solid, dotted]
        labels = []
        # handles.append(mpatches.Patch(fill=False, linestyle=None))
        # handles.append(mpatches.Patch(fill=False, linestyle=None))
        labels.append(f"significant transcripts ({num_sig_transcripts})")
        labels.append(f"all CAFE transcripts ({all_transcripts})")
        plt.legend(handles, labels, loc = "upper center", fontsize = fs)
        plt.title(f"{species} transcript surroundings {num_bp} bp up and downstream", fontsize = fs*1.25)
    
    plt.xlabel(f"basepairs upstream and downstream from transcript", fontsize = fs)

    plt.ylabel(f"percent of transcripts in which this base is a repeat", fontsize = fs)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x > 99 and x<1 else f'{int(x)}%'))

    plt.tight_layout()
    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)



def filepaths():

    repeats_dir = "/Users/milena/work/repeatmasking/repeat_gffs/"
    repeats_out = {
        "A_obtectus" : f"{repeats_dir}A_obtectus_assembly_genomic.fna.out",
        "A_verrucosus" : f"{repeats_dir}A_verrucosus_assembly_genomic.fna.out",
        "B_siliquastri" : f"{repeats_dir}B_siliquastri_assembly_genomic.fna.out",
        # "C_analis" : f"{repeats_dir}C_analis_assembly_genomic.fna.out",
        "C_chinensis" : f"{repeats_dir}C_chinensis_assembly_genomic.fna.out",
        "C_maculatus" : f"{repeats_dir}C_maculatus_assembly_genomic.fna.out",
        "C_septempunctata" : f"{repeats_dir}C_septempunctata_assembly_genomic.fna.out",
        "D_melanogaster" : f"{repeats_dir}D_melanogaster_assembly_genomic.fna.out",
        "D_ponderosae" : f"{repeats_dir}D_ponderosae_assembly_genomic.fna.out",
        "I_luminosus" : f"{repeats_dir}I_luminosus_assembly_genomic.fna.out",
        "P_pyralis" : f"{repeats_dir}P_pyralis_assembly_genomic.fna.out",
        "R_ferrugineus" : f"{repeats_dir}R_ferrugineus_assembly_genomic.fna.out",
        "T_castaneum" : f"{repeats_dir}T_castaneum_assembly_genomic.fna.out",
        "T_molitor" : f"{repeats_dir}T_molitor_assembly_genomic.fna.out",
        "Z_morio" : f"{repeats_dir}Z_morio_assembly_genomic.fna.out",
    }

    repeats_dir_work = "/Users/miltr339/work/repeatmasking/repeat_gffs/"
    repeats_out_work = {
        "A_obtectus" : f"{repeats_dir_work}A_obtectus_masking.ori.out",
        "A_verrucosus" : f"{repeats_dir_work}A_verrucosus_masking.ori.out",
        "B_siliquastri" : f"{repeats_dir_work}B_siliquastri_masking.ori.out",
        "C_analis" : f"{repeats_dir_work}C_analis_masking.ori.out",
        "C_chinensis" : f"{repeats_dir_work}C_chinensis_masking.ori.out",
        "C_maculatus" : f"{repeats_dir_work}C_maculatus_superscaffolded_masking.ori.out",
        "C_septempunctata" : f"{repeats_dir_work}C_septempunctata_masking.ori.out",
        "D_melanogaster" : f"{repeats_dir_work}D_melanogaster_masking.ori.out",
        "D_ponderosae" : f"{repeats_dir_work}D_ponderosae_masking.ori.out",
        "I_luminosus" : f"{repeats_dir_work}I_luminosus_masking.ori.out",
        "P_pyralis" : f"{repeats_dir_work}P_pyralis_masking.ori.out",
        "R_ferrugineus" : f"{repeats_dir_work}R_ferrugineus_masking.ori.out",
        "T_castaneum" : f"{repeats_dir_work}T_castaneum_masking.ori.out",
        "T_molitor" : f"{repeats_dir_work}T_molitor_masking.ori.out",
        "Z_morio" : f"{repeats_dir_work}Z_morio_masking.ori.out",
    }

    orthoDB_annot_dir = "/Users/milena/work/orthoDB_annotations/"
    orthoDB_annotations = {
        "A_obtectus" : f"{orthoDB_annot_dir}A_obtectus_orthoDB_filtered.gff",
        "A_verrucosus" : f"{orthoDB_annot_dir}A_verrucosus_orthoDB_filtered.gff",
        "B_siliquastri" : f"{orthoDB_annot_dir}B_siliquastri_orthoDB_filtered.gff",
        "C_analis" : f"{orthoDB_annot_dir}C_analis_orthoDB_filtered.gff",
        "C_chinensis" : f"{orthoDB_annot_dir}C_chinensis_orthoDB_filtered.gff",
        "C_maculatus" : f"{orthoDB_annot_dir}C_maculatus_orthoDB_filtered.gff",
        "C_septempunctata" : f"{orthoDB_annot_dir}C_septempunctata_s_orthoDB_filtered.gff",
        "D_melanogaster" : f"{orthoDB_annot_dir}D_melanogaster_orthoDB_filtered.gff",
        "D_ponderosae" : f"{orthoDB_annot_dir}D_ponderosae_orthoDB_filtered.gff",
        "I_luminosus" : f"{orthoDB_annot_dir}I_luminosus_orthoDB_filtered.gff",
        "P_pyralis" : f"{orthoDB_annot_dir}P_pyralis_orthoDB_filtered.gff",
        "R_ferrugineus" : f"{orthoDB_annot_dir}R_ferrugineus_orthoDB_filtered.gff",
        "T_castaneum_s" : f"{orthoDB_annot_dir}T_castaneum_s_orthoDB_filtered.gff",
        "T_molitor" : f"{orthoDB_annot_dir}T_molitor_orthoDB_filtered.gff",
        "Z_morio" : f"{orthoDB_annot_dir}Z_morio_orthoDB_filtered.gff",
    }

    orthoDB_annot_dir_work = "/Users/miltr339/work/orthoDB_annotations/"
    orthoDB_annotations_work = {
        "A_obtectus" : f"{orthoDB_annot_dir_work}A_obtectus_braker_isoform_filtered.gff",
        "A_verrucosus" : f"{orthoDB_annot_dir_work}A_verrucosus_braker_isoform_filtered.gff",
        "B_siliquastri" : f"{orthoDB_annot_dir_work}B_siliquastri_braker_isoform_filtered.gff",
        "C_analis" : f"{orthoDB_annot_dir_work}C_analis_braker_isoform_filtered.gff",
        "C_chinensis" : f"{orthoDB_annot_dir_work}C_chinensis_braker_isoform_filtered.gff",
        "C_maculatus" : f"{orthoDB_annot_dir_work}C_maculatus_superscaffolded_annotation_isoform_filtered.gff",
        "C_septempunctata" : f"{orthoDB_annot_dir_work}C_septempunctata_braker_isoform_filtered.gff",
        "D_melanogaster" : f"{orthoDB_annot_dir_work}D_melanogaster_braker_isoform_filtered.gff",
        "D_ponderosae" : f"{orthoDB_annot_dir_work}D_ponderosae_braker_isoform_filtered.gff",
        "I_luminosus" : f"{orthoDB_annot_dir_work}I_luminosus_braker_isoform_filtered.gff",
        "P_pyralis" : f"{orthoDB_annot_dir_work}P_pyralis_braker_isoform_filtered.gff",
        "R_ferrugineus" : f"{orthoDB_annot_dir_work}R_ferrugineus_braker_isoform_filtered.gff",
        "T_castaneum" : f"{orthoDB_annot_dir_work}T_castaneum_braker_isoform_filtered.gff",
        "T_molitor" : f"{orthoDB_annot_dir_work}T_molitor_braker_isoform_filtered.gff",
        "Z_morio" : f"{orthoDB_annot_dir_work}Z_morio_braker_isoform_filtered.gff",
    }

    orthogroups_native = "/Users/miltr339/work/orthofinder/native_orthogroups/N0.tsv"
    orthogroups_orthoDB = "/Users/miltr339/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # orthogroups_native = "/Users/milena/work/orthofinder/native_orthogroups/N0.tsv"
    # orthogroups_orthoDB = "/Users/milena/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    sig_native = "/Users/miltr339/Box Sync/code/CAFE/native_from_N0_Base_family_results.txt"
    sig_orthoDB = "/Users/miltr339/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"

    return repeats_out, repeats_out_work, orthoDB_annotations, orthoDB_annotations_work, orthogroups_native, orthogroups_orthoDB, sig_native, sig_orthoDB



def tables_filepaths():
    work_out_dir = "/Users/miltr339/work/PhD_code/"
    sig_before_transcript = {
        "A_obtectus" : f"{work_out_dir}A_obtectus_cumulative_repeats_before_sig_transcripts.txt",
        "A_verrucosus" : f"{work_out_dir}A_verrucosus_cumulative_repeats_before_sig_transcripts.txt",
        "B_siliquastri" : f"{work_out_dir}B_siliquastri_cumulative_repeats_before_sig_transcripts.txt",
        "C_analis" : f"{work_out_dir}C_analis_cumulative_repeats_before_sig_transcripts.txt",
        "C_chinensis" : f"{work_out_dir}C_chinensis_cumulative_repeats_before_sig_transcripts.txt",
        "C_maculatus" : f"{work_out_dir}C_maculatus_cumulative_repeats_before_sig_transcripts.txt",
        "C_septempunctata" : f"{work_out_dir}C_septempunctata_cumulative_repeats_before_sig_transcripts.txt",
        "D_melanogaster" : f"{work_out_dir}D_melanogaster_cumulative_repeats_before_sig_transcripts.txt",
        "D_ponderosae" : f"{work_out_dir}D_ponderosae_cumulative_repeats_before_sig_transcripts.txt",
        "I_luminosus" : f"{work_out_dir}I_luminosus_cumulative_repeats_before_sig_transcripts.txt",
        "P_pyralis" : f"{work_out_dir}P_pyralis_cumulative_repeats_before_sig_transcripts.txt",
        "R_ferrugineus" : f"{work_out_dir}R_ferrugineus_cumulative_repeats_before_sig_transcripts.txt",
        "T_castaneum" : f"{work_out_dir}T_castaneum_cumulative_repeats_before_sig_transcripts.txt",
        "T_molitor" : f"{work_out_dir}T_molitor_cumulative_repeats_before_sig_transcripts.txt",
        "Z_morio" : f"{work_out_dir}Z_morio_cumulative_repeats_before_sig_transcripts.txt",
    }
    sig_after_transcript = {
        "A_obtectus" : f"{work_out_dir}A_obtectus_cumulative_repeats_after_sig_transcripts.txt",
        "A_verrucosus" : f"{work_out_dir}A_verrucosus_cumulative_repeats_after_sig_transcripts.txt",
        "B_siliquastri" : f"{work_out_dir}B_siliquastri_cumulative_repeats_after_sig_transcripts.txt",
        "C_analis" : f"{work_out_dir}C_analis_cumulative_repeats_after_sig_transcripts.txt",
        "C_chinensis" : f"{work_out_dir}C_chinensis_cumulative_repeats_after_sig_transcripts.txt",
        "C_maculatus" : f"{work_out_dir}C_maculatus_cumulative_repeats_after_sig_transcripts.txt",
        "C_septempunctata" : f"{work_out_dir}C_septempunctata_cumulative_repeats_after_sig_transcripts.txt",
        "D_melanogaster" : f"{work_out_dir}D_melanogaster_cumulative_repeats_after_sig_transcripts.txt",
        "D_ponderosae" : f"{work_out_dir}D_ponderosae_cumulative_repeats_after_sig_transcripts.txt",
        "I_luminosus" : f"{work_out_dir}I_luminosus_cumulative_repeats_after_sig_transcripts.txt",
        "P_pyralis" : f"{work_out_dir}P_pyralis_cumulative_repeats_after_sig_transcripts.txt",
        "R_ferrugineus" : f"{work_out_dir}R_ferrugineus_cumulative_repeats_after_sig_transcripts.txt",
        "T_castaneum" : f"{work_out_dir}T_castaneum_cumulative_repeats_after_sig_transcripts.txt",
        "T_molitor" : f"{work_out_dir}T_molitor_cumulative_repeats_after_sig_transcripts.txt",
        "Z_morio" : f"{work_out_dir}Z_morio_cumulative_repeats_after_sig_transcripts.txt",
    }
    
    all_before_transcript = {
        "A_obtectus" : f"{work_out_dir}A_obtectus_cumulative_repeats_before_all_transcripts.txt",
        "A_verrucosus" : f"{work_out_dir}A_verrucosus_cumulative_repeats_before_all_transcripts.txt",
        "B_siliquastri" : f"{work_out_dir}B_siliquastri_cumulative_repeats_before_all_transcripts.txt",
        "C_analis" : f"{work_out_dir}C_analis_cumulative_repeats_before_all_transcripts.txt",
        "C_chinensis" : f"{work_out_dir}C_chinensis_cumulative_repeats_before_all_transcripts.txt",
        "C_maculatus" : f"{work_out_dir}C_maculatus_cumulative_repeats_before_all_transcripts.txt",
        "C_septempunctata" : f"{work_out_dir}C_septempunctata_cumulative_repeats_before_all_transcripts.txt",
        "D_melanogaster" : f"{work_out_dir}D_melanogaster_cumulative_repeats_before_all_transcripts.txt",
        "D_ponderosae" : f"{work_out_dir}D_ponderosae_cumulative_repeats_before_all_transcripts.txt",
        "I_luminosus" : f"{work_out_dir}I_luminosus_cumulative_repeats_before_all_transcripts.txt",
        "P_pyralis" : f"{work_out_dir}P_pyralis_cumulative_repeats_before_all_transcripts.txt",
        "R_ferrugineus" : f"{work_out_dir}R_ferrugineus_cumulative_repeats_before_all_transcripts.txt",
        "T_castaneum" : f"{work_out_dir}T_castaneum_cumulative_repeats_before_all_transcripts.txt",
        "T_molitor" : f"{work_out_dir}T_molitor_cumulative_repeats_before_all_transcripts.txt",
        "Z_morio" : f"{work_out_dir}Z_morio_cumulative_repeats_before_all_transcripts.txt",
    }
    all_after_transcript = {
        "A_obtectus" : f"{work_out_dir}A_obtectus_cumulative_repeats_after_all_transcripts.txt",
        "A_verrucosus" :f"{work_out_dir}A_verrucosus_cumulative_repeats_after_all_transcripts.txt",
        "B_siliquastri" :f"{work_out_dir}B_siliquastri_cumulative_repeats_after_all_transcripts.txt",
        "C_analis" :f"{work_out_dir}C_analis_cumulative_repeats_after_all_transcripts.txt",
        "C_chinensis" :f"{work_out_dir}C_chinensis_cumulative_repeats_after_all_transcripts.txt",
        "C_maculatus" :f"{work_out_dir}C_maculatus_cumulative_repeats_after_all_transcripts.txt",
        "C_septempunctata" :f"{work_out_dir}C_septempunctata_cumulative_repeats_after_all_transcripts.txt",
        "D_melanogaster" :f"{work_out_dir}D_melanogaster_cumulative_repeats_after_all_transcripts.txt",
        "D_ponderosae" :f"{work_out_dir}D_ponderosae_cumulative_repeats_after_all_transcripts.txt",
        "I_luminosus" :f"{work_out_dir}I_luminosus_cumulative_repeats_after_all_transcripts.txt",
        "P_pyralis" :f"{work_out_dir}P_pyralis_cumulative_repeats_after_all_transcripts.txt",
        "R_ferrugineus" :f"{work_out_dir}R_ferrugineus_cumulative_repeats_after_all_transcripts.txt",
        "T_castaneum" :f"{work_out_dir}T_castaneum_cumulative_repeats_after_all_transcripts.txt",
        "T_molitor" :f"{work_out_dir}T_molitor_cumulative_repeats_after_all_transcripts.txt",
        "Z_morio" :f"{work_out_dir}Z_morio_cumulative_repeats_after_all_transcripts.txt",
    }

    repeats_tables ="/Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables/"
    threshold_before_transcript = {
        'A_obtectus' : f'{repeats_tables}A_obtectus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'A_verrucosus' : f'{repeats_tables}A_verrucosus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'B_siliquastri' : f'{repeats_tables}B_siliquastri_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'C_analis' : f'{repeats_tables}C_analis_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'C_chinensis' : f'{repeats_tables}C_chinensis_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'C_maculatus' : f'{repeats_tables}C_maculatus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'C_septempunctata' : f'{repeats_tables}C_septempunctata_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'D_melanogaster' : f'{repeats_tables}D_melanogaster_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'D_ponderosae' : f'{repeats_tables}D_ponderosae_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'I_luminosus' : f'{repeats_tables}I_luminosus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'P_pyralis' : f'{repeats_tables}P_pyralis_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'R_ferrugineus' : f'{repeats_tables}R_ferrugineus_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'T_castaneum' : f'{repeats_tables}T_castaneum_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'T_molitor' : f'{repeats_tables}T_molitor_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
        'Z_morio' : f'{repeats_tables}Z_morio_cumulative_repeats_before_sig_transcripts_90th_GF_size_percentile.txt',
    }
    threshold_after_transcript = {
        'A_obtectus' : f'{repeats_tables}A_obtectus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'A_verrucosus' : f'{repeats_tables}A_verrucosus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'B_siliquastri' : f'{repeats_tables}B_siliquastri_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'C_analis' : f'{repeats_tables}C_analis_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'C_chinensis' : f'{repeats_tables}C_chinensis_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'C_maculatus' : f'{repeats_tables}C_maculatus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'C_septempunctata' : f'{repeats_tables}C_septempunctata_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'D_melanogaster' : f'{repeats_tables}D_melanogaster_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'D_ponderosae' : f'{repeats_tables}D_ponderosae_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'I_luminosus' : f'{repeats_tables}I_luminosus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'P_pyralis' : f'{repeats_tables}P_pyralis_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'R_ferrugineus' : f'{repeats_tables}R_ferrugineus_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'T_castaneum' : f'{repeats_tables}T_castaneum_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'T_molitor' : f'{repeats_tables}T_molitor_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
        'Z_morio' : f'{repeats_tables}Z_morio_cumulative_repeats_after_sig_transcripts_90th_GF_size_percentile.txt',
    }

    return sig_before_transcript, sig_after_transcript, all_before_transcript, all_after_transcript, threshold_before_transcript, threshold_after_transcript


if __name__ == "__main__":

    # read infile paths 
    repeats_out, repeats_out_work, orthoDB_annotations, orthoDB_annotations_work, orthogroups_native, orthogroups_orthoDB, sig_native, sig_orthoDB = filepaths()

    ##########################################################
    ######### make TE abundance tables for plotting ##########
    ##########################################################
    
    ## Tables are written to csv files since computing them takes a bit
    repeats_tables = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_tables/"
    ## uv run python3 for quicker runtimes

    

    ## test the function about the GF size filtering
    if False:
        all_species = list(repeats_out.keys())
        for species in all_species:
            sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
            orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_orthoDB_list, species=species)
            OG_id_list = filter_sig_OGs_by_size(orthoDB_orthogroups, species, q=90)

            print(f"{species}: \n\tunfiltered list:{sig_orthoDB_list[0:6]} ({len(sig_orthoDB_list)})\n\tfiltered list: {OG_id_list[0:6]} ({len(OG_id_list)})")

    
    ########
    ## tables for significant transcripts
    ########
    size_percentile_threshold = 90
    
    if False:
        all_species = list(repeats_out.keys())
        # failed = ['A_verrucosus', 'C_chinensis', 'D_ponderosae', 'I_luminosus', 'R_ferrugineus', 'T_molitor', 'Z_morio']
        failed = []
        all_species = ['C_maculatus']
        for species in all_species:
            print(f"orthoDB {species}: ")
            
            sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
            orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_orthoDB_list, species=species)
            sig_OGs_size_filtered = filter_sig_OGs_by_size(orthoDB_orthogroups=orthoDB_orthogroups, species=species, q=size_percentile_threshold)
            if True:
            #try:
                before_transcript, after_transcript = make_cumulative_TE_table(orthogroups_orthoDB, n=10000, species=species, repeats_annot_path=repeats_out_work[species], genome_annot_path=orthoDB_annotations_work[species], sig_orthogroups=sig_OGs_size_filtered)
                gff.write_dict_to_file(before_transcript, f"{repeats_tables}{species}_cumulative_repeats_before_sig_transcripts_{size_percentile_threshold}th_GF_size_percentile.txt")
                gff.write_dict_to_file(after_transcript, f"{repeats_tables}{species}_cumulative_repeats_after_sig_transcripts_{size_percentile_threshold}th_GF_size_percentile.txt")
            # except: 
            #     failed.append(species)
        print(f"failed species: {failed}")


    #########
    ## tables for all CAFE transcripts
    #########
    if False:
        all_species = list(repeats_out.keys())
        # failed = ['A_verrucosus', 'C_chinensis', 'D_ponderosae', 'I_luminosus', 'R_ferrugineus', 'T_molitor', 'Z_morio']
        failed = []
        for species in all_species:
            print(f"orthoDB {species}: ")
            ## uv run python3 for quicker runtimes
            sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
            try:
                # before_transcript, after_transcript = make_cumulative_TE_table(orthogroups_orthoDB, n=50, species=species, repeats_annot_path=repeats_out[species], genome_annot_path=orthoDB_annotations[species], sig_orthogroups=sig_orthoDB_list)
                before_transcript, after_transcript = make_cumulative_TE_table(orthogroups_orthoDB, n=10000, species=species, repeats_annot_path=repeats_out_work[species], genome_annot_path=orthoDB_annotations_work[species])
                gff.write_dict_to_file(before_transcript, f"{species}_cumulative_repeats_before_all_transcripts.txt")
                gff.write_dict_to_file(after_transcript, f"{species}_cumulative_repeats_after_all_transcripts.txt")
            except: 
                failed.append(species)
        print(f"failed species: {failed}")



    ######################################################
    ############ plot above computed tables ##############
    ######################################################

    sig_before_transcript, sig_after_transcript, all_before_transcript, all_after_transcript, threshold_before_transcript, threshold_after_transcript = tables_filepaths()

    if True:
        repeats_plots = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/repeats_cumulative_around_transcripts/"
        all_species = list(repeats_out.keys())
        failed_species = []
        for species in all_species:
            # species = "B_siliquastri"
            print(f"\n\n ------------------------------------------------------------- \nplot {species}")
            # get total number of transcripts that are part of significantly rapidly evolving orthogroups in this species
            sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)

            orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_orthoDB_list, species=species)
            sig_OGs_size_filtered = filter_sig_OGs_by_size(orthoDB_orthogroups=orthoDB_orthogroups, species=species, q=size_percentile_threshold)
            
            all_transcript_IDs = get_sig_transcripts(orthoDB_orthogroups)
            num_sig_transcripts_ = len(all_transcript_IDs)
            # print(f"\t{num_sig_transcripts_} significant transcripts according to reading the CAFE and orthoDB output")
            
            all_transcript_IDs = make_cumulative_TE_table(orthogroups_orthoDB, n=10000, species=species, repeats_annot_path=repeats_out_work[species], genome_annot_path=orthoDB_annotations_work[species], sig_orthogroups=sig_OGs_size_filtered, count_transcripts=True)
            num_sig_transcripts = len(all_transcript_IDs)
            if num_sig_transcripts_ != num_sig_transcripts:
                print(f"\t!!! {num_sig_transcripts} significant transcrips according to making the table\n\t!!! {num_sig_transcripts_} transcripts according to CAFE and orthoDB output\n")
                failed_species.append(species)

            # plot_TE_abundance(threshold_before_transcript[species], threshold_after_transcript[species], sig_transcripts = num_sig_transcripts, filename=f"{repeats_plots}cumulative_repeat_presence_around_transcripts_sig_only_{species}_{size_percentile_threshold}th_percentile_GF_size.png")

            orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_orthoDB, all_orthogroups_list, species=species)
            all_transcript_IDs = get_sig_transcripts(orthoDB_orthogroups)
            num_all_transcripts_ = len(all_transcript_IDs)
            # print(f"\t{num_all_transcripts_} CAFE transcripts according to reading the CAFE and orthoDB output")
            # write table to output file
            all_transcript_IDs = make_cumulative_TE_table(orthogroups_orthoDB, n=10000, species=species, repeats_annot_path=repeats_out_work[species], genome_annot_path=orthoDB_annotations_work[species], count_transcripts=True)
            num_all_transcripts = len(all_transcript_IDs)
            if num_all_transcripts_ != num_all_transcripts:
                print(f"\t!!! {num_all_transcripts} CAFE transcrips according to making the table\n\t!!! {num_all_transcripts_} transcripts according to CAFE and orthoDB output")
                failed_species.append(species)

            plot_TE_abundance(threshold_before_transcript[species], threshold_after_transcript[species], sig_transcripts = num_sig_transcripts, all_before_filepath=all_before_transcript[species], all_after_filepath=all_after_transcript[species], all_transcripts=num_all_transcripts, filename=f"{repeats_plots}cumulative_repeat_presence_around_transcripts_sig_and_all_{species}_{size_percentile_threshold}th_percentile_GF_size.png", legend=False)
            # break
        print(f"\n-----------------------------------------------------\nspecies where the transcript numbers don't match: {failed_species}")



