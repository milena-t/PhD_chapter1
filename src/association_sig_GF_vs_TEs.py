"""
Associate the surroundings of genes that are part of rapidly evolving gene families
for TE presence via a stacked barplot. Each bar is a base, and the coverage is how many TEs of each class are
annotated on that base
"""

import parse_gff as gff
import parse_repeats as repeats
import parse_orthogroups as OGs

from tqdm import tqdm

def get_all_repeat_categories(repeats_dict):
    all_contigs = []
    for contig in repeats_dict.keys():
        all_contigs.append(list(set([repeat_instance.repeat_category for repeat_instance in repeats_dict[contig]])))
    all_types = [rep_type for contig_list in all_contigs for rep_type in contig_list]
    all_types = list(set(all_types))
    return(all_types)




def make_cumulative_TE_table(orthogroups_path:str, n:int, species:str, repeats_annot_path:str, genome_annot_path:str, sig_orthogroups = []):
    """
    This function only analyzes one species at a time!
    make a table for the surrounding n bases upstream and downstream of each gene where each base is a row 
    and each column is the sum of how often this base is annotated as the TE-category across all transcripts
    optionally: provide a sig_transcripts list to not include all transcripts that are in orthogroups_path
    """
    orthoDB_orthogroups = OGs.parse_orthogroups_dict(orthogroups_path, sig_orthogroups, species=species)
    genes_dict = gff.parse_gff3_general(genome_annot_path, verbose=False, keep_feature_category=gff.FeatureCategory.Transcript)
    transcripts_ids = list(genes_dict.keys())
    repeats_dict = repeats.parse_repeats_repeatmasker_outfile(repeats_annot_path, verbose=False)
    repeats_categories = get_all_repeat_categories(repeats_dict=repeats_dict)
    repeats_categories.sort()
    print(repeats_categories)

    before_transcript = { cat : [0]*n for cat in repeats_categories} ## dict with lists for each category from 0 to n, each sub-list covering all the categories at base n
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
    all_transcript_IDs = []
    for transcripts_list in orthoDB_orthogroups.values():
        for transcript_id in transcripts_list:
            transcript_id = transcript_id[:-2] # remove the "_1" suffix
            all_transcript_IDs.append(transcript_id)
    
    for transcript_id in tqdm(all_transcript_IDs):
        try:
            transcript = genes_dict[transcript_id]
        except:
            raise RuntimeError(f"{transcript_id} can not be found in {genome_annot_path}")
        # start and end of the interval surrounding this transcript
        int_start = transcript.start - n
        int_stop = transcript.end + n

        # start filling first half before the coding region
        repeat_before_transcript = (repeat for repeat in repeats_dict[transcript.contig] if repeat.start < transcript.start or repeat.stop >int_start or (repeat.start < int_start and repeat.stop >int_stop))
        for index, base in enumerate(range(int_start, transcript.start)):
            for repeat in repeat_before_transcript:
                if base >= repeat.start and base <= repeat.stop:
                    before_transcript[repeat.repeat_category][index] += 1

        # fill out the second half after the coding region
        repeat_after_transcript = (repeat for repeat in repeats_dict[transcript.contig] if repeat.start < int_stop or repeat.stop > transcript.end or (repeat.start < transcript.end and repeat.stop > int_stop))
        for index, base in enumerate(range(transcript.end, int_stop)):
            for repeat in repeat_after_transcript:
                if base >= repeat.start and base <= repeat.stop:
                    after_transcript[repeat.repeat_category][index] += 1
            
    return before_transcript, after_transcript



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
        # '' : "#2A618D" , #lapis lazuli
        # dark red-brown
        'Low_complexity' : "#3A3335" , # Jet 
        'Satellite' : "#564D4F" , #Wenge 
        'Simple_repeat' : "#827376" , #Taupe gray
    }


if __name__ == "__main__":
    
    repeats_dir = "/Users/milena/work/repeatmasking/repeat_gffs/"
    repeats_out = {
        "A_obtectus" : f"{repeats_dir}A_obtectus_assembly_genomic.fna.out",
        "A_verrucosus" : f"{repeats_dir}A_verrucosus_assembly_genomic.fna.out",
        "B_siliquastri" : f"{repeats_dir}B_siliquastri_assembly_genomic.fna.out",
        "C_analis" : f"{repeats_dir}C_analis_assembly_genomic.fna.out",
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

    orthogroups_native = "/Users/milena/work/orthofinder/native_orthogroups/N0.tsv"
    orthogroups_orthoDB = "/Users/milena/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    sig_native = "/Users/milena/Box Sync/code/CAFE/native_from_N0_Base_family_results.txt"
    sig_orthoDB = "/Users/milena/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"


    species = "A_obtectus"
    print(f"orthoDB {species}: ")
    for species in repeats_out.keys():
        ## uv run python3 for quicker runtimes
        sig_orthoDB_list, all_orthogroups_list = OGs.get_sig_orthogroups(sig_orthoDB)
        before_transcript, after_transcript = make_cumulative_TE_table(orthogroups_orthoDB, n=50, species=species, repeats_annot_path=repeats_out[species], genome_annot_path=orthoDB_annotations[species], sig_orthogroups=sig_orthoDB_list)
        gff.write_dict_to_file(before_transcript, f"{species}_cumulative_repeats_before_sig_transcripts.txt")
        gff.write_dict_to_file(after_transcript, f"{species}_cumulative_repeats_after_sig_transcripts.txt")
    
