# use the cmac curated TE library to do a blast search on the orthoDB (unfiltered) annotation to see if the filtered-out transcripts have more TEs than expected

import matplotlib.pyplot as plt
from Bio import SeqIO
import typing as tp
from dataclasses import dataclass
import pandas as pd
# dataclass generates stuff like __init__() for classes automatically 


@dataclass
class TEMatch:
    TE_id:str
    seq_ident:float
    evalue:float
    bitscore:float
    TE_type:str

@dataclass
class Transcript:
    ID:str
    TE_matches:tp.List[TEMatch]

    def addTE(self, newTE:TEMatch):
        self.TE_matches.append(newTE)

    
    def filter_by_seqident(self, min_ident = 80):
        """ 
        Filter the TE_matches list to only contain matches that have a sequence identity above min_ident
        """
        filtered_matches = []
        for match in self.TE_matches:
            if match.seq_ident>min_ident:
                filtered_matches.append(match)
        self.TE_matches = filtered_matches
    
    @property
    def show(self):
        print(self.ID)
        match_types = list(set([match.TE_type for match in self.TE_matches]))
        print(f"There are {len(self.TE_matches)} hits of  {len(match_types)} different type(s) of TEs")
        for match_type in match_types:
            print(f" -> TE type: {match_type}")
            for match in self.TE_matches:
                if match.TE_type == match_type:
                    print(f"\t* {match.TE_id} : seq ident = {match.seq_ident}%")




@dataclass
class TranscriptMatch:
    Transcript_id:str
    seq_ident:float
    evalue:float
    bitscore:float
    TE_type:str

@dataclass
class TE_category:
    ID:str
    Transcript_matches:tp.Dict[str,tp.List[TranscriptMatch]] #{ TE_id : TranscriptMatch }
    """ dictionary of { TE_id : [TranscriptMatch1, TranscriptMatch2, ... ]}"""

    def addTranscript(self, newTranscript:TranscriptMatch):
        transcript_id = newTranscript.Transcript_id
        if transcript_id not in self.Transcript_matches:
            self.Transcript_matches[transcript_id] = []
        self.Transcript_matches[transcript_id].append(newTranscript)

    @property
    def num_transcripts(self):
        return(len(self.Transcript_matches))
    


    def remove_transcripts(self, transcript_ids_list, verbose = False):
        """ 
        filter all transcripts to only keep ones also represented in the filtered proteinsets. 
        Returns a Dictionary that can be used in TE_category.Transcript_matches 
        """
        filtered_dict = {}
        if verbose:
            len_before = len(self.Transcript_matches)
        for transcript_match, match_object in self.Transcript_matches.items():
            if transcript_match in transcript_ids_list:
                filtered_dict[transcript_match] = match_object
        if verbose:
            print(f"\t{self.ID}:\t\t {len_before} Transcripts before filtering and {len(filtered_dict)} Transcripts after filtering")
        # return {self.ID : filtered_dict}
        return(filtered_dict)
    
    @property
    def print_num_transcript_matches(self):
        # sometimes one transcripts has multiple hits of the same kind of TE in it. show all these cases
        print(self.ID)
        num_single = 0
        num_multiple = 0
        multiple_IDs_list = []
        for transcript_id, num_matches in self.Transcript_matches.items():
            if len(num_matches) == 1:
                num_single+=1
            else:
                num_multiple+=1
                multiple_IDs_list.append(transcript_id)

        print_stop_ind = 5
        if len(multiple_IDs_list)<print_stop_ind:
            print_stop_ind = len(multiple_IDs_list)
        
        print(f"\t{num_single} transcripts with one transcript match")
        if num_multiple>0:
            print(f"\t{num_multiple} transcripts with multiple transcript matches, e.g. {multiple_IDs_list[:print_stop_ind]}")
        else:
            print(f"\t{num_multiple} transcripts with multiple transcript matches")
        

def get_transctipt_list(proteinfasta_filepath):
    seq_name = []
    for record in SeqIO.parse(proteinfasta_filepath, "fasta"):
        seq_name.append(record.id)
    return seq_name

def read_blast_output(blast_filepath):
    """
    read blast output into two dictionaries: blast_outfile_transcripts, blast_outfile_TE_types
    That are made up of the Transcript and TE_Category dataclasses respectively
    """
    blast_outfile_transcripts = {}
    blast_outfile_TE_types = {}
    with open(blast_filepath, "r") as blast_file:
        for blast_line in blast_file.readlines():
            blast_line = blast_line.strip().split("\t")

            TE_string = blast_line[0]
            transcript = blast_line[1]
            seq_ident = float(blast_line[2])
            evalue = float(blast_line[10])
            bitscore = float(blast_line[11])
            TE_type=TE_string.split("#")[-1]

            TE_match = TEMatch(TE_id=TE_string, seq_ident=seq_ident, evalue=evalue, bitscore=bitscore, TE_type=TE_type)
            if transcript not in blast_outfile_transcripts:
                blast_outfile_transcripts[transcript] = Transcript(ID=transcript, TE_matches=[TE_match])
            else:
                blast_outfile_transcripts[transcript].addTE(TE_match)

            transcript_match = TranscriptMatch(Transcript_id=transcript, seq_ident=seq_ident, evalue=evalue, bitscore=bitscore, TE_type=TE_type)
            if TE_type not in blast_outfile_TE_types:
                blast_outfile_TE_types[TE_type] = TE_category(ID=TE_type, Transcript_matches={transcript : [TranscriptMatch]})
            else:
                blast_outfile_TE_types[TE_type].addTranscript(transcript_match)

    return(blast_outfile_transcripts, blast_outfile_TE_types)


def get_TE_transcripts_filtered_stats(blast_hits, unfiltered_fasta_filepath, filtered_fasta_filepath):
    # blast hits are the output from read_blast_output()

    # TODO I don't think this is right


    TE_transcripts_unfiltered = list(blast_hits.keys())
    filtered_all_transcripts = get_transctipt_list(filtered_fasta_filepath)
    unfiltered_all_transcripts = get_transctipt_list(unfiltered_fasta_filepath)
    TE_transcripts_filtered = list(set(filtered_all_transcripts) & set(TE_transcripts_unfiltered))

    unfiltered_ratio = (len(TE_transcripts_unfiltered)/len(unfiltered_all_transcripts))*100
    filtered_ratio = (len(TE_transcripts_filtered)/len(filtered_all_transcripts))*100

    print(f"\t{len(unfiltered_all_transcripts)} orthoDB transcripts, {len(TE_transcripts_unfiltered)} of them have TE sequences in them ({unfiltered_ratio:.2f}%)")
    print(f"\t{len(filtered_all_transcripts)} orthoDB transcripts, {len(TE_transcripts_filtered)} of them have TE sequences in them ({filtered_ratio:.2f}%)")


def plot_TE_histogram(TE_categories_unfiltered_dict, TE_categories_filtered_dict = {}, TE_categories_native_dict = {}, plot_filename = "", min_abundance = 100, verbose = False):

    # plot all species at once in a grid
    plot_legend = False
    if len(TE_categories_filtered_dict)>0:
        cols = len(TE_categories_unfiltered_dict)# +1 # add one empty space for the legend
        fig_len = (cols+1)*8
        # plot_legend = True
    else:
        cols = len(TE_categories_unfiltered_dict)
        fig_len = cols*8
    rows = 1
    fig, axes = plt.subplots(rows, cols, figsize=(fig_len, 8))

    if verbose:
        print(f"plot in {rows} x {cols} grid")
    colors = {"unfiltered" : "#4d7298", "filtered" : "#F2933A", "native" : "#b82946"}

    fs_main = 30
    plt.rcParams.update({'font.size': fs_main})

    # Loop over each file path and corresponding subplot axis
    for idx, species_name in enumerate(TE_categories_unfiltered_dict.keys()):
        
        # fontsize (set up for changes per loop later)
        fs = fs_main
        fs_modified = fs_main +10 # modified font size 

        TE_categories_unfiltered_species = TE_categories_unfiltered_dict[species_name]
        
        if verbose:
            print(species_name)
            print(f"{len(TE_categories_unfiltered_species.keys())} TE categories in the unfiltered transcripts")

        TE_amounts_unfiltered_dict = {}
        TE_amounts_filtered_dict = {}
        TE_amounts_native_dict = {}
        TE_ids_list = []

        for TE_id in TE_categories_unfiltered_species.keys():
            if TE_categories_unfiltered_species[TE_id].num_transcripts > min_abundance:
                TE_amounts_unfiltered_dict[TE_id] = TE_categories_unfiltered_species[TE_id].num_transcripts
                TE_ids_list.append(TE_id)

        if len(TE_categories_filtered_dict)>0:
            TE_categories_filtered_species = TE_categories_filtered_dict[species_name]
            if verbose:
                print(f"{len(TE_categories_filtered_species.keys())} TE categories in the filtered transcripts")

        if len(TE_categories_native_dict)>0:
            TE_categories_native_species = TE_categories_native_dict[species_name]
            if verbose:
                print(f"{len(TE_categories_native_species.keys())} TE categories in the native transcripts")

            # get TE abundance for all the categories
            # still check for IDs and abundance with unfiltered data so that the columns are the same
            for TE_id in TE_categories_unfiltered_species.keys():
                if TE_categories_unfiltered_species[TE_id].num_transcripts > min_abundance:
                    try:
                        TE_amounts_filtered_dict[TE_id] = TE_categories_filtered_species[TE_id].num_transcripts
                    except:
                        TE_amounts_filtered_dict[TE_id] = 0 
                    try:
                        TE_amounts_native_dict[TE_id] = TE_categories_native_species[TE_id].num_transcripts
                    except:
                        TE_amounts_native_dict[TE_id] = 0 

            if verbose:
                print(f"native TEs: {TE_amounts_native_dict}")
                print(f"filtered TEs: {TE_amounts_filtered_dict}")
        if verbose:
            print(f"unfiltered TEs: {TE_amounts_unfiltered_dict}")
            print()

        # Calculate row and column indices for the current subplot
        row = idx // cols
        col = idx % cols
        if verbose:
            print(f"\tin position {row+1},{col+1}: \t{species_name}")

        # make barplots

        # TE_amounts_unfiltered_dict = {TE_id : num for TE_id, num in TE_amounts_unfiltered_dict.items() if num>min_abundance}
        # TE_ids_list = list(TE_amounts_unfiltered_dict.keys())
        axes[col].bar(TE_ids_list, [TE_amounts_unfiltered_dict[TE_id] for TE_id in TE_ids_list], width = 1, color = colors["unfiltered"], label = "unfiltered TEs")

        if len(TE_categories_filtered_dict)>0:
            # axes[col].bar(TE_ids_list, [TE_amounts_filtered_dict[TE_id]   for TE_id in TE_ids_list if TE_amounts_unfiltered_dict[TE_id]-TE_amounts_filtered_dict[TE_id] != 0 and TE_amounts_unfiltered_dict[TE_id] > min_abundance], width = 1, color = colors["filtered"],   label = "filtered TEs")
            axes[col].bar(TE_ids_list, [TE_amounts_filtered_dict[TE_id]   for TE_id in TE_ids_list], width = 1, color = colors["filtered"],   label = "filtered TEs")
        
        if len(TE_categories_native_dict)>0:
            # axes[col].bar(TE_ids_list, [TE_amounts_filtered_dict[TE_id]   for TE_id in TE_ids_list if TE_amounts_unfiltered_dict[TE_id]-TE_amounts_filtered_dict[TE_id] != 0 and TE_amounts_unfiltered_dict[TE_id] > min_abundance], width = 1, color = colors["filtered"],   label = "filtered TEs")
            axes[col].bar(TE_ids_list, [TE_amounts_native_dict[TE_id]   for TE_id in TE_ids_list], width = 1, color = colors["native"],   label = "native TEs", alpha = 0.5)


        axes[col].set_title(f'{species_name}')
        axes[col].set_xlabel('')
        axes[col].set_ylabel('')
        # reduce font size if there's too many axis tick labels
        if len(TE_ids_list) > 30:
            fsm = fs_modified-29
            if verbose:
                print(f"long reduced fontsize in {species_name}")
        elif len(TE_ids_list) > 29:
            fsm = fs_modified-26
        elif len(TE_ids_list) > 25:
            fsm = fs_modified-25
        elif len(TE_ids_list) > 10:
            fsm = fs_modified-20
        axes[col].tick_params(axis='both', which='major', labelsize=fsm)
        axes[col].tick_params("x", rotation=90)
        axes[col].tick_params("y", labelsize = fs)
        
        if idx == 0:
            axes[col].set_ylabel("number of transcripts", fontsize = fs)


    # add legend to the final plot: 
    if len(TE_categories_filtered_dict)>0 and plot_legend:
        last_row = (idx + 1) // cols
        last_col = (idx + 1) % cols
        handles = [
            plt.Line2D([0], [0], color=colors["unfiltered"], lw=6, label="unfiltered TEs"),
            plt.Line2D([0], [0], color=colors["filtered"], lw=6, label="filtered TEs"),
        ]

        axes[last_col].legend(handles=handles, loc='center left', fontsize=fs_main, markerscale=1.5)
        axes[last_col].axis('off')  # Turn off axes for the legend-only plot

        for i in range(idx+1,rows*cols):
            last_row = i // cols
            last_col = i % cols
            axes[last_col].axis('off')  # Turn off axes for all empty plot spaces

    # Set a single x-axis label for all subplots
    x_label = f"only TEs that are present in more than {min_abundance} transcripts included"
    # fig.text(0.5, 0.04, x_label, ha='center', va='center', fontsize=fs_main)
    # fig.text(0.5, 0.96, x_label, ha='center', va='center', fontsize=fs_main)

    
    # Adjust layout to prevent overlap
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    # Show the plot
    # plt.show()
    if len(plot_filename)== 0:
        plot_filename = "TE_abundance_bruchids_histogram.png"
    plt.savefig(plot_filename, dpi = 300, transparent = True)
    print(f"plot saved in current working directory as: {plot_filename}")



def filter_TEs_by_fasta(TE_category_by_species, filtered_proteinfasta_filepath, verbose = False):

        TE_categories_filtered_dict = {}
        species_transcripts_list = get_transctipt_list(filtered_proteinfasta_filepath)
        for TE_category_str in TE_category_by_species.keys():
            TE_categories_filtered_dict[TE_category_str] = TE_category(ID = TE_category_str, Transcript_matches=TE_category_by_species[TE_category_str].remove_transcripts(species_transcripts_list, verbose = verbose))

        return(TE_categories_filtered_dict) 
    

def filter_fasta_by_TE_matches(blast_outfile, proteinfasta_file, min_ident = 90, verbose=False):
    """
    Filter a proteinfasta file to remove all entries that have matches with TE families with a sequence identity above min_ident
    """
    transcripts, TE_categories = read_blast_output(blast_outfile)

    transcripts_filtered = []
    for transcript_id in transcripts.keys():
        transcripts[transcript_id].filter_by_seqident(min_ident = min_ident)
        if len(transcripts[transcript_id].TE_matches) >0:
            transcripts_filtered.append(transcript_id)
    if verbose:
        print(f"there are {len(transcripts_filtered)} transcripts with TE-hits after filtering for sequence identity above {min_ident}%\n")   

    output_filename = proteinfasta_file.split(".fa")[0] + "_TE_filtered.fa"
    with open(output_filename, "w") as output:
        for record in SeqIO.parse(proteinfasta_file, "fasta"):
            if record.id not in transcripts_filtered:
                SeqIO.write(record, output, "fasta")
    
    if verbose:
        print(f"filtered output file saved in {output_filename}")


def get_filtering_stats(blast_outfile, proteinfasta_file, min_ident = 90, verbose = False):
    """
    Get all the numbers around the filtering for TE-containing proteins
    """
    transcripts, TE_categories = read_blast_output(blast_outfile)


    transcripts_filtered = []
    for transcript_id in transcripts.keys():
        transcripts[transcript_id].filter_by_seqident(min_ident = min_ident)
        if len(transcripts[transcript_id].TE_matches) >0:
            transcripts_filtered.append(transcript_id)
    if verbose:
        print(f"{len(transcripts_filtered)} transcripts with TE-hits above {min_ident}%")   

    output_filename = proteinfasta_file.split(".fa")[0] + "_TE_filtered.fa"
    
    all_prot = 0
    filtered_prot = 0
    with open(output_filename, "w") as output:
        for record in SeqIO.parse(proteinfasta_file, "fasta"):
            all_prot+=1
            if record.id not in transcripts_filtered:
                filtered_prot+=1
    if verbose:
        print(f"proteinfasta is downfiltered from {all_prot} to {filtered_prot} sequences")

    return([all_prot, filtered_prot, len(transcripts_filtered)])


if __name__ == "__main__":

    ## Bruchids TEs
    if False: 
        # blast search from /Users/milena/Box Sync/code/annotation_pipeline/filtering_annotation/blast_repeat_library.sh

        proteinfasta_dir = "/Users/milena/work/orthoDB_proteinseqs/"
        unfiltered_orthoDB_proteinfasta = {
            "A_obtectus" : f"{proteinfasta_dir}A_obtectus_filtered_proteinfasta.fa",
            "A_verrucosus" : f"{proteinfasta_dir}A_verrucosus_filtered_proteinfasta.fa",
            "B_siliquastri" : f"{proteinfasta_dir}B_siliquastri_filtered_proteinfasta.fa",
            "C_chinensis" : f"{proteinfasta_dir}C_chinensis_filtered_proteinfasta.fa",
            "C_maculatus" : f"{proteinfasta_dir}C_maculatus_filtered_proteinfasta.fa",
            "C_septempunctata" : f"{proteinfasta_dir}C_septempunctata_filtered_proteinfasta.fa",
            "D_melanogaster" : f"{proteinfasta_dir}D_melanogaster_filtered_proteinfasta.fa",
            "D_ponderosae" : f"{proteinfasta_dir}D_ponderosae_filtered_proteinfasta.fa",
            "I_luminosus" : f"{proteinfasta_dir}I_luminosus_filtered_proteinfasta.fa",
            "P_pyralis" : f"{proteinfasta_dir}P_pyralis_filtered_proteinfasta.fa",
            "R_ferrugineus" : f"{proteinfasta_dir}R_ferrugineus_filtered_proteinfasta.fa",
            "T_castaneum" : f"{proteinfasta_dir}T_castaneum_filtered_proteinfasta.fa",
            "T_molitor" : f"{proteinfasta_dir}T_molitor_filtered_proteinfasta.fa",
            "Z_morio" : f"{proteinfasta_dir}Z_morio_filtered_proteinfasta.fa"
        }

        proteinfasta_filt_dir = "/Users/milena/work/orthoDB_proteinseqs/overlap_filtered_proteinseqs/"
        filtered_orthoDB_proteinfasta = {
            "A_obtectus" : f"{proteinfasta_filt_dir}A_obtectus_filtered_proteinfasta_overlap_filtered.fa",
            "A_verrucosus" : f"{proteinfasta_filt_dir}A_verrucosus_filtered_proteinfasta_overlap_filtered.fa",
            "B_siliquastri" : f"{proteinfasta_filt_dir}B_siliquastri_filtered_proteinfasta_overlap_filtered.fa",
            "C_chinensis" : f"{proteinfasta_filt_dir}C_chinensis_filtered_proteinfasta_overlap_filtered.fa",
            "C_maculatus" : f"{proteinfasta_filt_dir}C_maculatus_filtered_proteinfasta_overlap_filtered.fa",
            "C_septempunctata" : f"{proteinfasta_filt_dir}C_septempunctata_filtered_proteinfasta_overlap_filtered.fa",
            "D_melanogaster" : f"{proteinfasta_filt_dir}D_melanogaster_filtered_proteinfasta_overlap_filtered.fa",
            "D_ponderosae" : f"{proteinfasta_filt_dir}D_ponderosae_filtered_proteinfasta_overlap_filtered.fa",
            "I_luminosus" : f"{proteinfasta_filt_dir}I_luminosus_filtered_proteinfasta_overlap_filtered.fa",
            "P_pyralis" : f"{proteinfasta_filt_dir}P_pyralis_filtered_proteinfasta_overlap_filtered.fa",
            "R_ferrugineus" : f"{proteinfasta_filt_dir}R_ferrugineus_filtered_proteinfasta_overlap_filtered.fa",
            "T_castaneum" : f"{proteinfasta_filt_dir}T_castaneum_filtered_proteinfasta_overlap_filtered.fa",
            "T_molitor" : f"{proteinfasta_filt_dir}T_molitor_filtered_proteinfasta_overlap_filtered.fa",
            "Z_morio" : f"{proteinfasta_filt_dir}Z_morio_filtered_proteinfast_overlap_filtereda.fa"
        }

        TE_dir = "/Users/milena/work/bruchids_TE_analysis/"
        blast_results_e30 = {
            "C_chinensis" : f"{TE_dir}chinensis_blast_cmac_repeats_evalue30.out",
            "C_maculatus" : f"{TE_dir}maculatus_blast_cmac_repeats_evalue30.out",
            "A_obtectus" : f"{TE_dir}obtectus_blast_cmac_repeats_evalue30.out",
            "B_siliquastri" : f"{TE_dir}siliquastri_blast_cmac_repeats_evalue30.out",
        }

        blast_results_native = {
            "C_chinensis" : f"{TE_dir}native_chinensis_blast_cmac_repeats_evalue30.out",
            "C_maculatus" : f"{TE_dir}native_maculatus_blast_cmac_repeats_evalue30.out",
            "A_obtectus" : f"{TE_dir}native_obtectus_blast_cmac_repeats.out",
            "B_siliquastri" : f"{TE_dir}native_siliquastri_blast_cmac_repeats.out"
        }


        # aobt_blast_out = read_blast_output(blast_results_e30["A_obtectus"])
        # cmac_blast_out = read_blast_output(blast_results_e30["C_maculatus"])
        # cchi_blast_out = read_blast_output(blast_results_e30["C_chinensis"])
        # bsil_blast_out = read_blast_output(blast_results_e30["B_siliquastri"])
        
        TE_categories_by_species_unfiltered = {}
        TE_categories_by_species_filtered = {}
        TE_categories_by_species_native = {}

        for species, blast_filepath in blast_results_e30.items():
            print(f"\n===> {species}")

            blast_dict, TE_dict_native = read_blast_output(blast_results_native[species])
            TE_categories_by_species_native[species] = TE_dict_native
            print(f"{len(TE_dict_native)} different types of TEs in native dict")

            blast_dict, TE_dict = read_blast_output(blast_results_e30[species])
            TE_dict_unfiltered = TE_dict # to make sure that there's an unfiltered version that does not get overwritten by filter_TEs_by_fasta()
            TE_categories_by_species_unfiltered[species] = TE_dict_unfiltered
            print(f"{len(TE_dict)} different types of TEs in unfiltered dict")

            TE_categories_by_species_filtered[species] = filter_TEs_by_fasta(TE_dict, filtered_orthoDB_proteinfasta[species], verbose = False)
            # filter_TEs_by_fasta(TE_dict, filtered_orthoDB_proteinfasta[species], verbose = False)
            # TE_categories_by_species_filtered[species] = TE_dict
            #print(TE_categories_by_species_filtered[species]["LTR"])
            non_empty_TEs = [TE_id# TE_id 
                            for TE_id 
                            in TE_categories_by_species_filtered[species].keys() 
                            if TE_categories_by_species_filtered[species][TE_id].num_transcripts > 0
                            ]
            print(f"{len(non_empty_TEs)} different types of TEs in filtered dict")

            # get_TE_transcripts_filtered_stats(blast_dict, unfiltered_fasta_filepath=unfiltered_orthoDB_proteinfasta[species], filtered_fasta_filepath=filtered_orthoDB_proteinfasta[species])

            # for TE_value in TE_dict.values():
            #     TE_value.print_num_transcript_matches
        
        # print(f"plot unfiltered")
        # plot_TE_histogram(TE_categories_by_species_unfiltered, plot_filename = "TE_abundance_orthoDB_unfiltered_bruchids_histogram.png",  min_abundance = 50)
        # print(f"\nplot filtered")
        # plot_TE_histogram(TE_categories_by_species_unfiltered, TE_categories_by_species_filtered, plot_filename = "TE_abundance_orthoDB_filtered_bruchids_histogram.png", min_abundance=50)
        print(f"\nplot all")
        plot_TE_histogram(TE_categories_by_species_unfiltered, TE_categories_by_species_filtered, TE_categories_by_species_native, plot_filename = "TE_abundance_orthoDB_all_annot_bruchids_histogram_small_font.png", min_abundance=50)


    ## ALL TEs
    if True:
        # rsync -azP "milenatr@rackham.uppmax.uu.se:/proj/naiss2023-6-65/Milena/annotation_pipeline/repeat_libraries_blast/*" .
        all_libraries_dir = "/Users/milena/work/repeat_libraries_blast/"

        blast_outfiles = {
            "A_obtectus" : f"{all_libraries_dir}A_obtectus_blast_repeat_family.out",
            "A_verrucosus" : f"{all_libraries_dir}A_verrucosus_blast_repeat_family.out",
            "B_siliquastri" : f"{all_libraries_dir}B_siliquastri_blast_repeat_family.out",
            "C_analis" : f"{all_libraries_dir}C_analis_blast_repeat_family.out",
            "C_chinensis" : f"{all_libraries_dir}C_chinensis_blast_repeat_family.out",
            "C_maculatus" : f"{all_libraries_dir}C_maculatus_blast_repeat_family.out",
            "C_septempunctata" : f"{all_libraries_dir}C_septempunctata_blast_repeat_family.out",
            "D_melanogaster" : f"{all_libraries_dir}D_melanogaster_blast_repeat_family.out",
            "D_ponderosae" : f"{all_libraries_dir}D_ponderosae_blast_repeat_family.out",
            "I_luminosus" : f"{all_libraries_dir}I_luminosus_blast_repeat_family.out",
            "P_pyralis" : f"{all_libraries_dir}P_pyralis_blast_repeat_family.out",
            "R_ferrugineus" : f"{all_libraries_dir}R_ferrugineus_blast_repeat_family.out",
            "T_castaneum" : f"{all_libraries_dir}T_castaneum_blast_repeat_family.out",
            "T_molitor" : f"{all_libraries_dir}T_molitor_blast_repeat_family.out",
            "Z_morio" : f"{all_libraries_dir}Z_morio_blast_repeat_family.out"
        }

        proteinfasta_dir = "/Users/milena/work/orthoDB_proteinseqs/"
        proteinfasta_files = {
            "A_obtectus" : f"{proteinfasta_dir}A_obtectus_filtered_proteinfasta.fa",
            "A_verrucosus" : f"{proteinfasta_dir}A_verrucosus_filtered_proteinfasta.fa",
            "B_siliquastri" : f"{proteinfasta_dir}B_siliquastri_filtered_proteinfasta.fa",
            "C_analis" : f"{proteinfasta_dir}C_analis_filtered_proteinfasta.fa",
            "C_chinensis" : f"{proteinfasta_dir}C_chinensis_filtered_proteinfasta.fa",
            "C_maculatus" : f"{proteinfasta_dir}C_maculatus_filtered_proteinfasta.fa",
            "C_septempunctata" : f"{proteinfasta_dir}C_septempunctata_filtered_proteinfasta.fa",
            "D_melanogaster" : f"{proteinfasta_dir}D_melanogaster_filtered_proteinfasta.fa",
            "D_ponderosae" : f"{proteinfasta_dir}D_ponderosae_filtered_proteinfasta.fa",
            "I_luminosus" : f"{proteinfasta_dir}I_luminosus_filtered_proteinfasta.fa",
            "P_pyralis" : f"{proteinfasta_dir}P_pyralis_filtered_proteinfasta.fa",
            "R_ferrugineus" : f"{proteinfasta_dir}R_ferrugineus_filtered_proteinfasta.fa",
            "T_castaneum" : f"{proteinfasta_dir}T_castaneum_filtered_proteinfasta.fa",
            "T_molitor" : f"{proteinfasta_dir}T_molitor_filtered_proteinfasta.fa",
            "Z_morio" : f"{proteinfasta_dir}Z_morio_filtered_proteinfasta.fa",
        }

        filter_stats = {}
        for species_name in blast_outfiles.keys():
            filter_fasta_by_TE_matches(blast_outfiles[species_name], proteinfasta_files[species_name], verbose=True)
            filter_stats[species_name] = get_filtering_stats(blast_outfiles[species_name], proteinfasta_files[species_name])
        
        filter_df = pd.DataFrame.from_dict(filter_stats, orient="index", columns=["before_filtering", "after_filtering", "num_filtered"])
        print(filter_df)
            
