"""
Look at all the transcripts in each orthogroup and their position
Are they on the same or different contig? how close/distant are they?
"""
import parse_gff as gff
import parse_orthogroups as OGs
import matplotlib.pyplot as plt
from Bio import SeqIO


def filepaths_native():
    native_annot_dir = "/Users/miltr339/work/native_annotations/all_native_annot/"
    native_annotations = {
        "A_obtectus" : f"{native_annot_dir}A_obtectus_annotation_isoform_filtered.gff",
        "A_verrucosus" : f"{native_annot_dir}A_verrucosus_annotation_isoform_filtered.gff",
        "B_siliquastri" : f"{native_annot_dir}B_siliquastri_annotation_isoform_filtered.gff",
        "C_analis" : f"{native_annot_dir}C_analis_annotation_isoform_filtered.gff",
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

    orthogroups_native = "/Users/miltr339/work/orthofinder/native_orthogroups/N0.tsv"
    # orthogroups_native = "/Users/milena/work/orthofinder/native_orthogroups/N0.tsv"
    sig_native = "/Users/miltr339/Box Sync/code/CAFE/native_from_N0_Base_family_results.txt"

    return native_annotations,orthogroups_native,sig_native,native_proteinseqs

def filepaths_orthoDB():
    orthoDB_annot_dir = "/Users/miltr339/work/orthoDB_annotations/"
    orthoDB_annotations = {
        "A_obtectus" : f"{orthoDB_annot_dir}A_obtectus_braker_isoform_filtered.gff",
        "A_verrucosus" : f"{orthoDB_annot_dir}A_verrucosus_braker_isoform_filtered.gff",
        "B_siliquastri" : f"{orthoDB_annot_dir}B_siliquastri_braker_isoform_filtered.gff",
        "C_analis" : f"{orthoDB_annot_dir}C_analis_braker_isoform_filtered.gff",
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
        "C_analis" : f"{orthoDB_proteinseqs_dir}C_analis_filtered_proteinfasta_TE_filtered.fa",
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

    orthogroups_orthoDB = "/Users/miltr339/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # orthogroups_orthoDB = "/Users/milena/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    
    sig_orthoDB = "/Users/miltr339/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"

    return orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs


def get_flybase_IDs(orthogroup_dict, drosophila_gff_path, outfile_name:str = "native_sig_OGs_flybase_IDs.tsv", OGs_list = []):
    """
    get all the flybase IDs from the native drosophila annotation.
    This only works with the native annotations! the orthoDB annotations don't have functional annotations to extract flybase ID from
    If you only want to include a specific set of orthogroups then include the list in OGs_list (maybe only the largest or something)
    """

    drosophila_attributes_dict = gff.parse_gff3_for_attributes(drosophila_gff_path)
    not_found_sig_IDs = []
    # print(f"{len(drosophila_attributes_dict)} transcripts in drosophila")

    with open(outfile_name, "w") as outfile:
        outfile.write("Orthogroup_ID\ttranscript_ID_native\tFlybase\n")

        for OG_id, transcripts_list in orthogroup_dict.items():
            # for weird parsing stuff i did like a year ago the transcript IDs in the native drosophila annotation have leading "__" that should be removed
            # also remove the tailing "_1" 
            transcripts_list = [transcript.replace("__", "")[:-2] for transcript in transcripts_list]
            
            if OGs_list != [] and OG_id not in OGs_list:
                continue
            
            for transcript in transcripts_list:
                try:
                    attributes = drosophila_attributes_dict[transcript]
                except:
                    # raise RuntimeError(f"the transcript ID {transcript} listed in the orthofinder output does not appear in the annotation")
                    attributes = {}
                    not_found_sig_IDs.append(transcript)
                
                try:
                    flybase = attributes["Dbxref"]
                except:
                    flybase = "None"

                outfile_string = f"{OG_id}\t{transcript}\t{flybase}\n"
                outfile.write(outfile_string)
    
    print(f"flybase IDs written to {outfile_name} in current working directory.")
    print(f"{len(not_found_sig_IDs)} transcript IDs from the significant orthogroups not found in the annotation")
    return not_found_sig_IDs
        

def make_proteinfasta_from_orthogroup(orthogroups_dict, proteinfasta_reference, orthogroups_to_include:list=[], outfile_name = "Dmel_transcripts_from_sig_OGs.fasta", species = "D_melanogaster"):
    """
    make a proteinfasta for transcripts from significant orthogroups. 
    orthogroups dict assumes a dict already resolved by species (probably drosophila)
    """

    proteinfasta = {record.id : record for record in SeqIO.parse(proteinfasta_reference,"fasta")}
    filtered_fasta = []

    for OG_id, transcripts in orthogroups_dict.items():
        if orthogroups_to_include != [] and OG_id not in orthogroups_to_include:
            continue
        for transcript in transcripts:
            transcript = f"{species}_{transcript}"
            try: 
                record = proteinfasta[transcript]
            except:
                raise RuntimeError(f"{transcript} not found in reference proteinfasta file {proteinfasta_reference}")
            
            record.id = f"{record.id}_{OG_id}"
            filtered_fasta.append(record)

    SeqIO.write(filtered_fasta, outfile_name, "fasta")
    print(f"proteinfasta written to {outfile_name}")
    return filtered_fasta




if __name__ == "__main__":

    orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs = filepaths_orthoDB()
    native_annotations,orthogroups_native,sig_native,native_proteinseqs = filepaths_native()
    
    ## Get stuff from native with functional annotations
    if False:

        native_sig_list, native_all_list =OGs.get_sig_orthogroups(sig_native)
        print(f"{len(native_sig_list)} significant orthogroups : {native_sig_list[0:10]}...")
        native_sig_OGs_dict = OGs.parse_orthogroups_dict(orthogroups_native, sig_list = native_sig_list, species="D_melanogaster")
        print(f"{len(native_sig_OGs_dict)} orthogroups in sig dict")

        native_sig_all_species = OGs.parse_orthogroups_dict(orthogroups_native, sig_list = native_sig_list)
        native_large_OGs = OGs.get_orthogroup_sizes(native_sig_all_species, q=95)
        large_OG_IDs = list(native_large_OGs.keys())

        get_flybase_IDs(native_sig_OGs_dict, native_annotations["D_melanogaster"], OGs_list=large_OG_IDs)

    # get proteinfasta for orthoDB annotations
    if True:
        orthoDB_sig_list, orthoDB_all_list =OGs.get_sig_orthogroups(sig_orthoDB)
        print(f"{len(orthoDB_sig_list)} significant orthogroups : {orthoDB_sig_list[0:10]}...")
        orthoDB_sig_OGs_dict = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_list = orthoDB_sig_list, species="D_melanogaster")
        print(f"{len(orthoDB_sig_OGs_dict)} orthogroups in sig dict")

        orthoDB_sig_all_species = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_list = orthoDB_sig_list)
        orthoDB_large_OGs = OGs.get_orthogroup_sizes(orthoDB_sig_all_species, q=95)
        large_OG_IDs = list(orthoDB_large_OGs.keys())
        
        make_proteinfasta_from_orthogroup(orthoDB_sig_OGs_dict, orthoDB_proteinseqs["D_melanogaster"], orthogroups_to_include=large_OG_IDs)