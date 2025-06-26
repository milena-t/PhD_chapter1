"""
Look at all the transcripts in each orthogroup and their position
Are they on the same or different contig? how close/distant are they?
"""
import parse_gff as gff
import parse_orthogroups as OGs
from Bio import SeqIO
import requests as req
from tqdm import tqdm
import plot_basics as my_plotting


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

    # orthogroups_orthoDB = "/Users/miltr339/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # orthogroups_orthoDB = "/Users/milena/work/orthofinder/orthoDB_uniform_masking_orthogroups/N0.tsv"
    # sig_orthoDB = "/Users/miltr339/Box Sync/code/CAFE/orthoDB_TE_filtered_Base_family_results.txt"
    
    orthogroups_orthoDB = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/N0.tsv"
    sig_orthoDB = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/CAFE_uniform_Base_Family_results.txt"

    return orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs


def make_table_with_flybase_functions(orthogroup_dict_species, drosophila_gff_path, outfile_name:str = "native_sig_OGs_flybase_IDs.tsv", OGs_list = [], orthogroups_dict_all = {}, CAFE_results_path = "", get_gene_functions_from_API = True, david_gene_groups = {}, david_functions = {}, species_tree = ""):
    """
    get all the flybase IDs from the orthogroups_dict from the native drosophila annotation. 
    orthogroups_dict_species = OGs.parse_orthogroups_dict(..., species = species) so that there is no nested dict, but it's only transcript IDs from one species in a list
    orthogroups_dict_all = OGs.parse_orthogroups_dict(..., species = "") so that there is a nested dict with all species and their separated transcripts. is only used for OG and GF sizes
    CAFE_results_path:str for the CAFE results includes the cafe p-value

    also gives you the option to access the flybase API to get functional information based on the gene ID
        --> this takes a long time! the function runs very fast otherwise if you don't use the API

    if you did a DAVID analysis of the flybase IDs, you can run this function again with a dict with the gene groups generated with parse_david_gene_groups_file()
    I have added manual functional annotation to the david output that can also be added as a column
    """

    drosophila_attributes_dict = gff.parse_gff3_for_attributes(drosophila_gff_path)
    not_found_sig_IDs = []
    unassigned_FB_IDs = [] # count IDs not in DAVID analysis if that is included
    # print(f"{len(drosophila_attributes_dict)} transcripts in drosophila")

    CAFE_results = {}
    if CAFE_results_path != "":
        CAFE_results = OGs.parse_CAFE_output(CAFE_results_path)

    if orthogroups_dict_all != {}:
        OGs_all_list = list(orthogroups_dict_all.keys())
    
    if species_tree == "":
        # get complete species list:
        species_list = []
        for OG_id in OGs_all_list:
            species_list.extend(list(orthogroups_dict_all[OG_id].keys()))
        species_list = list(set(species_list))
        print(f"Species included in the input file: {species_list}")
        # assert "T_molitor" in species_list
    else:
        species_names_unsorted = my_plotting.plot_tree_manually(species_tree)
        species_coords_sorted = sorted(list(species_names_unsorted.keys()))
        species_list = [species_names_unsorted[species_coord] for species_coord in species_coords_sorted]
        print(f"Species sorted according to phylogenetic tree: {species_list}")


    flybase_url = "https://api.flybase.org/api/v1.0/gene/summaries/auto/" ## add gene ID afterwards

    # write headers according to input files
    outfile_name = f"/Users/miltr339/work/PhD_code/PhD_chapter1/data/{outfile_name}"
    with open(outfile_name, "w") as outfile:
        if orthogroups_dict_all =={} and david_gene_groups == {} and david_functions == {}:
            outfile.write("Orthogroup_ID\ttranscript_ID_native\tFlybase\tFlybase_summary\tCAFE_p-value\n")
        elif orthogroups_dict_all !={} and david_gene_groups == {} and david_functions == {}:
            species_header = "\t".join([f"{species}" for species in species_list])
            outfile.write(f"Orthogroup_ID\t{species_header}\tmax_delta_GF\ttranscript_ID_native\tFlybase\tFlybase_summary\tCAFE_p-value\n")
        elif orthogroups_dict_all !={} and david_gene_groups != {} and david_functions == {}:
            species_header = "\t".join([f"{species}" for species in species_list])
            outfile.write(f"Orthogroup_ID\tGene_Group\tGene_Name\t{species_header}\tmax_delta_GF\ttranscript_ID_native\tFlybase\tFlybase_summary\tCAFE_p-value\n")
        elif orthogroups_dict_all !={} and david_gene_groups != {} and david_functions != {}:
            species_header = "\t".join([f"{species}" for species in species_list])
            outfile.write(f"Orthogroup_ID\tGene_Group\tGroup_function\tGene_Name\t{species_header}\tmax_delta_GF\ttranscript_ID_native\tFlybase\tFlybase_summary\tCAFE_p-value\n")


        for OG_id, transcripts_list in tqdm(orthogroup_dict_species.items()):
            # for weird parsing stuff i did like a year ago the transcript IDs in the native drosophila annotation have leading "__" that should be removed
            # also remove the tailing "_1" 
            transcripts_list = [transcript.replace("__", "")[:-2] for transcript in transcripts_list]
            
            if OGs_list != [] and OG_id not in OGs_list:
                continue
            
            if orthogroups_dict_all != {}:
                OG_species = orthogroups_dict_all[OG_id]
                GF_size = []
                for species in species_list:
                    try:
                        GF_size.append(len(OG_species[species]))
                    except:
                        GF_size.append(int(0))
                
                delta_GF = int(max(GF_size)) - int(min(GF_size))
                GF_size = [str(size) for size in GF_size]
                GF_size = "\t".join(GF_size) + f"\t{delta_GF}"
            
            for transcript in transcripts_list:
                try:
                    attributes = drosophila_attributes_dict[transcript]
                except:
                    raise RuntimeError(f"the transcript ID {transcript} listed in the orthofinder output does not appear in the annotation")
                    attributes = {}
                    not_found_sig_IDs.append(transcript)
                
                flybase = attributes["Dbxref"].split(":FBgn")[-1]
                flybase = f"FBgn{flybase}".split(",")[0]
                if get_gene_functions_from_API:

                    response = req.get(f"{flybase_url}{flybase}")

                    if response.status_code == 200:
                        api_data = response.json()
                        flybase_summary = api_data["resultset"]["result"][0]["summary"]
                        flybase_summary = flybase_summary.strip() # sometimes they end with a tab to fuck with me
                        # print(flybase_summary)
                    else:
                        # print(f"{flybase_url}{flybase}")
                        flybase_summary = "None"
                        # raise RuntimeError(f"flybase didn't work for {OG_id}, {flybase}")
                else:
                    flybase_summary = "None"

                try:
                    cafe_p = CAFE_results[OG_id]
                except:
                    cafe_p = "None"

                # Orthogroup_ID transcript_ID_native Flybase Flybase_summary CAFE_p-value max_delta_GF
                if orthogroups_dict_all =={} and david_gene_groups =={} and david_functions == {}:
                    outfile_string = f"{OG_id}\t{transcript}\t{flybase}\t{flybase_summary}\t{cafe_p}\n"
                elif orthogroups_dict_all !={} and david_gene_groups =={} and david_functions == {}:
                    outfile_string = f"{OG_id}\t{transcript}\t{flybase}\t{flybase_summary}\t{cafe_p}\t{GF_size}\n"
                elif orthogroups_dict_all !={} and david_gene_groups !={}:
                    try:
                        gene_group, gene_name = david_gene_groups[flybase]
                        gene_group = ",".join(gene_group)
                    except:
                        gene_group = "None"
                        gene_name = "None"
                        unassigned_FB_IDs.append(flybase)
                        # try to parse the gene name from the flybase API summary
                        if "The gene" in flybase_summary[0:10]:
                            gene_name =flybase_summary.split("The gene ")[-1]
                            gene_name = gene_name.split(" is referred to")[0]
                            gene_name = f"{gene_name} (from API summary)"
                    try:
                        group_function = david_functions[gene_group]
                    except:
                        group_function = "None"
                    # Orthogroup_ID Gene_Group Gene_Name {species_header} transcript_ID_native Flybase Flybase_summary CAFE_p-value max_delta_GF 
                    outfile_string = f"{OG_id}\t{gene_group}\t{group_function}\t{gene_name}\t{GF_size}\t{transcript}\t{flybase}\t{flybase_summary}\t{cafe_p}\n"
                outfile.write(f"{outfile_string}")
    if len(unassigned_FB_IDs):
        print(f"{len(unassigned_FB_IDs)} Flybase Gene IDs in the orthogroups were not found in David")
    print(f"flybase IDs written to {outfile_name} in the data directory.")
    
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
            if not transcript:
                continue
            transcript = f"{species}_{transcript}"
            try: 
                record = proteinfasta[transcript]
            except:
                # continue ## TODO it seems like some of the transcripts in here should have been removed by TE_filtering?
                raise RuntimeError(f"{transcript} not found in reference proteinfasta file {proteinfasta_reference}")
            
            record.id = f"{record.id}_{OG_id}"
            filtered_fasta.append(record)

    SeqIO.write(filtered_fasta, outfile_name, "fasta")
    print(f"proteinfasta written to {outfile_name}")
    return filtered_fasta


def parse_blast_outfile(blast_filepath, min_seq_ident = 90):
    """
    makes a dictionary with only transcript IDs like:
    {
        orthoDB_query1 : [ native_hit1 , native_hit2 , ... ] ,
        orthoDB_query2 : [ native_hit1 , native_hit2 , ... ] ,
    }
    """

    out_dict = {}
    with open(blast_filepath, "r") as blast_file:
        for line in blast_file.readlines():
            line = line.strip().split("\t")
            
            orthoDB_query = line[0]
            OG_id = orthoDB_query.split("_")[-1]
            native_ref = line[1]
            native_trans_ID = native_ref.replace("D_melanogaster_", "")
            #native_trans_ID = native_trans_ID[:-2] #remove tailing "_1"
            seq_ident = line[2]
            if float(seq_ident) < float(min_seq_ident):
                continue
            
            try:
                out_dict[OG_id].append(native_trans_ID)
            except:
                out_dict[OG_id] = [native_trans_ID]
    
    return(out_dict)


def parse_david_gene_groups_file(david_gene_groups_filepath:str):
    """
    parse the gene groups output from DAVID functional analysis to sort the flybase IDs into functional gene groups
    out_dict = {
        flybase_ID : [ [Gene_group_1] ,  Gene_name ],
        flybase_ID : [ [Gene_group_1, Gene_group_2] ,  Gene_name ],
        ...
    }
    """
    out_dict = {}
    count_gene_groups = 0
    with open(david_gene_groups_filepath, "r") as gene_groups_file:
        current_gene_group = ""
        for line_string in gene_groups_file.readlines():
            line_string = line_string.strip()
            if not line_string:
                continue
            line = line_string.split("\t")
            if "Gene Group" in line[0]:
                current_gene_group = line[0]
                count_gene_groups += 1
            elif "FLYBASE_GENE_ID" in line[0]:
                continue
            elif "FBgn" in line[0]:
                flybase_ID = line[0]
                if flybase_ID not in out_dict:
                    out_dict[flybase_ID] = [[current_gene_group], line[1]]
                else:
                    out_dict[flybase_ID][0].append(current_gene_group)

            else:
                raise RuntimeError(f"Line: '{line_string}' could not be parsed")
    print(f"{count_gene_groups} Gene Groups in file.")
    return out_dict


def parse_david_group_functions(david_gene_groups_filepath:str):
    """
    I manually functionally annotated the gene groups file and here i will parse those into a dict
    { gene group : functional annotation }
    """
    out_dict = {}
    count_gene_groups = 0
    with open(david_gene_groups_filepath, "r") as gene_groups_file:
        for line_string in gene_groups_file.readlines():
            line_string = line_string.strip()
            if not line_string:
                continue
            line = line_string.split("\t")
            if "Gene Group" not in line[0]:
                continue
            current_gene_group = line[0]
            count_gene_groups += 1
            if "Function: " in line[-1]:
                try:
                    function_String, function_annot = line[-1].split("Function: ")
                except:
                    raise RuntimeError(f"{line[-1]} not correctly parsed by function")
                
                out_dict[current_gene_group] = function_annot
    return out_dict
            


if __name__ == "__main__":

    orthoDB_annotations, orthogroups_orthoDB, sig_orthoDB, orthoDB_proteinseqs = filepaths_orthoDB()
    native_annotations, orthogroups_native, sig_native, native_proteinseqs = filepaths_native()
    dmel_unfiltered_annot = "/Users/miltr339/work/native_annotations/d_melanogaster_NOT_isoform_filtered.gff"    
    david_gene_groups_path = "/Users/miltr339/Box Sync/thesis writing/Milena chapter1/Sig OG Flybase IDs/DAVID-FunctionalClustering_FBgenes.txt"
    david_gene_groups_path_home = "/Users/milena/Box Sync/thesis writing/Milena chapter1/Sig OG Flybase IDs/DAVID-FunctionalClustering_FBgenes.txt"
    tree_path = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_native/SpeciesTree_native_only_species_names.nw"
    ## Get stuff from native with functional annotations
    if False:
        print(f"\n\tnative")
        native_sig_list, native_all_list =OGs.get_sig_orthogroups(sig_native)
        native_sig_OGs_dict = OGs.parse_orthogroups_dict(orthogroups_native, sig_list = native_sig_list, species="D_melanogaster")
        if len(native_sig_list) != len(native_sig_OGs_dict):
            print(f"{len(native_sig_list)} significant orthogroups : {native_sig_list[0:5]}...")
            print(f"{len(native_sig_OGs_dict)} orthogroups in sig dict")

        native_sig_all_species = OGs.parse_orthogroups_dict(orthogroups_native, sig_list = native_sig_list)
        native_large_OGs = OGs.get_orthogroup_sizes(native_sig_all_species, q=0)
        large_OG_IDs = list(native_large_OGs.keys())

        num_transcripts = 0
        for og, tr_list in native_sig_OGs_dict.items():
            num_transcripts += len(tr_list)

        not_found_transcripts = get_flybase_IDs(native_sig_OGs_dict, native_annotations["D_melanogaster"], OGs_list=large_OG_IDs, orthogroups_dict_all = native_sig_all_species, CAFE_results_path = sig_native)
        if len(not_found_transcripts)>0:
            print(f"{len(not_found_transcripts)} (of {num_transcripts}) transcripts from orthoDB not found in annotation: {not_found_transcripts}")


    # get stuff for orthoDB annotations
    if True:
        print(f"\n\torthoDB")
        orthoDB_sig_list, orthoDB_all_list =OGs.get_sig_orthogroups(sig_orthoDB)
        orthoDB_sig_OGs_dict = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_list = orthoDB_sig_list, species="D_melanogaster")
        if len(orthoDB_sig_list) != len(orthoDB_sig_OGs_dict):
            print(f"{len(orthoDB_sig_list)} significant orthogroups : {orthoDB_sig_list[0:5]}...")
            print(f"{len(orthoDB_sig_OGs_dict)} orthogroups in sig dict")

        orthoDB_sig_all_species = OGs.parse_orthogroups_dict(orthogroups_orthoDB, sig_list = orthoDB_sig_list)
        orthoDB_large_OGs = OGs.get_orthogroup_sizes(orthoDB_sig_all_species, q=0)
        large_OG_IDs = list(orthoDB_large_OGs.keys())
        if len(orthoDB_sig_all_species) != len(large_OG_IDs):
            print(len(orthoDB_sig_all_species))
            print(len(large_OG_IDs))
        
        # make_proteinfasta_from_orthogroup(orthoDB_sig_OGs_dict, orthoDB_proteinseqs["D_melanogaster"], orthogroups_to_include=large_OG_IDs)

        # blast the orthoDB proteins against the native ones
        """
        BLASTp search of the above created proteins against the native Dmel annotation to get the flybase IDs
       
        makeblastdb -in D_melanogaster.faa -dbtype prot     # native annotation for reference db
        blastp -query PhD_chapter1/Dmel_transcripts_from_sig_OGs.fasta -db /Users/miltr339/work/native_proteinseqs/D_melanogaster.faa -out /Users/miltr339/work/PhD_code/Dmel_oDB_vs_nat.out -outfmt 6 -num_threads 5 -evalue 1e-10
        """
        if True:
            blast_outfile = "/Users/miltr339/work/PhD_code/Dmel_oDB_vs_nat.out"
            blast_out_dict = parse_blast_outfile(blast_outfile)

            num_transcripts = 0
            for og, tr_list in blast_out_dict.items():
                num_transcripts += len(tr_list)

            david_gene_groups_dict = parse_david_gene_groups_file(david_gene_groups_path)
            david_gene_groups_function = parse_david_group_functions(david_gene_groups_path)

            ## TODO fix this somehow 
            # orthoDB_sig_all_species = OGs.parse_orthogroups_dict(orthogroups_orthoDB)
            not_found_transcripts = make_table_with_flybase_functions(blast_out_dict, dmel_unfiltered_annot, outfile_name = "orthoDB_sig_OGs_flybase_IDs_with_group_function.tsv", 
                orthogroups_dict_all = orthoDB_sig_all_species, 
                CAFE_results_path = sig_orthoDB, 
                get_gene_functions_from_API=False,
                david_gene_groups = david_gene_groups_dict, 
                david_functions = david_gene_groups_function, 
                species_tree = tree_path
            )

            if len(not_found_transcripts)>0:
                print(f"{len(not_found_transcripts)} (of {num_transcripts}) transcripts from orthoDB not found in annotation: {not_found_transcripts}")


