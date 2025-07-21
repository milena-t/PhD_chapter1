import pandas as pd
import re
from collections import Counter
import parse_gff as gff

# add a functional description column to the Orthogroups.GeneCount.tsv output to be able to run CAFE
# also make sure that the species tree names match the columns in the file


def read_orthogroups_file(path_to_orthogroups_file):
    orthogroups_dict = {} # dictionary containing {orthogroup_ID : [list, of, orthogroup, member, gene_IDs]}
    with open(path_to_orthogroups_file, "r") as orthogroups_file:
        for orthogroup in orthogroups_file:
            orthogroup = orthogroup.strip().split(":")
            orthogroups_dict[orthogroup[0]] = orthogroup[1].strip().split(" ")
    return orthogroups_dict


def get_tree_name_from_header(header_str, species_names):
    header_found = False
    for species in species_names:
        if species in header_str:
            header_found = True
            return(species)
    if not header_found:
        return(header_str)


def modify_orthogroups_from_tsv(orthogroups_filepath, outfile_name, species_names):
    orthogroups_df = pd.read_csv(orthogroups_filepath, sep="\t")
    outfile_path = "/".join(orthogroups_filepath.split("/")[:-1])+"/"+outfile_name

    # remove the "Total" column
    orthogroups_df.drop('Total', axis=1, inplace=True)
    
    # modify header names
    old_header_names = orthogroups_df.columns.tolist()
    species_names.append("Orthogroup")
    rename_dict = { old_name : get_tree_name_from_header(old_name, species_names) for old_name in old_header_names} 
    orthogroups_df.rename(columns=rename_dict, inplace = True)

    # add Desc column
    orthogroups_df.insert(loc = 0, column = 'Desc', value = ["(null)"] * orthogroups_df.shape[0])

    # write to file
    orthogroups_df.to_csv(outfile_name, sep='\t', index = False)
    print(f"modified tsv file for CAFE written to: {outfile_path}")


def parse_table_from_txt(orthogroups_filepath, outfile_name, species_names = ""): 
    orthogroups_dict = read_orthogroups_file(orthogroups_filepath)
    outfile_path = "/".join(orthogroups_filepath.split("/")[:-1])+"/"+outfile_name

    # make a list of dictionaries to populate the pandas dataframe
    # each dictionary has "Desc" (which is always "(null)"), "Orthogroup" (which is the orthogroup ID) and all of the species
    table_rows = []
    for orthogroup_ID, transcripts_list in orthogroups_dict.items():
        # count the number of times each species appears
        transcripts_list = [gff.split_at_second_occurrence(transcript) for transcript in transcripts_list]        
        transcripts_dict = Counter(transcripts_list)
        transcripts_dict = {species : int(num) for species, num in transcripts_dict.items()}

        #add the stuff to make CAFE work
        transcripts_dict = {"Desc" : "(null)" , "Orthogroup" : orthogroup_ID, **transcripts_dict} # do it this way to get the order right. not sure if it matters but just in case
        # print(transcripts_dict)
        table_rows.append(transcripts_dict)
        # break
    print(table_rows[0])

    table_df = pd.DataFrame(table_rows)
    # Fill NaN in numeric columns and convert to integers
    numeric_cols = table_df.select_dtypes(include=['number']).columns
    table_df[numeric_cols] = table_df[numeric_cols].fillna(0).astype(int)
    print(table_df)

    if len(species_names)>0:
        # modify header names to make sure they match the species name
        old_header_names = table_df.columns.tolist()
        species_names.append("Orthogroup")
        rename_dict = { old_name : get_tree_name_from_header(old_name, species_names) for old_name in old_header_names}
        table_df.rename(columns=rename_dict, inplace = True)

    table_df.to_csv(outfile_name, sep='\t', index = False)
    print(f"modified tsv file for CAFE written to: {outfile_name}")
    return table_df


def count_csv_items(cell):
    if pd.isna(cell):
        return 0
    return len(str(cell).split(','))

def modify_orthogroups_from_N0(orthogroups_filepath, species_names, outfile_name = "CAFE_input_orthoDB_TE_filtered.tsv"):
    orthogroups_df = pd.read_csv(orthogroups_filepath, sep="\t")
    outfile_path = "/".join(orthogroups_filepath.split("/")[:-1])
    outfile_name = outfile_path+"/"+outfile_name

    # remove the unnecessary column
    orthogroups_df.drop('Gene Tree Parent Clade', axis=1, inplace=True) # this can not be used to determine the MRCA of the gene family! see documentation for the hierarchical orthogroups output file
    orthogroups_df.drop('OG', axis=1, inplace=True) 

    orthogroups_df.iloc[:, 1:] = orthogroups_df.iloc[:, 1:].map(count_csv_items)

    # Fill NaN in numeric columns and convert to integers
    # numeric_cols = orthogroups_df.select_dtypes(include=['number']).columns
    # orthogroups_df[numeric_cols] = orthogroups_df[numeric_cols].fillna(0).astype(int)
    
    # modify header names
    old_header_names = orthogroups_df.columns.tolist()
    rename_dict = { old_name : get_tree_name_from_header(old_name, species_names) for old_name in old_header_names} # since the "OG" header is not in the species names list it will just leave that
    orthogroups_df.rename(columns=rename_dict, inplace = True)
    
    # add Desc column
    orthogroups_df.insert(loc = 0, column = 'Desc', value = ["(null)"] * orthogroups_df.shape[0])
    
    # write to file
    orthogroups_df.to_csv(outfile_name, sep='\t', index = False)
    print(f"modified tsv file for CAFE written to: {outfile_name}")


def filter_CAFE_input_file(cafe_filepath:str, outgroup_header_name = "D_melanogaster", max_GF_size = 100, outfile_path = ""):
    """
    Filter the CAFE input file to 
    * only include species present at the root of the tree (here in D. melanogaster)
    * only inlcude orthogroups where all gene families are below a certain size threshold
    """
    if outfile_path == "":
        outfile_path = cafe_filepath.replace(f".tsv", "")
        outfile_path = f"{outfile_path}_filtered.tsv"

    no_outgroup = 0
    too_big_OG = 0
    with open(cafe_filepath, "r") as cafe_infile, open(outfile_path, "w") as cafe_outfile:
        cafe_lines = cafe_infile.readlines()
        header = cafe_lines[0]
        cafe_outfile.write(header)
        
        header = header.strip().split("\t")
        outgroup_index = header.index(outgroup_header_name)
        
        for line_raw in cafe_lines[1:]:
            line = line_raw.strip().split("\t")

            if line[outgroup_index] == "0":
                # print(f"no Dmel --> {line_raw}")
                no_outgroup += 1
                continue

            if any(int(GF_size) > max_GF_size for GF_size in line[2:]):
                print(f" too big GF --> {line_raw}")
                too_big_OG += 1
                continue
            
            cafe_outfile.write(line_raw)

    print(f"{no_outgroup} OGs not at the root of the tree, \n{too_big_OG} orthogroups with gene families larger than {max_GF_size}\noutfile written to {outfile_path}")
                
        
            

            

if __name__ == '__main__':

    try:
        species_names = gff.make_species_order_from_tree("/Users/miltr339/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw")
    except:
        species_names = gff.make_species_order_from_tree("/Users/milena/Box Sync/code/annotation_pipeline/annotation_scripts_ordered/14_species_orthofinder_tree.nw")

    orthoDB_orthogroups = "/Users/milena/work/chapter1_final_postprocessing/orthofinder_uniform/N0.tsv"
    native_orthogroups = "/Users/milena/work/chapter1_final_postprocessing/orthofinder_native/N0.tsv"
    
    # modify_orthogroups_from_N0(native_orthogroups, species_names=species_names, outfile_name="CAFE_input_native_from_N0.tsv")
    # modify_orthogroups_from_N0(orthoDB_orthogroups, species_names=species_names, outfile_name="CAFE_input_orthoDB_from_N0.tsv")
    
    filter_CAFE_input_file(cafe_filepath = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/orthofinder_uniform/CAFE_input_orthoDB_from_N0.tsv", outgroup_header_name = "D_melanogaster", max_GF_size = 100, )