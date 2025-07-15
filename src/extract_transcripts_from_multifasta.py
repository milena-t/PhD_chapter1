"""
extract a fasta file with all transcripts in an input list from a larger multifasta file
"""
from Bio import SeqIO
from parse_gff import read_dict_from_file

in_list = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/functional_annot_eval/aobt_expansion_GF.txt"
in_fasta = "/Users/miltr339/work/orthoDB_proteinseqs_TE_filtered/A_obtectus_filtered_proteinfasta_TE_filtered.fa"
out_path = "/Users/miltr339/work/orthofinder/orthoDB_uniform_masking_orthogroups/"


def filter_fasta_by_header(in_fasta:str, headers:list[str], out_fasta:str):
    out_fasta_list = []
    for record in SeqIO.parse(in_fasta, "fasta"):
        
        if record.id in headers:
            out_fasta_list.append(record)

    SeqIO.write(out_fasta_list, out_fasta, "fasta")
    print(f" --> filtered proteinfasta written to {out_fasta}")
    return out_fasta


if __name__ == "__main__":
    
    orthogroups_dict = read_dict_from_file(in_list)

    for orthogroup, transcript_list in orthogroups_dict.items():
        filter_fasta_by_header(in_fasta=in_fasta, headers=transcript_list, out_fasta=f"{out_path}Aobt_{orthogroup}_proteins.fasta")


 