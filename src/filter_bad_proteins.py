from Bio import SeqIO
import sys

# for the gene family analysis, you should filter out proteins that don't start with M or have stop-codons in the middle
# only works for proteinfasta! 

# for use on uppmax:
# module load bioinfo-tools augustus/3.3.3-CGP Biopython/1.73-foss-2019a

def filter_bad_proteins(input_proteinfasta, output_filtered_proteinfasta):
    # initialize variables
    count_non_start_M = 0
    count_internal_stop = 0
    count_good_records = 0
    filtered_records = [] 

    # loop through each sequence in the input multifasta file
    for record in SeqIO.parse(input_proteinfasta, "fasta"):
        transcript = record.seq
        transcript_modif = transcript[0:len(transcript)-1] # cut off last stop codon

        # count transcripts that don't start with M
        if transcript[0]!="M":
            count_non_start_M += 1
        # count transcripts that have an internal stop codon
        if "*" in transcript_modif:
            count_internal_stop += 1

        # write good records to good-records-list
        if transcript[0]=="M" and "*" not in transcript_modif: 
            filtered_records.append(record)
            count_good_records +=1
    
    # print output stats
    print(input_proteinfasta)
    print("non-M start codons: ", count_non_start_M, "("+str(round(count_non_start_M/count_good_records, 4))+ " %)")
    print("internal stop codons: ", count_internal_stop, "("+str(round(count_internal_stop/count_good_records, 4))+ " %)")
    print("number of normal records: ", count_good_records)
    print(" ----------------- ")
    # write good-records-list to output file
    SeqIO.write(filtered_records, output_filtered_proteinfasta, "fasta")
    return(count_non_start_M, count_internal_stop, count_good_records)      


if len(sys.argv) == 3:
    filter_bad_proteins(sys.argv[1], sys.argv[2])
else:
    print("python3 filter_bad_proteins.py input_proteinfasta.fna output_filtered_proteinfasta.fna")