#!/bin/bash -l
# Compare the position of genes between native and orthoDB annotations of the same species. 
# check for overlap because of this scenario:
# 
# native gene:      --------------*******----------------------***********---------------   (single multi-exon gene with exons *** and introns ---)
# 
# orthoDB genes:                  *******                      ***********                  (three single-exon genes)
# 
# If the orthoDB annotation has many single-exon genes because of this, then native genes have on average more than one genes on their position in the orthoDB annotation

#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 30:00
#SBATCH -J gene_overlap
#SBATCH -o gene_overlap.log

module load bioinfo-tools
module load BEDTools/2.31.1

# native_path=/Users/milena/work/a_obtectus/a_obtectus_native_isoform_filtered.gff
# orthoDB_path=/Users/milena/work/a_obtectus/a_obtectus_orthodb_isoform_filtered_scaffold_name_corrected.gff # identical to check if it works
# outfile=/Users/milena/work/a_obtectus/obtectus_comparison.tsv 

# native_path=/Users/milena/work/c_maculatus/cmac_Lu2024_native_isoform_filtered.gff
# orthoDB_path=/Users/milena/work/c_maculatus/cmac_Lu2024_orthoDB_isoform_filtered_scaffold_name_corrected.gff
# outfile=/Users/milena/work/c_maculatus/cmac_Lu2024_overlap.tsv

native_path=$2
orthoDB_path=$1
outfile=$3

# native_path=/Users/milena/work/c_maculatus/cmac_Lu2024_native_isoform_filtered.gff
# orthoDB_path=/Users/milena/work/c_maculatus/cmac_Lu2024_orthoDB_isoform_filtered_scaffold_name_corrected.gff
# outfile=/Users/milena/work/c_maculatus/cmac_Lu2024_overlap.tsv


echo "native_path: ${native_path}"
echo "orthoDB_path: ${orthoDB_path}"
echo "outfile: ${outfile}"

tmp_native=$(basename "$native_path")
tmp_native="${tmp_native%.*}.tmp.gff"
tmp_orthodb=$(basename "$orthoDB_path")
tmp_orthodb="${tmp_orthodb%.*}.tmp.gff"

echo "temporary native: ${tmp_native}"
echo "temporary orthoDB: ${tmp_orthodb}"

# grep "\tgene\t" $native_path > ${tmp_native}.2 
# awk -F'\t' '$4 ~ /^[0-9]+$/' ${tmp_native}.2 > $tmp_native # apparently there's a line with a negative start position for ID=snap_masked-@001138F_arrow_arrow-processed-gene-0.54
# grep "\tgene\t" $orthoDB_path > ${tmp_orthodb}.2
# awk -F'\t' '$4 ~ /^[0-9]+$/' ${tmp_orthodb}.2 > $tmp_orthodb


### prepare input files:
# filter only for the transcripts since otherwise you compute the overlap stats of every line because it ignores nested features.

# grep -P "mRNA\t" $native_path > ${tmp_native}.2  #OBS! on linux you need the -P flag for \t, on macOS it's not necessary
grep "\ttranscript\t" $native_path > ${tmp_native}.2  
# if the tag for the transcripts is not "mRNA"
if [ ! -s "${tmp_native}.2" ]; then
    # If the file size is 0, try again with "transcript"
    # grep -P "transcript\t" $native_path > ${tmp_native}.2 
    grep "transcript\t" $native_path > ${tmp_native}.2 
    echo "native transcript flag: \"transcript\""
else
    echo "native transcript flag: \"mRNA\""
fi

filesize_tmp=$(stat -c%s "${tmp_native}.2")
echo "The size of '${tmp_native}.2' is $filesize_tmp in bytes"

# some files have negative transcript positions somehow? filter them out
awk -F'\t' '$4 ~ /^[0-9]+$/' ${tmp_native}.2  > $tmp_native

filesize_tmp=$(stat -c%s "${tmp_native}")
echo "The size of '${tmp_native}' is $filesize_tmp in bytes"
head -n 20 ${tmp_native}


# filter uniform annotations the same as above, except the transcript tag is always "transcript" so no testing is necessary
# grep -P "\ttranscript\t" $orthoDB_path > ${tmp_orthodb}.2
grep "\ttranscript\t" $orthoDB_path > ${tmp_orthodb}.2
awk -F'\t' '$4 ~ /^[0-9]+$/' ${tmp_orthodb}.2 > $tmp_orthodb

filesize_tmp=$(stat -c%s "${tmp_orthodb}")
echo "The size of '${tmp_orthodb}' is $filesize_tmp in bytes"
head -n 20 ${tmp_orthodb}


### compute the number of overlapping genes with bedtools, see documentation for details.
# "Each feature in A is compared to B in search of overlaps" -> A is the query, and B is the reference.


# get the number of overlaps that native annotations have in the orthodb annotations
bedtools intersect -c -a $tmp_native -b $tmp_orthodb | awk -F'\t' '{print $3, $(NF-1), $NF}' OFS='\t' > "${outfile%.*}_numbers_only_native_query.txt" # gets the last two columns

# the other way around:
bedtools intersect -c -a $tmp_orthodb -b $tmp_native | awk -F'\t' '{print $3, $(NF-1), $NF}' OFS='\t' > "${outfile%.*}_numbers_only_orthodb_query.txt" # gets the last two columns

# get the ids and positions of all overlapping features with positions
# like this:
# native_feature 0 10 orthodb_feature_a 0 3 
# native_feature 0 10 orthodb_feature_b 6 9
bedtools intersect -a $tmp_native -b $tmp_orthodb -wa -wb  > "${outfile%.*}_complete.txt"

rm $tmp_native
rm ${tmp_native}.2
rm $tmp_orthodb
rm ${tmp_orthodb}.2