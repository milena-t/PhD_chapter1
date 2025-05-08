#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:30:00
#SBATCH -J eval_annot_agat
#SBATCH -o annotation_evaluation_orthoDB_with_AGAT.log


module load bioinfo-tools AGAT/1.3.2

ANNOT_DIRS=$(echo /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/*)

for ANNOT_DIR in $ANNOT_DIRS
do
    echo "############"
    echo $ANNOT_DIR
    ANNOT_GTF=$ANNOT_DIR/braker/braker.gtf
    FILTERED_GTF=$ANNOT_DIR/braker/braker_isoform_filtered.gff
    
    agat_sp_statistics.pl -gff $ANNOT_GTF -o $ANNOT_DIR/braker/agat_unfiltered_stats.txt
    agat_sp_statistics.pl -gff $FILTERED_GTF -o $ANNOT_DIR/braker/agat_filtered_stats.txt

done

# extract relevant information


ANNOT_DIRS=$(echo /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/*)
echo "unfiltered_stats"> annot_stats_agat.txt
echo "filtered_stats" >> annot_stats_agat.txt

for ANNOT_DIR in $ANNOT_DIRS
do
    echo "############" >> annot_stats_agat.txt
    echo $ANNOT_DIR >> annot_stats_agat.txt
    UNFILTERED_STATS=$ANNOT_DIR/braker/agat_unfiltered_stats.txt
    FILTERED_STATS=$ANNOT_DIR/braker/agat_filtered_stats.txt
    
    grep "Number of gene" -m 1 $UNFILTERED_STATS >> annot_stats_agat.txt # only get first occurence since it repeats all the statistics for only the longest isoforms again
    grep "Number of gene" -m 1 $FILTERED_STATS >> annot_stats_agat.txt
    grep "Number of transcript" -m 1 $UNFILTERED_STATS >> annot_stats_agat.txt
    grep "Number of transcript" -m 1 $FILTERED_STATS >> annot_stats_agat.txt
    grep "mean exons per transcript" -m 1 $UNFILTERED_STATS >> annot_stats_agat.txt
    grep "mean exons per transcript" -m 1 $FILTERED_STATS >> annot_stats_agat.txt
    grep "mean transcript length (bp)"  -m 1 $UNFILTERED_STATS >> annot_stats_agat.txt
    grep "mean transcript length (bp)"  -m 1 $FILTERED_STATS >> annot_stats_agat.txt
    grep "mean exon length (bp)"  -m 1 $UNFILTERED_STATS >> annot_stats_agat.txt
    grep "mean exon length (bp)"  -m 1 $FILTERED_STATS >> annot_stats_agat.txt
    grep "Number gene overlapping" -m 1 $UNFILTERED_STATS >> annot_stats_agat.txt
    grep "Number gene overlapping" -m 1 $FILTERED_STATS >> annot_stats_agat.txt
    grep "Number of single exon gene" -m 1 $UNFILTERED_STATS >> annot_stats_agat.txt
    grep "Number of single exon gene" -m 1 $FILTERED_STATS >> annot_stats_agat.txt

done