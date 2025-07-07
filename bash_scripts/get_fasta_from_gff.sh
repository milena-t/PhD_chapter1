#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J make_filtered_transcripts
#SBATCH -o make_filtered_transcripts.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

# use gffread to extract the protein coding sequences
# -M :  cluster the input transcripts into loci, discarding "duplicated" transcripts (those with the same exact introns and fully contained or equal boundaries)
# -x :  write a FASTA file with spliced CDS for each GFF transcript

module load bioinfo-tools gffread/0.12.7 samtools/1.20 emboss/6.6.0

# ANNOT_DIRS=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation
ANNOT_DIRS=/proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lome_superscaffolded_comparison/Cmac_*

for ANNOT_DIR in $(echo $ANNOT_DIRS) # add string to specify species (or species subsets) around wildcard
do 
    cd $ANNOT_DIR
    ANNOT_GTF="${ANNOT_DIR}/Cmac_*/braker/braker.gtf"
    FILTERED_GTF="${ANNOT_GTF%.*}_isoform_filtered.gff" # originally gtf but keep_longest_isoform.pl automatically returns gff version 3 
    ANNOT_TRANSCRIPTS=isoform_filtered_transcripts.fna
    ANNOT_PROTEINS=isoform_filtered_proteins.faa
    ASSEMBLY=assembly_genomic.fna.masked

    # braker doesn't like spaces in the contig names so it replaces them with underscores. 
    # The below command makes them match the first column in the gtf file again.
    # run only once
    # sed -i 's/ /_/g' $ASSEMBLY 

    # for the callosobruchuses there is an additional replacement necessary
    # In analis and chinensis, the fasta headers look like this: >31|quiver but the annotation looks for this 31_quiver
    # maculatus has a longer header but also some "|" that are replaced by braker in the annotation 
    # also run only once
    #sed -i 's/|/_/g' $ASSEMBLY

    echo $(pwd)
    echo $ASSEMBLY

    SPECIES_NAME=$(basename "$ANNOT_DIR")
    echo $SPECIES_NAME

    # index assemblies (greatly decreases computing time, and won't work for the more fragmented callosobruchus assemblies otherwise)
    samtools faidx $ASSEMBLY
    # extract transcript sequences
    gffread -M -x $ANNOT_TRANSCRIPTS -g $ASSEMBLY $ANNOT_GTF
    # change fasta headers to include species names
    sed -i "s/>/>${SPECIES_NAME}_/g" $ANNOT_TRANSCRIPTS
    # translate transcript sequences
    transeq -sequence $ANNOT_TRANSCRIPTS -outseq $ANNOT_PROTEINS


    ls -lh $ANNOT_TRANSCRIPTS
    echo "###########################################"

done

# get number of transcripts in all output files:
# for transcripts in $(echo */isoform_filtered_transcripts.faa) ; do echo $transcripts ; grep ">" $transcripts | wc -l ; done

# rename transcripts and link them in a folder for orthofinder
# for SPECIES in /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/*
# do 
#     SPECIES_NAME=$(basename "$SPECIES")
#     sed -i "s/>/>${SPECIES_NAME}_/g" $SPECIES_NAME/isoform_filtered_transcripts.faa 
#     echo "done $SPECIES_NAME"
#     ln -s $SPECIES/isoform_filtered_transcripts.faa /proj/naiss2023-6-65/Milena/gene_family_analysis/proteinseqs_braker_all_species_annot/${SPECIES_NAME}_transcripts.fa 
# done

# for SPECIES in /proj/naiss2023-6-65/Milena/annotation_pipeline/all_proteinrefs_annotation/annotation_species/* ; do SPECIES_NAME=$(basename "$SPECIES") ; ln -s $SPECIES/isoform_filtered_transcripts.faa /proj/naiss2023-6-65/Milena/gene_family_analysis/proteinseqs_braker_all_species_annot/${SPECIES_NAME}_transcripts.fa ; done