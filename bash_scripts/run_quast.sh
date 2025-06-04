#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 10:00
#SBATCH -J quast
#SBATCH -o quast.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools
module load python3
module load quast/5.0.2

QUAST_PATH=/proj/naiss2023-6-65/Milena
RAGTAG_ANALIS_PATH=/proj/naiss2023-6-65/Milena/ragtag_assembly_improvement/ragtag_output_analis
RAGTAG_ASSEMBLY=ragtag.scaffold.fasta

NCBI_DIR=/proj/naiss2023-6-65/Milena/coleoptera_sequences/sequence_downloads
OUTPUT_DIR=/proj/naiss2023-6-65/Milena/data_processing/assembly_evaluation_scripts_and_results

W_DIR=$(pwd)

# for GCA_DIR in $(echo "GCA_000355655.1_ncbi_download GCA_000390285.2_ncbi_download GCA_000500325.2_ncbi_download GCA_000648695.2_ncbi_download GCA_000699045.2_ncbi_download GCA_001412225.1_ncbi_download GCA_002938485.2_ncbi_download GCA_003013835.2_ncbi_download GCA_004193795.1_ncbi_download GCA_008802855.1_ncbi_download GCA_011009095.1_ncbi_download GCA_014462685.1_ncbi_download GCA_015345945.1_ncbi_download GCA_020466585.2_ncbi_download GCA_020466635.2_ncbi_download GCA_027724725.1_ncbi_download GCA_027725215.1_ncbi_download GCA_029955315.1_ncbi_download GCA_963669975.1_ncbi_download"); \
# do \
# cd $NCBI_DIR/$GCA_DIR/ncbi_dataset/data/${GCA_DIR%_ncbi_download} ; \
# echo "\n ==================================== \n"
# echo $GCA_DIR/ncbi_dataset/data/${GCA_DIR%_ncbi_download} ; \
# ASSEMBLY_FASTA=$(ls GCA*)
# echo $ASSEMBLY_FASTA
# 
# mkdir $OUTPUT_DIR/${ASSEMBLY_FASTA%_genomic.fna}
# echo $OUTPUT_DIR/${ASSEMBLY_FASTA%_genomic.fna}
# 
# python3 $QUAST_PATH/quast.py -o $OUTPUT_DIR/${ASSEMBLY_FASTA%_genomic.fna} -t 4 $ASSEMBLY_FASTA
# 
# cd $NCBI_DIR ; \
# done
# 
# cd $W_DIR 



# echo "analis quast"
# python3 /proj/naiss2023-6-65/Milena/quast.py -o analis_quast -t 10 /proj/naiss2023-6-65/Milena/coleoptera_sequences/c_analis/analis_from_uppmax.fasta
# echo "analis quast done"
# echo "chinensis quast"
# python3 /proj/naiss2023-6-65/Milena/quast.py -o chinensis_quast -t 10 /proj/naiss2023-6-65/Milena/coleoptera_sequences/c_chinensis/chinensis_from_uppmax.fasta
# echo "chinensis quast done"

CMAC_superscaffolded=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus_superscaffolded/assembly_genomic.fna.masked
OUT_DIR_CMAC=/proj/naiss2023-6-65/Milena/data_processing/assembly_evaluation_scripts_and_results/CMAC_superscaffolded
python3 /proj/naiss2023-6-65/Milena/quast.py -o $OUT_DIR_CMAC -t 10 $CMAC_superscaffolded

# python3 /proj/naiss2023-6-65/Milena/quast.py -o maculatus_from_ENA_quast -t 4 /proj/naiss2023-6-65/Milena/coleoptera_sequences/c_maculatus/Cmac_from_ENA_GCA_951848785.1.fasta.masked
# python3 /proj/naiss2023-6-65/Milena/quast.py -o maculatus_from_uppmax_quast -t 4 /proj/naiss2023-6-65/Milena/coleoptera_sequences/c_maculatus/C_mac_Ys_hifiasm.fasta
# python3 /proj/naiss2023-6-65/Milena/quast.py -o obtectus_from_ENA_quast -t 4 /proj/naiss2023-6-65/Milena/coleoptera_sequences/a_obtectus/a_obtectus_ENA.fasta
# python3 /proj/naiss2023-6-65/Milena/quast.py -o obtectus_from_uppmax_quast -t 4 /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/A_obtectus_old/assembly_genomic.fna