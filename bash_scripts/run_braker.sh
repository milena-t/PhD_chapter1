[milenatr@rackham3 annotation_pipeline]$ cat braker3_singularity_with_RNAseq_in_SNIC_TMP.sh 
#!/bin/bash -l

# example run:
# sbatch --job-name="C_septempunctata" --output="C_septempunctata_braker.out" /proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_all_species_proteinseqs.sh C_septempunctata assembly_genomic.fna
# sbatch --job-name="C_maculatus" --output="C_maculatus_Lu2024_braker.out" braker3_singularity_all_species_proteinseqs.sh C_maculatus_Lu2024 /proj/naiss2023-6-65/Milena/annotation_pipeline/Cmac_Lu2024_comparison/braker/GCA_040182625.1_Cmac_2024_genomic.fna
# sbatch --job-name="D_ponderosae" --output="D_ponderosae_braker.out" /proj/naiss2023-6-65/Milena/annotation_pipeline/braker3_singularity_all_species_proteinseqs.sh D_ponderosae /proj/naiss2023-6-65/Milena/coleoptera_sequences/d_ponderosae/GCA_020466635.2_Dpon_M_20191212v2_genomic.fna


#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J braker_%j
#SBATCH -o braker_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user milena.trabert@ebc.uu.se

module load bioinfo-tools
    # be careful to not have any other braker dependency modules or something loaded when running the container, because they go with different perl versions
    # If the perl version is incompatible with the one in the container there will be a fatal error
    # the safest option is probably to not have any modules loaded at all



if [ $# -lt 2 ]; then
    echo "Usage: $0 name_of_species_dir path_to_masked_assembly."
    echo "you have $#"
    exit 1
fi



SPECIES=$1
ASSEMBLY_MASKED=$2
PROTEIN_DATA=$3
# $4 is optional for RNA reference data by SRR number 


# run the script from this directory. It's species-specific
# export wd=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/$SPECIES
export home_wd=${PWD}/${SPECIES}

if [ -d ${home_wd} ]; then
    echo "Working directory already exists: ${home_wd}"
else
    mkdir $home_wd
    echo "created directory: $home_wd"
fi

export wd=${SNIC_TMP}/${SPECIES}
#ASSEMBLY_MASKED=/proj/naiss2023-6-65/Milena/coleoptera_sequences/c_chinensis/chinensis_from_uppmax.fasta.masked

if [ -d ${wd} ]; then
    echo "Working directory in temporary directory already exists: ${wd}"
else
    mkdir $wd
    echo "created directory in SNIC_TMP: $wd"
fi

cd $wd



# link braker.sif file
ln -s /proj/naiss2023-6-65/Milena/annotation_pipeline/braker3.sif braker3.sif

# check if the augustus_config direcotry exists,
export AUGUSTUS_CONFIG_PATH=${wd}/augustus_config
if [ -d ${wd}/augustus_config ]; then
    echo "Augustus_config already exists: ${wd}/augustus_config/species"
else
    echo "Augustus config does not exist, create it and change write permissions"
    module load augustus/3.5.0-20231223-33fc04d # so that the source command works
    source $AUGUSTUS_CONFIG_COPY
    chmod a+w -R ${wd}/augustus_config/species
    module unload augustus/3.5.0-20231223-33fc04d # same as above, some weird shit with conflicting perl versions
    echo "augustus config path in the function: ${AUGUSTUS_CONFIG_PATH}"
fi

# export PROTEIN_REF_ALL_SPECIES=/proj/naiss2023-6-65/Milena/annotation_pipeline/all_proteinrefs_annotation/orthoDB_and_species_proteins.fa

export ETP=/sw/bioinfo/GeneMark-ETP/1.02-20231213-dd8b37b/rackham/bin
# the braker example for using the container references this variable in the GENEMARK_PATH flag,
# Just from the name I assume it's genemark-ETP and not ES

# there should not already be an existing braker output directory in the working directory, otherwise there will be an error that it can't create the genemark-es ouptut file
if [ -d ${wd}/braker ]; then
  rm -r ${wd}/braker
  echo "removed preexisting output directory at: ${wd}/braker"
else
  echo "no existing directory at: ${wd}/braker, proceed"
fi


# old header:
# singularity exec -B ${PWD}:${PWD} braker3.sif braker.pl \
# new: Bind the working directory to ensure it's accessible within the container
if [ $# -eq 4 ]; then
    echo "You have included a third command line argument that is assumed to contain SRA-ids for species-specific RNAseq data"
    SRA_IDS=$4
    singularity exec -B ${wd}:${wd} braker3.sif braker.pl \
        --genome=${ASSEMBLY_MASKED} \
        --prot_seq $PROTEIN_DATA \
        --rnaseq_sets_ids=$SRA_IDS \
        --threads 20 \
        --GENEMARK_PATH=${ETP}/gmes \
        --AUGUSTUS_CONFIG_PATH=${wd}/augustus_config \
        --useexisting
else
    singularity exec -B ${wd}:${wd} braker3.sif braker.pl \
        --genome=${ASSEMBLY_MASKED} \
        --prot_seq $PROTEIN_DATA \
        --threads 20 \
        --GENEMARK_PATH=${ETP}/gmes \
        --AUGUSTUS_CONFIG_PATH=${wd}/augustus_config \
        --useexisting
fi



echo "move SNIC_TMP directory to our storage"

if [ -d ${home_wd}/braker ]; then
  rm -r ${home_wd}/braker
  echo "removed preexisting output directory at: ${home_wd}/braker"
else
  echo "no existing directory at: ${home_wd}/braker"
fi
mv $wd/braker $home_wd/braker

echo "done!"