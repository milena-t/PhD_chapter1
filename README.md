# Uniform Genome annotations for comparative genomics

This is all the code that was used in the first chapter of my thesis where I create uniform de-novo genome annotations to investigate gene family evolution in coleoptera. I implement and evaluate a uniform pipeline and remove as much technical bias as possible. The important aspects are:

* start with repeatmodeller and repeatmasker, even if the assembly is already masked
* exclude species-specific RNAseq (only orthoDB protein evidence)
* filter annotated protein sequences in each species with 

# Pipeline

## Repeatmasking

Use repeatmodeller and repeatmasker to re-mask all assemblies. There isn't really parameters to adjust, so just run as-is, `bash_scripts/repeatmasking.sh`. Keep track of the repeat families faster that repeatmodeller gives because it's necessary for filtering later.

## Genome annotation with BRAKER

Use BRAKER3 with the orthoDB arthropoda dataset to uniformly annotate all genomes. Since I don't use RNAseq, technically BRAKER2 is run but I use the BRAKER3 docker container ([here](https://hub.docker.com/r/teambraker/braker3)). Since the Uppmax cluster doesn't allow docker containers/images, I used a singularity container to make a singularity image of the docker container (see `src/make_singularity_image_braker3`). I use the orthoDB Arthropoda proteinfasta provided by the BRAKER authors for evidence (v12 [here](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/)). 

### Running BRAKER

Since BRAKER can be a little difficult to run, I have written a wrapper script that does all the setup. It is extremely specific to using the singularity image on uppmax with the scratch storage, so likely not useful for anyone on a different system. for Uppsala University people, the script is here: `bash_scripts/run_braker.sh`, you pass four command line arguments when running it with sbatch: 
* species name
* assembly path
* proteinfasta path
* optional: list of SRR numbers of RNAseq runs (separated by commas with no space, directly in the command line not path to a separate file)


## Annotation post-processing
I do the standard QC ckecks and post-processing that is necessary to do gene family evolution analysis. You can do some evaluation statistics with AGAT before and after filtering the gff to see effects with this script: `bash_scripts/annotation_statistics.sh`.

`bash_scripts/isoform_filter_gff.sh`:
1. Resolve overlapping genes into different isoforms of the same gene
2. Only keep the longest isoform of any gene

`bash_scripts/get_fasta_from_gff.sh`
1. Index the assembly (to reduce compute hours)
2. extract transcript sequences
3. translate transcript sequences

Filter protein sequences that don't have correct start/stop codons: `bash_scripts/filter_bad_proteins.sh`.

## Filter for repeats

I do a blast search of all the proteins extracted above against the repeat families fround by repeatmodeller to identify which protein sequences are likely repeats that were misannotated: `bash_scripts/blast_repeat_families_against_proteins.sh`

## Annotation evaluation

At this point I evaluate the annotation for several aspects, which is mainly facilitated by reading it into a data structure defined in `src/parse_gff.py`. This uses the third column of the gff to 

And I compare it to native annotations for the same assembly from NCBI
* `plotting/compare_gene_positions_histograms.py`
* 