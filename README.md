# Uniform Genome annotations for comparative genomics
This is all the code that was used in the first chapter of my thesis where I create uniform de-novo genome annotations to investigate gene family evolution in coleoptera. I implement and evaluate a uniform pipeline and remove as much technical bias as possible. The important aspects are:

* start with repeatmodeller and repeatmasker, even if the assembly is already masked
* exclude species-specific RNAseq (only orthoDB protein evidence)
* filter annotated protein sequences in each species with 

# Pipeline

## Repeatmasking
Use repeatmodeller and repeatmasker to re-mask all assemblies. There isn't really parameters to adjust, so just run as-is. Keep track of the repeat families faster that repeatmodeller gives because it's necessary for filtering later.

## Genome annotation with BRAKER
Use BRAKER3 with the orthoDB arthropoda dataset to uniformly annotate all genomes. Since I don't use RNAseq, technically BRAKER2 is run but I use the BRAKER3 docker container ([here](https://hub.docker.com/r/teambraker/braker3)). Since the Uppmax cluster doesn't allow docker containers/images, I used a singularity container to make a singularity image of the docker container (see `src/make_singularity_image_braker3`). I use the orthoDB Arthropoda proteinfasta provided by the BRAKER authors for evidence (v12 [here](https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/)). 

Since BRAKER can be a little difficult to run, I have written a wrapper script that does all the setup. It is extremely specific to using the singularity image on uppmax with the scratch storage, so likely not useful for anyone on a different system. for Uppsala University people, the script is here: `bash_scripts/run_braker.sh`, you pass four command line arguments when running it with sbatch: 
* species name
* assembly path
* proteinfasta path
* optional: list of SRR numbers of RNAseq runs (separated by commas with no space, directly in the command line not path to a separate file)


### Annotation post-processing
I do the standard QC ckecks and post-processing that is necessary to do gene family evolution analysis. 
1. 
