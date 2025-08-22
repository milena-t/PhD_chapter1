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

### Filter for repeats

I do a blast search of all the proteins extracted above against the repeat families fround by repeatmodeller to identify which protein sequences are likely repeats that were misannotated: `bash_scripts/blast_repeat_families_against_proteins.sh`

### Annotation evaluation

At this point I evaluate the annotation for several aspects, which is mainly facilitated by reading it into a data structure defined in `src/parse_gff.py` (filtering done only to the protein sequences is not taken into account! I am not back-transforming these edits onto the annotation). This uses the third column of the gff to infer the type of feature, which in braker is mainly `gene` as the top feature, then `transcript` which is a child feature to `gene`, and `exon` (and `intron`) which is a child feature to `transcript`. Other annotation pipelines use other tags, and sometimes the meaning isn't quite clear, e.g. I am unsure if `mRNA`, `rRNA`, `tRNA` and others are sub-categories of what braker calls `transcript`. I assigned all unfamiliar feature categories in the most simple way (so assuming that they are subcategories of `transcript`), but this might miss some nuances I am unfamiliar with, see the source code for how exactly I assign each feature type. This also only covers the ones I encountered in the native annotations of the species included here, there are certainly others that I did not include, but adding them would be pretty easy.

I compare my uniform annotations to native annotations for the same assembly from NCBI
* `src/plotting/compare_gene_positions_histograms.py` (both annotations are based on the assembly, I check how the positions of annotated features overlap between both)
* `src/plotting/plot_basics.py` (protein length histograms and number of transcripts)
* `src/calculate_single_exon_stats.py` (make a file with single-exon stats that can be plotted in `plot_basics.py`)
  
## Orthogroup identification with orthofinder

I run orthofinder (`bash_scripts/orthofinder.sh`) with the processed proteinfasta files and no tree.

### Filter orthogroups

### Orthofinder evaluation

Maybe add the blast thing idea i had if it goes fast

## CAFE5: identify significantly rapidly evolving orthogroups

## DAVID gene family clustering

## Genome characteristics correlations

### Genome-wide

### close to transcripts of interest