# my bash scripts
This file contains all bash scripts. Both the standalone ones and the ones that are used to run python scripts

## getting protein sets
### `run_braker.sh`
Runs braker from a singularity image (convert from the docker container that the braker people publish), and takes command line input. Internally sorts out all the uppmax requirements for braker and mounts the image to the scratch storage. After the braker run it moves the results to a directory of the same name as the provided species name. If a directory of that name already exists it will be overwritten!

### `filter_gff.sh`
After braker returns a gtf file, this file will be filtered to remove overlapping genes and only keep the longest isoform of every gene. A new gff file is created.

### `get_fasta_from_gff.sh`
After the gff is filtered, extract a fasta file of all annotated genes, one proteinfasta (.faa) and one nucleotidefasta (.fna)

### `annotation_statistics.sh`
Use the AGAT package to extract some basic annotation statistics.

### `filter_bad_proteins.sh`
calls the src/filter_bad_proteins.py script which filters out proteins that don't have correct start and stop codons. Does not work on nucleotide fasta, only proteinfasta