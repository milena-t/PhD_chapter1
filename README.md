# PhD chapter1
This is all the code that was used in the first chapter of my thesis where i investigate the relationship between genome size, repeat content, and number of genes. I use de-novo repeat and gene annotation to ensure a uniform pipeline and remove any technical bias from my estimates

# Pipeline

## Repeatmasking
Use repeatmodeller and repeatmasker to re-mask all assemblies

## Genome annotation
Use BRAKER3 with the orthoDB arthropoda dataset to uniformly annotate all genomes. no RNAseq!
