# How to create the singularity image to run braker3

## singularity and docker

Uppmax runs singularity images, not docker images, since singularity requires less permissions and is therefore better from a security perspective or something.
For the same security reasons, it also only allows you to run singularity images, not actually build singularity images.
Singularity is not available outside linux though, and the online build tool sucks. 
Therefore, I am running singularity inside a docker container (docker2sif.Dockerfile) to build a singularity image on my laptop, that I can then upload to uppmax and run.

## actually building the container

The other .sh files are referenced by docker2sif.Dockerfile internally and should stay where they are.

To run, open the docker desktop application.
The container can be built and run in the command line with:

```
docker build -t braker3singularitybuild -f /Users/miltr339/Box\ Sync/code/annotation_pipeline/make_singularity_image_braker3/docker2sif.Dockerfile .
```

Then run the container by pressing the "play" button in the desktop GUI, and specify a container name, like `make_braker3_singularityimage` 

## getting the image out of the container and onto uppmax

You can copy the image to the host with:

```
docker cp make_braker3_singularityimage:/root/braker3.sif /path/on/host/braker3.sif
```
(be careful to specify the correct container name!)

And then copy it to uppmax like this:
```
rsync -azP braker3.sif milenatr@rackham.uppmax.uu.se:/home/milenatr/private/genome_analysis_2024/results/braker_results/braker3_test/braker3.sif 
```

# How to run the container on uppmax

I need to mount all that shit correctly to the container which is fucking miserable. The below command is modified from the braker test files for the runs with the container, with a few unnecessary options removed.

```
export wd=/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_chinensis
singularity exec -B ${PWD}:${PWD} braker3.sif braker.pl --genome=/proj/naiss2023-6-65/Milena/coleoptera_sequences/c_chinensis/chinensis_from_uppmax.fasta.masked --workingdir=${wd} --GENEMARK_PATH=${ETP}/gmes --threads 1
```

$PWD is a default environmental variable of the current working directory. use one above the assembly and the proteinfiles, and then refer to them by full path. Use then working directory also below this directory
$wd is set by the user before