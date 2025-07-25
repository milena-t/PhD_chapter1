# Functional annotation of rapidly evolving orthogroups


The goal here is to identify interesting functional annotation orthogroups that can be summarized in a table in the main text. Probably like 8-15 or so. Ones that i definitely want to mention:

* Whatever is going on with the giant ones in A. obtectus
    * Transposases misannotated in flybase
* detoxification stuff about cytochrome P450
    * some interesting stuff potentially about polyethylene breakdown in tenebrioids?
* Some stuff about the genes annotated with sexual reproduction
    * Nothing concrete
* Some stuff about the olfactory genes
    * the chitin stuff is development, not odorant binding
* The fluorescence stuff in i. luminosius and the other glowing one
    * esterase and luciferase
* Growth factor and development
    * the chitin stuff that is not pheromone sensing

Other considerations:
* T. molitor and Z. morio are of particular interest because they are of industrial relevance, both as feed insects for agriculture and exotic pets, and as human food (alternative to other animal proteins). (They also use BRAKER2 and then TSEBRA manually with the default weights!! they also get high gene numbers and lots of single exon genes, especially z. morio)
    * T. molitor RNA: SRR18735292
    * z. morio RNA: SRR18735291

I am doing stuff with the table in `/PhD_chapter1/src/functional_annotation_eval.py`. 

## CAFE significant orthogroups

CAFE is a maximum likelihood method and needs to be run multiple times to ensure convergence. I have run CAFE 20 times, and the lambda values all converge nicely to ca. 0.4067, and each run has between 530 and 570 significant orthogorups. I want to inlcude only orthogroups where I can be certain they are significant, so I decided to only include ones significant in all runs, which are 496 orthogroups (out of 8315 considered.)

<p><img src="../CAFE_convergence/runs_sig_OGs_hist.svg" alt="Sig OG abundance histogram" width=30%></p>

Other pre-processing I did to the CAFE input:

1. only include orthogroups that have at least one member in *D. melanogaster* because this means they are for sure only present at the root of the tree.
2. only inlcude orthogroups with fewer than 100 gene family members in any species. This excludes (after the first step) only two more, which have huge expansions in *A. obtectus* and are explained below.


## large expansion in *A. obtectus*

A. obtectus has two large gene families, N0.HOG0000035 and N0.HOG0000014, which are twice as large as the third largest one. The flybase IDs of these orthogroups have no functional annotation, and are not associated with any DAVID Gene Family Cluster. The correlation with genome size and repeat content is *significant* in both cases. Considering that most of the other orthogroups investigated here are not significant in their association, maybe that means this is not a biological signal for selection based on function, but instead the expansion is caused by drift through TE activity (and whatever the genome size does)?

I used the [ncbi conserved domain search](https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi) on all the proteins in these orthogroups, and for both of them most proteins match transposases, and sometimes at the c-terminus zinc-finger related domains (which is also common for transposases), which probably means that this is a false positive gene. I also then checked the uniform *D. melanogaster* orthologs and the corresponding flybase IDs (FBgn0037633 for N0.HOG0000035 and FBgn0033750 for N0.HOG0000014) and all of those also follow the pattern of a very good transposase hit and a less good zinc-finger hit at the C-terminus. Therefore this is likely a false positive.


<p><img src="Aobt_expansion_GF_sizes.svg" alt="A. obtectus expansion" width=45%></p>


## Detoxification

Gene Family Cluster 1 is the main detoxification group, with all the cytochrome P450 in it, but Gene Family Cluster 17 which has aldehyde deoxide in it also has something to do with detoxification. All of these orthogroups are very quiet in *Bruchinae* and *Curculionidae* (*D. ponderosae* and *R. ferrugineus*). Notable expansions are especially *Elateriforma* (*I. luminosus*) which are fluorescent, and *Tenebrioids*. The two functional anotations in this category are:

* cytochrome P450 (CYP) is detoxification and we think it could be related to host plant adaptation?
* Aldehyde oxidase is a part of the acetaldehyde metabolic process, which breaks down acetaldehyde into less toxic substances

*Tenebrioid* insects such as *T. molitor* and *Z. morio* can break down plasctics in their gut as larvae. This is possible through specific microorganisms in their gut working together with insect-encoded proteins that are secreted (the degradation rate of isolated microorganisms is significantly lower). [Source](https://link.springer.com/article/10.1007/s10924-023-03029-z#Sec11) for all that follows. The first step is polyethylene being oxidated by Cytochrome p450 (and others), then more chemistry (involving aldehyde dehydrogenase, which is not the same as aldehyde oxidase!), until it is broken down by lipid metabolism in the cells. According to the flybase annotations, the expression profiles of the Gene Family Cluster 1 orthogroups have their peak mostly during larval, and pupae stages, but not exclusively (some are adult sex specific). The differential expression analysis does not provide super strong evidence for CYP though.


### Gene Family Cluster 1: Cytochrome P450

All Orthogroups in this Gene Family Cluster are very insignificant for the correlation with both repeat content and genome size. I am wondering if this is a good thing, since the other functional stuff below is not significantly correlated, and e.g. the false-positive annotations for the large orthogroups in *A. obtectus* are.

This proliferation is unsurprising (expected even?), especially in insects ([source](https://link.springer.com/article/10.1186/1471-2164-14-174), *T. castaneum* already has more duplications in the CYP superfamily than *B. mori* or *D. melanogaster*). Also detoxification, yes, but they also talk about it in relation to larval development and egg maturation. Some subfamilies (existing but not found as significant here) have also been linked to insecticide resistance in *T. castaneum*. Cytochrome P450 according to [this paper](https://link.springer.com/article/10.1186/1471-2164-14-174#Sec2) is also important in mitochondria, but Gene Family Cluster 33, which is mitochondrial translation shows no interesting dynamics.


### Gene Family Cluster 3: lipid metabolic process

Since the lipid metabolism in the cells is the end of the polyethylene breakdown, this Gene Family Cluster is also relevant, because it's the last step in the digestion process in the cells. N0.HOG0000085 shows an expansion in *Z. morio*, which also relates to the polyethylene synthesis described above. The other two are pretty boring


### Gene Family Cluster 18: aldehyde oxidase

This Gene Family Cluster contains only aldehyde oxidase. The only rapidly evolving orthogroup is N0.HOG0000669, which is insignificant for the correlation with repeat content and also with genome size. It is expressed during all developmental stages and adults of both sexes. There is what seems like a significant expansion in *I. luminosus*, but *I. luminosus* has a shit assembly so take that with a grain of salt.

### Plot

<p>
<img src="detoxificatoin_clusters.svg" alt="Gene Family Cluster 17" width=45%>
</p>


Aldehyde is a byproduct of alcohol metabolism (in humans?) and is broken down into less toxic molecules by aldehyde oxidase. In insects, aldehyde can come from other sources, such as pheromones, and "*Long-chain, unsaturated alcohol and aldehyde compounds are common  female-produced  sex  pheromones (REF)*" ([source](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0067794)). They conclude that "*\[The\] results suggest that an antennae-specific aldehyde oxidasefrom  the  navel  orangeworm,  AtraAOX2,  might  be  involved  indegradation  of  host  plant  volatile  compounds  and  pheromone*", which helps both with distinguishing plant compounds from the pheromones they are actually interested in, and also protects them from volatile plant compounts and potentially even pesticides.

Quote from the abstract: "*Our data suggest AtraAOX2 might be involved in degradation of a diversity of aldehydes including sex pheromones, plant-derived semiochemicals, and chemical cues for oviposition sites. Additionally, AtraAOX2 could protect the insect's olfactory system from xenobiotics, including pesticides that might reach the sensillar lymph surrounding the olfactory receptor neurons.*" (I don't understand the methods or biochemistry in general enough to judge how accurate this is). 


## Sexual reproduction 

Here, the relevant Gene Family Clusters are these, and all orthogroups are insiginificant for correlations with GS or repeat content. It seems like immunity and sexual reproduction are tightly linked.

* Gene Family Cluster 5, Function: protease inhibitor (immunity, reproduction)
* Gene Family Cluster 8, Function: immunity and sexual reproduction
* Gene Family Cluster 9, Function: sexual reproduction
  
Similar to the detoxification, which means not a lot going on in Bruchinae and Curculionidae, but there's one massively expanding in *Z. morio* (N0.HOG0000401), whose expression peak is in adult males, and is predicted to be in the ER according to flybase, the molecular function is just protein homodimerization. Otherwise there is not a lot of info, even in *D. melanogaster*. I have a hard time finding good papers for this based on this extremely vague information.  The other orthogroups are also not super helpful, none show very interesting dynamics that I could read into.

<p>
<img src="sexual_reproduction.svg" alt="Gene Family Cluster 5,8,9" width=45%>
</p>


## Pheromone sensing

This also involves cuticular proteins due to how insects actually do the sensing, so here I think these are the relevant Gene Family Clusters:

* Gene Family Cluster 7, Function: odorant binding
* Gene Family Cluster 20, Function: transmembrane transport (olfactory) 
* Gene Family Cluster 30, Function: pheromone sensing

Chitin can also be related to pheromone sensing due to cuticular hydrocarbons. No orthogroup that is significant in any CAFE run is annotated with cuticular hydrocarbons, and the others with chitin/cuticular functions are expresssed only in early development and regulate chitin formation.

### Gene Family Cluster 7, 20 and 30 

None of them have a significant correlation with GS or TE content except N0.HOG0001445 in cluster 30.
Cluster 20, ATPase-coupled transmembrane transporter mostly is not annotated on flybase (N0.HOG0000761 is annotated with "response to toxic substance"). According to [this paper](https://www.mdpi.com/2075-4450/15/12/1016), it's expressed in antennae in drosophila and might be related to odorant processing.

<p>
<img src="pheromone_sensing_clusters.svg" alt="Pheromone sensing" width=45%>
</p>

The highlight is N0.HOG0000037 for Gene Family Cluster 30 and N0.HOG0000056 for Gene Family Cluster 7, which make the "M" shape with the peaks at *A. verrucosus* and *Z. morio* or *T. molitor* respectively. Both of them have their expression peak in adult males.

N0.HOG0000038 (GG30) and N0.HOG0000436 (GG7) have the high peak in *I. luminosus*.

## Development

### Gene Family Cluster 15, 24, 

Chitin related and cuticular protein. No significant correlation with repeats or GS not even in the large cluster 26, except N0.HOG0000108 in Cluster 24

There is again expansions in *Z. morio*, for both categories. for GC16, it's N0.HOG0001194 and N0.HOG0003035, the former (higher peak) is only expressed in early embryonic development, but the latter is an intercellular matrix component also expressed in adults, so it might have something to do with odorants. All of Gene Family Cluster 23 is stuff only expressed during larval development and is responsible for cuticle development, nothing to do with odorant receptors. For this reason I decided to investigate other Gene Family Clusters that also have function in early development, such as Gene Family cluster 26 and 13. TODO

Also this time there's stuff going on in bruchids! GG23:N0.HOG0000108 (cluster 24) and N0.HOG0000307 (cluster 15) have expansions in bruchids. However, both of these are involved in cuticular development during larval and pupae stages, where their expression peak also happens, so this is not related to adult pheromone sensing.


### Gene Family Cluster 11 and 26

Gene Family cluster 26 is annnotated as the recently identified Adenosine deaminase-related growth factor in drosophila ([paper](https://www.sciencedirect.com/science/article/abs/pii/S0378111901007624)). Maybe they also have something to do with the early development? It is always expressed at low levels, but peak expression is in adult males. Gene Cluster 26 is early development in general. All are insignificant for the correlation with GS or repeat content.

### Plots

<p>
<img src="early_development.svg" alt="early development" width=45%>
<img src="early_development_GF_cluster_26.svg" alt="cluster 26" width=45%>
</p>

## Fluorescence in *Elateriformia*

### Gene Family Cluster 4

Gene Family Cluster 4 (Esterase and mating behavior) expanded in *Elateriforma*, and is also common in the significantly rapidly expanding gene families, which is not the case in other species. It also expands in *Z. morio*, but there only 3 of the expanding gene families overall are in this group, and this species has a lot of expansions in general


### Acyl-CoA synthetase

According to [this paper](https://web.archive.org/web/20180722105054id_/https://www.biorxiv.org/content/biorxiv/early/2017/12/21/237586.full.pdf) (especially figure 3), the *Elateriforma* share a common genomic mechanism for their fluorescence made up of luciferases and their paralogs, which form a cluster of genes (that's what they call it, no clue what that means exactly. I assume a gene family since they talk about how it came to be through tandem duplication later) containing the luciferase and close relatives such as peroxisomal fatty acyl-CoA synthetase (PACS) and non-peroxisomal acyl-CoA synthetase (ACS). This is not represented in any Gene Family Cluster, but the gene family N0.HOG0000613 is annotated as *Acyl-CoA synthetase family member 2* (but there's no other family members annotated). Some other orthogroups have Acyl-CoA synthesis-related annotations (not oxidase/reductase/hydrolase!) in their API summaries:

* N0.HOG0000284 (maybe a good candidate for the table?)
* N0.HOG0000397 (pudgy in Dmel)

### Plots

Seems like these all show expansions in *Elateriforma*, cool!

<p>
<img src="fluorescence_elateriformia.svg" alt="elateriformia" width=45%>
</p>

# General "enrichment" of Gene Family Clusters in rapidly expanding gene families

I looked at all the gene families in a species, and selected the ones that are in the upper 5th size percentile of all gene families (in species with few genes, this is sometimes just with one or more members, so I also set a min. size limit of 2 members here). Then I checked which ones are significantly rapidly evolving according to CAFE, and of these which ones have been assigned to functional Gene Family Clusters by their *D. melanogaster* member.  This is only a very small subset of all gene families in the upper 5th size percentile. I have then compiled a file that shows how many of the above selected gene families fall under each functional Gene Family Cluster. I think this can get a more general idea of the kind of things that might be evolving in some of the species. The full file is here: `PhD_chapter1/data/functional_annot_eval/frequent_Gene_Groups_in_expanding_GFs.txt` and I am going to show the output for all species below, only including Gene Family Clusters that appear at least four times and are not uncharacterized or unannotated.

These are my conclusions: 

* Gene Family Cluster 25 (glycolysis and early development), and Gene Family Cluster 1 (detoxification, cytochrome P450) are basically in all of them.
* Gene Family Cluster 2 (proteolysis) is in all families except bruchids (where it only occurs in B. siliquastri)
* Tenebrioids all have Gene Family Cluster 30 (pheromone sensing)
* Elateriformia all have Gene Family Cluster 4: Esterase and mating behavior.


## Shortened output of functional summary

#### D_melanogaster

number of gene families with more than 2 members (upper 5th percentile) = 506 (of 8760 orthogroups). 25 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 11 GFs annotated as : glycolysis and early development (and other)
- 10 GFs annotated as : protein breakdown (proteolysis)
- 6 GFs annotated as : detoxification

### Elateriformia

#### I_luminosus

number of gene families with more than 3 members (upper 5th percentile) = 1276 (of 13090 orthogroups). 27 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 15 GFs annotated as : protein breakdown (proteolysis)
- 11 GFs annotated as : glycolysis and early development (and other)
- 6 GFs annotated as : detoxification
- 4 GFs annotated as : Esterase and mating behavior

#### P_pyralis

number of gene families with more than 3 members (upper 5th percentile) = 1213 (of 12644 orthogroups). 27 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 20 GFs annotated as : glycolysis and early development (and other)
- 11 GFs annotated as : protein breakdown (proteolysis)
- 6 GFs annotated as : detoxification
- 5 GFs annotated as : transcription regulation
- 4 GFs annotated as : odorant binding
- 4 GFs annotated as : Esterase and mating behavior

### Coccinellidae

#### C_septempunctata

number of gene families with more than 2 members (upper 5th percentile) = 968 (of 10892 orthogroups). 27 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 16 GFs annotated as : glycolysis and early development (and other)
- 5 GFs annotated as : protein breakdown (proteolysis)
- 4 GFs annotated as : detoxification


### Tenebrionidae

#### A_verrucosus

number of gene families with more than 2 members (upper 5th percentile) = 914 (of 12405 orthogroups). 27 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 28 GFs annotated as : glycolysis and early development (and other)
- 10 GFs annotated as : protein breakdown (proteolysis)
- 6 GFs annotated as : detoxification
- 6 GFs annotated as : pheromone sensing

#### T_castaneum

number of gene families with more than 2 members (upper 5th percentile) = 548 (of 11778 orthogroups). 31 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 27 GFs annotated as : glycolysis and early development (and other)
- 8 GFs annotated as : protein breakdown (proteolysis)
- 7 GFs annotated as : detoxification
- 6 GFs annotated as : pheromone sensing

#### Z_morio
number of gene families with more than 2 members (upper 5th percentile) = 1373 (of 13008 orthogroups). 34 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 28 GFs annotated as : glycolysis and early development (and other)
- 12 GFs annotated as : protein breakdown (proteolysis)
- 6 GFs annotated as : pheromone sensing
- 5 GFs annotated as : detoxification
- 4 GFs annotated as : chromatin organization and transcription regulation
- 4 GFs annotated as : neurological (mostly uncharacterized)

#### T_molitor
number of gene families with more than 2 members (upper 5th percentile) = 732 (of 11174 orthogroups). 34 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 19 GFs annotated as : glycolysis and early development (and other)
- 8 GFs annotated as : protein breakdown (proteolysis)
- 7 GFs annotated as : detoxification
- 6 GFs annotated as : pheromone sensing


### Curculionidae

#### D_ponderosae

number of gene families with more than 2 members (upper 5th percentile) = 839 (of 10658 orthogroups): 23 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 24 GFs annotated as : glycolysis and early development (and other)
- 10 GFs annotated as : protein breakdown (proteolysis)
- 5 GFs annotated as : detoxification

#### R_ferrugineus

number of gene families with more than 2 members (upper 5th percentile) = 600 (of 11320 orthogroups). 23 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 19 GFs annotated as : glycolysis and early development (and other)
- 7 GFs annotated as : protein breakdown (proteolysis)
- 6 GFs annotated as : detoxification


### Bruchinae

#### A_obtectus

number of gene families with more than 3 members (upper 5th percentile) = 1211 (of 12608 orthogroups). 20 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 16 GFs annotated as : glycolysis and early development (and other)
- 12 GFs annotated as : transcription regulation
- 7 GFs annotated as : detoxification
- 5 GFs annotated as : chromatin organization and transcription regulation

#### B_siliquastri

number of gene families with more than 2 members (upper 5th percentile) = 433 (of 11255 orthogroups). 18 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 13 GFs annotated as : glycolysis and early development (and other)
- 5 GFs annotated as : detoxification
- 4 GFs annotated as : protein breakdown (proteolysis)

#### C_chinensis

number of gene families with more than 3 members (upper 5th percentile) = 877 (of 12655 orthogroups). 25 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 18 GFs annotated as : glycolysis and early development (and other)
- 6 GFs annotated as : chromatin organization and transcription regulation
- 5 GFs annotated as : transcription regulation
- 4 GFs annotated as : detoxification

#### C_maculatus

number of gene families with more than 4 members (upper 5th percentile) = 1106 (of 14402 orthogroups). 22 unique Gene Family Clusters with functional annotations, these ones appear at least four times:
- 9 GFs annotated as : glycolysis and early development (and other)
- 5 GFs annotated as : detoxification




# Papers 

These are papers I didn't cite anywhere above but might still want to reference in my manuscript.
* tenebrioid assembly and gene family evolution: https://academic.oup.com/g3journal/article/13/6/jkad079/7099445
* beetles are super underrepresented in sequencing projects: https://doi.org/10.1073/pnas.2109019118


# TABLE

These are the orthogroups that should be included in a table in the paper

* N0.HOG0000035 
* N0.HOG0000014 
* N0.HOG0000140 
* N0.HOG0000085
* N0.HOG0000401 
* N0.HOG0000037
* N0.HOG0000056
* N0.HOG0000284