# Functional annotation of rapidly evolving orthogroups

The goal here is to identify interesting functional annotation orthogroups that can be summarized in a table in the main text. Probably like 8-15 or so. Ones that i definitely want to mention:

* Whatever is going on with the giant ones in A. obtectus
* detoxification stuff about cytochrome P450
* Some stuff about the genes annotated with sexual reproduction
* Some stuff about the olfactory genes
* The fluorescence stuff in i. luminosius and the other glowing one

Other considerations:
* T. molitor and Z. morio are of particular interest because they are of industrial relevance, both as feed insects for agriculture and exotic pets, and as human food (alternative to other animal proteins). (They also use BRAKER2 and then TSEBRA manually with the default weights!! they also get high gene numbers and lots of single exon genes, especially z. morio)
    * T. molitor RNA: SRR18735292
    * z. morio RNA: SRR18735291

I am doing stuff with the table in `/PhD_chapter1/src/functional_annotation_eval.py`. 

## large expansion in *A. obtectus*

A. obtectus has two large orthogroups that are twice as large as the third largest one, N0.HOG0000035 and N0.HOG0000014. The flybase IDs have no functional annotation, and are not associated with any DAVID gene group. The correlation with genome size and repeat content is *significant* in both cases. Considering that the others are not significant in their association, maybe that's not a biological signal for selection based on function, but instead the expansion is caused by drift through TE activity (and whatever the genome size does)?

<p><img src="/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/Aobt_expansion_GF_sizes.png" alt="A. obtectus expansion" width=45%></p>

## Detoxification

Gene Group 1 is the main detoxification group, with all the cytochrome P450 in it, but Gene Group 17 which has aldehyde deoxide in it also has something to do with detoxification. All of these orthogroups are very quiet in Bruchinae and Curculionidae (D. ponderosae and R. ferrugineus). Notable expansions are especially Elateriforma (I. luminosus) which are fluorescent, and Tenebrioids.

* cytochrome P450 is detoxification and we think it could be related to host plant adaptation?

* Aldehyde oxidase is a part of the acetaldehyde metabolic process, which breaks down acetaldehyde, a toxic byprooduct of alcohol metabolization, into less toxic substances

<p>
<img src="/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/Gene_Group_1_detoxofication_GF_sizes.png" alt="Gene group 1" width=45%>
<img src="/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/Gene_Group_17_detoxofication_GF_sizes.png" alt="Gene group 17" width=45%>
</p>

### Gene Group 1: Cytochrome P450

This gene group is very insignificant for the correlation with both repeat content and genome size. All orthogroups in this gene group are: 
* N0.HOG0000027
* N0.HOG0000059
* N0.HOG0000095
* N0.HOG0000204
* N0.HOG0001077
* N0.HOG0000140
* N0.HOG0000492
* N0.HOG0001030



### Gene Group 17: aldehyde oxidase

This gene group contains only aldehyde oxidase. The only rapidly evolving orthogroup is:
* N0.HOG0000669
which is also insignificant for the correlation with repeat content and also with genome size.

## Sexual reproduction

Here, the relevant Gene Groups are these:

* Gene Group 5, Function: protease inhibitor (immunity, reproduction)
* Gene Group 8, Function: immunity and sexual reproduction
* Gene Group 9, Function: sexual reproduction

Which are these orthogroups (all are unsignificant in GS and TE correlations):
* N0.HOG0000541 (Gene Group 5,Gene Group 8,Gene Group 25)
* N0.HOG0000775 (Gene Group 8)
* N0.HOG0000892 (Gene Group 8)
* N0.HOG0000401 (Gene Group 9) 
* N0.HOG0009002 (Gene Group 9)
  
similar to the detoxification, which means not a lot going on in Bruchinae and Curculionidae, but there's one massively expanding in Z. morio (N0.HOG0000401), whose expression peak is in adult males. TODO read about that one specifically? what is it about the reproduction in that species?

<p>
<img src="/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/Gene_Group_5_8_9_reproduction_GF_sizes.png" alt="Gene group 5,8,9" width=45%>
<img src="/Users/milena/work/PhD_chapter1_code/PhD_chapter1/data/functional_annot_eval/OG_N0.HOG0000401_reproduction_GF_sizes.png" alt="N0.HOG0000401" width=45%>
</p>



# Papers
* tenebrioid assembly and gene family evolution: https://academic.oup.com/g3journal/article/13/6/jkad079/7099445
* beetles are super underrepresented in sequencing projects: https://doi.org/10.1073/pnas.2109019118


# TABLE

| Orthogroup_ID | CAFE_p-value | Gene_Group | Group_function | repeat_correlation_slope | repeat_correlation_p-value | GS_correlation_slope | GS_correlation_p-value | Gene_Name | D_melanogaster | I_luminosus | P_pyralis | C_septempunctata | A_verrucosus | T_castaneum | Z_morio | T_molitor | D_ponderosae | R_ferrugineus | A_obtectus | B_siliquastri | C_chinensis | C_maculatus | max_delta_GF | transcript_ID_native | Flybase |
| N0.HOG0000035 | 0.0 | None | None | 0.07833822365005534 | 0.017242915938343744 | 0.004367727725001441 | 0.011561278395119738 | None | 1 | 11 | 23 | 19 | 5 | 1 | 15 | 10 | 3 | 4 | 135 | 1 | 32 | 17 | 134 | rna-NM_001275459.1 | FBgn0037633 |
| N0.HOG0000014 | 0.0 | None | None | 0.08252868448759042 | 0.029857058816618244 | 0.00544176454054848 | 0.003931094669538225 | None | 1 | 29 | 60 | 13 | 34 | 1 | 10 | 4 | 1 | 3 | 121 | 6 | 38 | 56 | 120 | rna-NM_136949.2 | FBgn0033750 |
| N0.HOG0000140 | 0.0 | Gene Group 1 | detoxification | 0.0074300433129667505 | 0.7286375887369507 | 0.00030145542591610025 | 0.7920270858115639 | Cytochrome P450 6w1(Cyp6w1) | 6 | 26 | 8 | 8 | 9 | 8 | 14 | 17 | 1 | 6 | 5 | 3 | 2 | 6 | 25 | rna-NM_001299216.1 | FBgn0033065 |
| N0.HOG0000401 | 0.0 | Gene Group 9 | sexual reproduction | -0.022203390866792034 | 0.3806646826028459 | -0.0013093512895652032 | 0.33042982105018115 | uncharacterized protein(CG34189) | 3 | 10 | 6 | 4 | 9 | 3 | 21 | 5 | 1 | 3 | 1 | 1 | 1 | 1 | 20 | rna-NM_001103877.2 | FBgn0085218 | 
