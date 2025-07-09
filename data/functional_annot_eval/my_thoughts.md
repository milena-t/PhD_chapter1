# Functional annotation of rapidly evolving orthogroups

The goal here is to identify interesting functional annotation orthogroups that can be summarized in a table in the main text. Probably like 8-15 or so. Ones that i definitely want to mention:

* Whatever is going on with the giant ones in A. obtectus
* Some stuff about the genes annotated with sexual reproduction
* Some stuff about the olfactory genes
* detoxification stuff about cytochrome b450
* The fluorescence stuff in i. luminosius and the other glowing one

I am doing stuff with the table in `/PhD_chapter1/src/functional_annotation_eval.py`. 

## large expansion in *A. obtectus*

A. obtectus has two large orthogroups that are twice as large as the third largest one, N0.HOG0000035 and N0.HOG0000014. The flybase IDs have no functional annotation, and are not associated with any DAVID gene group. The correlation with genome size and repeat content is significant in both cases.

# TABLE

Orthogroup_ID | CAFE_p-value | Gene_Group | Group_function | repeat_correlation_slope | repeat_correlation_p-value | GS_correlation_slope | GS_correlation_p-value | Gene_Name | D_melanogaster | I_luminosus | P_pyralis | C_septempunctata | A_verrucosus | T_castaneum | Z_morio | T_molitor | D_ponderosae | R_ferrugineus | A_obtectus | B_siliquastri | C_chinensis | C_maculatus | max_delta_GF | transcript_ID_native | Flybase | Flybase_summary
N0.HOG0000035 | 0.0 | None | None | 0.07833822365005534 | 0.017242915938343744 | 0.004367727725001441 | 0.011561278395119738 | None | 1 | 11 | 23 | 19 | 5 | 1 | 15 | 10 | 3 | 4 | 135 | 1 | 32 | 17 | 134 | rna-NM_001275459.1 | FBgn0037633 | This gene is referred to in FlyBase by the symbol Dmel\CG9839 (FBgn0037633). It is a protein_coding_gene from Dmel. It has 2 annotated transcripts and 2 polypeptides (1 unique). Gene sequence location is 3R:8803471..8805802. Its molecular function is unknown. The biological processes in which it is involved are not known. 5 alleles are reported. No phenotypic data is available. The phenotypic classes of alleles include: short lived; viable. Summary of modENCODE Temporal Expression Profile:  Temporal profile ranges from a peak of high expression to a trough of low expression.  Peak expression observed within 00-12 hour embryonic stages.
N0.HOG0000014 | 0.0 | None | None | 0.08252868448759042 | 0.029857058816618244 | 0.00544176454054848 | 0.003931094669538225 | None | 1 | 29 | 60 | 13 | 34 | 1 | 10 | 4 | 1 | 3 | 121 | 6 | 38 | 56 | 120 | rna-NM_136949.2 | FBgn0033750 | This gene is referred to in FlyBase by the symbol Dmel\CG13151 (FBgn0033750). It is a protein_coding_gene from Dmel. It has one annotated transcript and one polypeptide. Gene sequence location is 2R:12515375..12517284. Its molecular function is unknown. The biological processes in which it is involved are not known. 5 alleles are reported. No phenotypic data is available. The phenotypic class of alleles includes: viable. Summary of modENCODE Temporal Expression Profile:  Temporal profile ranges from a peak of moderate expression to a trough of low expression.  Peak expression observed within 00-18 hour embryonic stages, during late larval stages, at stages throughout the pupal period.
