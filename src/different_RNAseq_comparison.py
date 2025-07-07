"""
compare annotation statistics for the three versions of the annotation on the superscaffolded Cmac Lome assembly
"""

def filepaths_SE_stats():
    SE_stats_dict = {
        "Lome_single_exon_stats" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/Lome_single_exon_stats.txt",
        "Lu_single_exon_stats" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/Lu_single_exon_stats.txt",
        "SI_single_exon_stats" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/SI_single_exon_stats.txt",
        "Kaufmann_single_exon_stats" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/Kaufmann_nonsuperscaffoleded.txt"
    }
    return SE_stats_dict


def parse_SE_stats(filepath:str):
    """
    parse the single exon stats outfiles 
    """
    transcripts_string = "total number of transcripts"
    SE_string = "no. single exon transcripts"
    ME_String = "no. multi exon transcripts"
    SE_length_String = "single-exon transcripts with average length"
    ME_length_String = "multi-exon transcripts with average length"

if __name__ == "__main__":
    SE_stats_paths = filepaths_SE_stats()