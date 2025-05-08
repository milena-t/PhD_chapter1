from enum import Enum
import time
import pandas as pd
from tqdm import tqdm
import subprocess as sp
import tempfile
import os


class FeatureCategory(str,Enum):
    """
    All possible feature categories that can be present in an annotation 
    and example strings of how they could be represented in the gff file
    """
    Gene = "gene"
    CDS = "CDS"
    Exon = "exon"
    Transcript = "transcript"
    Intron = "intron"
    Start_codon = "start_codon"
    Stop_codon = "stop_codon"
    TP_UTR = "three_prime_UTR"
    FP_UTR = "five_prime_UTR"
    RNA = "mRNA" # this and below are specific to the native annotations only
    Region = "region"
    Repeat = "dispersed_repeat"


string_to_category = { # define aliases when some strings should map to the same feature category
    "gene": FeatureCategory.Gene,
    "pseudogene": FeatureCategory.Gene,
    "sequence_feature": FeatureCategory.Gene,
    "mobile_genetic_element": FeatureCategory.Gene,
    "CDS": FeatureCategory.CDS,
    "cds": FeatureCategory.CDS, 
    "exon": FeatureCategory.Exon,
    "transcript": FeatureCategory.Transcript,
    "primary_transcript": FeatureCategory.Transcript,
    "intron": FeatureCategory.Intron,
    "start_codon": FeatureCategory.Start_codon,
    "stop_codon": FeatureCategory.Stop_codon,
    "three_prime_UTR": FeatureCategory.TP_UTR,
    "five_prime_UTR": FeatureCategory.FP_UTR,
    "mRNA": FeatureCategory.Transcript, # I will also use mRNA as "transcript" because even though braker calls transcripts "transcript", 
                                        # other annotation pipelines often say "mRNA" instead when they mean the same thing. 
                                        # Correctly identifying the transcirpt-section of a feature is important for working with e.g. orthofinder output
                                        # which works with transcript IDs only, so if I want to also get transcript stats for native and breaker annotations, 
                                        # I need to group it with the "transcript" feature category instead of the "RNA" feature category.
    "RNA": FeatureCategory.RNA,
    "antisense_RNA": FeatureCategory.RNA,
    "miRNA": FeatureCategory.RNA,
    "piRNA": FeatureCategory.RNA,
    "rRNA": FeatureCategory.RNA,
    "snoRNA": FeatureCategory.RNA,
    "snRNA": FeatureCategory.RNA,
    "ncRNA": FeatureCategory.RNA,
    "tRNA": FeatureCategory.RNA,
    "lnc_RNA": FeatureCategory.RNA,
    "RNase_MRP_RNA": FeatureCategory.RNA,
    "RNase_P_RNA": FeatureCategory.RNA,
    "SRP_RNA": FeatureCategory.RNA,
    "dispersed_repeat": FeatureCategory.Repeat,
    "region": FeatureCategory.Region, # since it has no parent_ID
    "repeat_region": FeatureCategory.Region, # since it has no parent_ID
    
}

def categorize_string(s):
    """
    sort the category string into one of the existing feature categories
    """
    return string_to_category.get(s)


class Feature:
    """
    any feature that could be present in the gff file
    the gff columns are:
    contig,source,category_,start,stop,score,strandedness,frame,otherproperties
    """
    def __init__(self, feature_id:str, contig:str, category:FeatureCategory, start:int, end:int, strandedness:str, frame:str, parent_id=None, child_id=None):
        
        self.feature_id = feature_id

        self.contig = contig
        self.category = category 
        self.start = start
        self.end = end
        self.strandedness = strandedness
        self.frame = frame

        # since this is general for all features, parent ID and child ID are default None, since they're not always present
        self.parent_id = parent_id 
        self.child_id = child_id

    def __repr__(self):
        return "Feature"
    
    def __str__(self):
        return(
            f"""
            Feature ID: {self.feature_id}, Feature category: {self.category}
            \tParent ID: {self.parent_id}, 
            \tchild ID: {self.child_id}
            \tOn contig: {self.contig}, from {self.start} to {self.end} (on strand {self.strandedness})
            """
        )

if __name__ == "__main__":
    print(Feature("1", "contig1", FeatureCategory.Gene, 1, 100, "+", ".", parent_id=None, child_id= ["child1", "child2"]))