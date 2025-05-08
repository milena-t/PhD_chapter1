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
    """
    def __init__(self):
        pass