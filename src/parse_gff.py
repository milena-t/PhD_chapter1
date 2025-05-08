# gff parser and associated functions for evaluating annotations
# Currently only functional for gff3

# TODO implement the stuff with the isoforms. maybe make it general so there is an infinite amount of nested features that are parent and child at the same time

import typing as tp
from dataclasses import dataclass
# dataclass generates stuff like __init__() for classes automatically 
from enum import Enum
import time
import pandas as pd
from tqdm import tqdm
import subprocess as sp
import tempfile
import os

# create an enumeration class that contains all possible categories a feature can have
# strings are exactly spelled as they are in the gff and are therefore used to categorize the lines into feature categories
# extracted with awk -F'\t' '{print $3}' a_obtectus_orthodb_isoform_filtered.gff | sort | uniq 
class FeatureCategory(str,Enum):
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
    return string_to_category.get(s)


### Class for all features that have "gene" as a parent
@dataclass # for explanation see top where dataclass library is loaded or google
class Feature:
    category:FeatureCategory
    start:int
    stop:int

    @property # property is a function without arguments, and is used like Feature.length without "()" and returns an integer
    def length(self):
        return int(self.stop)-int(self.start) +1 
        # in the gff, start and stop positions are both included in the length of the feature, like [start : stop]
        # therefore, when calculating feature length, +1 is necessary
    
    contig:str
    source:str
    strandedness:str
    
    parent_id:str
    feature_id:str

    @property
    def show(self):
        print("Feature category: ", self.category)
        print("Feature ID: ", self.feature_id, ";  Parent ID: ", self.parent_id)
        print("Feature is on contig: ", self.contig)
        print("\tin positions ", str(self.start), " to ", str(self.stop), ";  on strand: ", self.strandedness)
        print("Annotation source: ", self.source)        



### class for "gene" features that have no parent and contain a list of child-features
@dataclass
class Gene:
    category:FeatureCategory
    start:int
    stop:int
    contig:str
    source:str
    strandedness:str
    feature_id:str

    @property
    def length(self):
        return int(self.stop)-int(self.start) +1 # same reason as in Feature class
    
    # a gene has different features (in the gff file they have the gene_id as their parent_id)
    # the features_attribute contains a dictionary with the different feature categories as types, and a list of class instances of those features as respective values
    # for example: {exon : [exon_class1, exon_class2, exon_class3]}
    # and the elements can then be accessed by number in the list
    # they have the above defined feature functions associated with them, so I can say gff_dict["g1"][exon][0].length
    features:dict[FeatureCategory,list[Feature]] 

    # function to add a new "child" to the gene class
    def addFeature(self,newfeature:Feature):
        feat_cat=newfeature.category
        if feat_cat not in self.features:
            self.features[feat_cat]=[]
        self.features[feat_cat].append(newfeature)

    def filterFeature(self, transcript_id, feature_list):
        # filter features to only keep the ones that contain a transcript_id as parent_id
        # filter all features in the feature list returned from getFeatures
        for feature_category in feature_list:
            if feature_category is not FeatureCategory.Transcript:
                features_filt = [feature for feature in self.features[feature_category] if feature.parent_id == transcript_id]
                self.features[feature_category] = features_filt
            elif feature_category is FeatureCategory.Transcript:
                features_filt = [feature for feature in self.features[feature_category] if feature.feature_id == transcript_id]
                self.features[feature_category] = features_filt
        return self.features
    
    @property
    def getFeatures(self):
        features_present = [feature_category for feature_category in FeatureCategory if feature_category in self.features]
        return features_present
    
    # get number of how many instances of each feature category are pesent in a gene
    @property
    def getFeatureStats(self):
        feature_stats = {}
        for feature in self.features:
            feature_stats[feature] = {len(self.features[feature])}
        return feature_stats

    
    @property
    def show(self):
        print(self.category)
        if self.category == FeatureCategory.Gene:
            print("Gene ID: ", self.feature_id)
        elif self.category == FeatureCategory.Transcript:
            print("Transcript ID: ", self.feature_id)
        print(f"on contig: {self.contig} ;  from {str(self.start)} to {str(self.stop)} ;  on strand: {self.strandedness}")
        print("Annotation source: ", self.source)
        print("Child features present:")
        feature_stats = self.getFeatureStats
        if len(feature_stats) <1:
            print("\tthis feature has no children. If it is a Gene, check the transcript ID to access exons and introns")
        else:
            for feature in feature_stats:
                print("\t", feature, " in ", str(feature_stats[feature]), " instances")


@dataclass
class Repeat:
    start:int
    stop:int
    contig:str
    feature_id:str
    repeat_category:str

    @property
    def length(self):
        return int(self.stop)-int(self.start) +1 # same reason as in Feature class
    
    @property
    def get_category(self):
        return [self.repeat_category, self.start, self.stop]
    
    # get number of how many instances of each feature category are pesent in a gene
    @property
    def getFeatureStats(self):
        feature_stats = {}
        for feature in self.features:
            feature_stats[feature] = {len(self.features[feature])}
        return feature_stats
    
    @property
    def show(self):
        print(self.category)
        print(f"repeat category: {self.repeat_category}, repeat family: {self.repeat_family}")
        print(f"on contig: {self.contig} ;  from {str(self.start)} to {str(self.stop)} ;  on strand: {self.strandedness}")
        print(f"Annotation source: {self.source}")
        

def parse_repeats_out_dict(filepath):
    """
    get a dictionary that pairs repeat_ID and repeat class from the .out file generated by repeatmasker
    """
    temp_file = tempfile.NamedTemporaryFile(delete=False) # some lines have " *" at the end of them and it sucks -> create a cleanded temporary file that has them removed
    command = ["sed", "s/ \*$//", filepath]
    with open(temp_file.name, 'w') as temp_output:
        sp.run(command, stdout=temp_output, check=True)

    out_headers = ["SW score", "perc div.", "perc del.", "perc ins.", "query sequence", "position in query begin", "position in query end", "position in query (left)", "matching repeat", "repeat class/family", "position in repeat begin", "position in repeat end", "ID"]
    repeats_df = pd.read_csv(temp_file.name, sep="\s+", comment="#", skiprows=3, names=out_headers, dtype=str)

    os.remove(temp_file.name)

    repeats_dict = repeats_df.set_index('ID')['matching repeat'].to_dict()
    return repeats_dict


def find_gene_by_transcript(gene_dict: dict[str, Gene], transcript_id: str):
    for gene_id, gene in gene_dict.items():
        if FeatureCategory.Transcript in gene.features:
            for transcript in gene.features[FeatureCategory.Transcript]:
                if transcript.feature_id == transcript_id:
                    return gene.feature_id
    return None

def find_gene_by_RNA(gene_dict: dict[str, Gene], transcript_id: str):
    for gene_id, gene in gene_dict.items():
        if FeatureCategory.RNA in gene.features:
            for transcript in gene.features[FeatureCategory.RNA]:
                if transcript.feature_id == transcript_id:
                    return gene.feature_id
    return None


### DEPRECATED!!!
def parse_repeats_gff(filepath_out, filepath = "", verbose = True, contig_list = []):
    """
    !!! DEPRECATED !!!
    parse the gff with all assigned repeats generated by repeatmasker.
    makes a dictionary with {contig : [list_of_repeats]} with the list containing instances of the repeats class
    optionally a contig list can be provided so that only contigs that are included in the list are parsed and returned

    filepath_out is the output from repeatmasker, filepath is inferred to be .out.gff when left empty
    """

    if len(filepath) == 0:
        filepath = filepath_out+".gff"

    if verbose:
        print(f"File that is being parsed: {filepath}")
        print(f"file that is used for repeat cateogries: {filepath_out}")
        start_time = time.perf_counter()
        if len(contig_list)>0:
            print(f"a contig list has been included, only contigs that are represented in the list will be parsed")

        print(f"\nparse repeat_categories_dict")

    repeat_categories_dict = parse_repeats_out_dict(filepath_out)
    
    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"done getting repeat categories, {execution_time:.2f} seconds")
        start_time = time.perf_counter()

    gfffile={}
    with open(filepath, "r") as file:
        linelist = file.readlines()

        # determine separator (" " or "=") so that I don't have to test in every line
        tail_line = linelist[-1].split("\t")[-1].split(";")[0]
        separator = " "
        if "=" in tail_line:
            separator = "="
        
        non_repeat_feature_count = 0

        print("parsing file...")
        for line in tqdm(linelist):

            # Skip empty lines or lines starting with a comment character
            line = line.strip()
            if not line or line.startswith('#'):
                continue


            # parse gff file
            # more info on file format and columns here: https://www.ensembl.org/info/website/upload/gff.html?redirect=no
            contig,source,category_,start,stop,score,strandedness,frame,otherproperties=[c for c in line.split("\t") if len(c)>0]

            if len(contig_list)>0 and contig not in contig_list:
                continue

            category = categorize_string(category_)

            if contig not in gfffile:
                gfffile[contig] = []

            # parse all the IDs from the last gff column
            properties={}
            for attr in otherproperties.strip().split(";"):
                attr = attr.strip()
                try:
                    key,value=attr.split(separator)[-2:]
                except:
                    continue

                #if "RNA-" in attr: # Sometimes, the parent of an exon is the transcript "mRNA", which itself is a child of the "gene" feature. Since I don't have double-nested fetaures like that I will just remove that level of separation and parse the parent ID from the child feature. The [:-2] is so that the "-m" or "-t" possibilities are both removed
                #    properties[key]=value.split("RNA-")[0][:-2] 
                if "," in attr:
                    properties[key]=value.split(",")[0]
                else:
                    properties[key]=value

            if category==FeatureCategory.Repeat:
                if "ID" not in properties:
                    raise RuntimeError(f"no id property found for gene in line: {line}")

                repeat_id = properties["ID"]

                rep_category = repeat_categories_dict[repeat_id].strip().split("/")[0]


                # rep_identification= otherproperties.strip().split(" ")[1]
                # try:
                #     rep_category, rep_family = rep_identification.strip('"').split(":")[-1].split("_")
                # except:
                #     rep_category = rep_identification.split(":")[-1].strip('"')
                #     rep_family = "not_assigned" # sometimes no repeat family is assigned for motifs like Motif:(ATT)n or Motif:A-rich
                
                newrepeat=Repeat(category,int(start),int(stop),contig,source,strandedness,feature_id=properties["ID"], repeat_category=rep_category, repeat_family="")
                # print(newgene)
                gfffile[contig].append(newrepeat)
                continue # move to the next line in the gff file
            else:
                non_repeat_feature_count += 1


            # count+=1
            # print(count)
    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tparsing time: {execution_time:.2f} seconds")
        if non_repeat_feature_count>0:
            print(f"\non-repeat features: {non_repeat_feature_count}")
        print(f"\ttotal number of contigs: {len(gfffile)}")
    return(gfffile)
###

def parse_repeats_repeatmasker_outfile(filepath_out, verbose = True, filter_overlap = False, contig_list = []):
    """
    parse the .out file with all assigned repeats generated by repeatmasker.
    makes a dictionary with {contig : [list_of_repeats]} with the list containing instances of the repeats class
    optionally a contig list can be provided so that only contigs that are included in the list are parsed and returned

    this function does not filter overlaps, because many levels of repeats can overlap through insertions into insertions into an original repeat.
    Overlaps can be filtered after parsing
    """

    base_outfile = False
    # check if the file is .out or .ori.out
    if filepath_out[-8:] == ".ori.out":
        base_outfile = True

    if verbose:
        print(f"File that is being parsed: {filepath_out}")
        if base_outfile:
            print(f"the outfile is the base file")
        else:
            print(f"weirdly formatted outfile, not .ori.out")
        start_time = time.perf_counter()
        if len(contig_list)>0:
            print(f"\t(a contig list has been included, only contigs that are represented in the list will be parsed)")


    gfffile={}

    #### parse dataframe from .out file
    ## clean the outfile first in a temporary file by removing the tailing " *" that some lines have 
    if base_outfile == False:
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        command = ["sed", "s/ \*$//", filepath_out]

        if verbose:
            print(f"* Prepare input file:\n\tcleaning (may take a few seconds if there is a lot of repeats): \n\t\t{filepath_out}, creating the temporary file {temp_file.name}")

        with open(temp_file.name, 'w') as temp_output:
            sp.run(command, stdout=temp_output, check=True)

        if verbose:
            end_time = time.perf_counter()
            execution_time = end_time - start_time
            print(f"\tdone cleaning ({execution_time:.2f} seconds), reading file ...")
            start_time = time.perf_counter()

        out_headers = ["SW score", "perc div.", "perc del.", "perc ins.", "contig", "position in query begin", "position in query end", "position in query (left)", "something unidentified", "matching repeat", "repeat class/family", "position in repeat begin", "position in repeat end", "position in repeat (left)", "ID"]
        repeats_df = pd.read_csv(temp_file.name, sep="\s+", comment="#", skiprows=3, names=out_headers, dtype=str, index_col=False)
        # print(repeats_df)
        # print(repeats_df["repeat class/family"])

        os.remove(temp_file.name)
        ####
    
    else:
        out_headers = ["SW score", "perc div.", "perc del.", "perc ins.", "contig", "position in query begin", "position in query end", "position in query (left)", "something unidentified", "matching repeat", "repeat class/family", "position in repeat begin", "position in repeat end", "position in repeat (left)", "ID"]
        repeats_df = pd.read_csv(filepath_out, sep="\s+", comment="#", names=out_headers, dtype=str, index_col=False)


    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tdone reading ({execution_time:.2f} seconds), filtering out contigs shorter than the window length (may take long for fragmented assemblies) ...")
        start_time = time.perf_counter()

    repeats_df = repeats_df[["contig", "position in query begin", "position in query end", "repeat class/family", "ID"]]

    repeats_before_filter = len(repeats_df)

    repeats_list:list[list[str]] = []
    if len(contig_list)>0:
        for contig in contig_list:
            repeats_contig = repeats_df[repeats_df["contig"] == contig]
            repeats_list.extend(repeats_contig.values.tolist())
            # if verbose:
            #     print(f"{contig} : {len(repeats_contig)} repeats, (in total: {len(repeats_list)})")
    else:
        repeats_list = repeats_df.values.tolist()

    if verbose and len(contig_list)>0:
        print(f"\t\t{repeats_before_filter} repeats before filtering\n\t\t{len(repeats_list)} repeats after filtering")

    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tdone filtering ({execution_time:.2f} seconds)\n")
        start_time = time.perf_counter()
        print("* parse gff...")
    
    # check for overlaps and filter them out
    prev_end = 0
    count_overlap = 0
    bp_overlap = 0

    for line in tqdm(repeats_list):

        # parse gff file
        # more info on file format and columns here: https://www.ensembl.org/info/website/upload/gff.html?redirect=no
        contig,start_original,stop_original,repeat_category,id=line
        
        # if repeat direction is inverse, flip it
        if int(stop_original)<int(start_original):
            start = int(stop_original)
            stop = int(start_original)
        else:
            stop = int(stop_original)
            start = int(start_original)
        
        # if there's an overlap
        if int(start) < int(prev_end):
            count_overlap+=1
            if filter_overlap:
                bp_overlap = bp_overlap + int(prev_end)-int(start)
                start = int(prev_end)
        prev_end = stop

        id = str(id)

        repeat_category = repeat_category.strip().split("/")[0]

        if len(contig_list)>0 and contig not in contig_list:
            continue

        # category = categorize_string(category_)

        if contig not in gfffile:
            gfffile[contig] = []

        if len(id)==0:
            raise RuntimeError(f"no id property found for repeat in line: {line}")
        
        newrepeat=Repeat(start,stop,contig,feature_id=id, repeat_category=repeat_category)
        gfffile[contig].append(newrepeat)

    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tparsing time: {execution_time:.2f} seconds")

        print(f"\ttotal number of contigs: {len(gfffile)}")

        print(f"\tnumber of repeats that overlap with a repeat ahead of them: {count_overlap}\n\t --> {bp_overlap} bp of overlapping repeats removed during parsing of the gff file")


    return(gfffile)

###
def parse_gff3_by_contig(filepath, verbose = True, featurecategory = FeatureCategory.Gene):
    """
    parse the gff into a dictionary by contig { contig : [ FeatureCategory.Gene ] } 
    This includes by default only genes, not any of the child features, 
    the featurecategory variable is optional and can be set to pick other features than genes, for example exons
    """

    if verbose:
        print(f"File that is being parsed: {filepath}")
        print(f"feature category that was set for the output: {featurecategory}")
        start_time = time.perf_counter()

    gfffile={}
    count = 0
    with open(filepath, "r") as file:
        linelist = file.readlines()

        # determine separator (" " or "=") so that I don't have to test in every line
        # sometimes in NCBI assemblies the last few lines are comments, don't use one of these
        if "#" not in linelist[-1]:
            tail_line = linelist[-1].split("\t")[-1].split(";")[0]
        else:
            tail_line = linelist[-10].split("\t")[-1].split(";")[0]
        separator = " "
        if "=" in tail_line:
            separator = "="
        
        count_mRNA = 0
        count_trans = 0
        count_gene = 0
        count_feature = 0
        for line in tqdm(linelist):

            # Skip empty lines or lines starting with a comment character
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # parse gff file
            # more info on file format and columns here: https://www.ensembl.org/info/website/upload/gff.html?redirect=no
            contig,source,category_,start,stop,score,strandedness,frame,otherproperties=[c for c in line.split("\t") if len(c)>0]
            # categorize feature with Enum class
            #category=FeatureCategory(category_)#.split()[-1])
            category = categorize_string(category_)

            
            if category != featurecategory:
                # There's only genes supposed to be included, skip everything that isn't a gene
                continue

            if category_ == "mRNA":
                count_mRNA+=1
            elif category_ == "transcript":
                count_trans +=1
            elif category_ == "gene":
                count_gene +=1


            # parse all the IDs from the last gff column
            properties={}
            for attr in otherproperties.strip().split(";"):
                attr = attr.strip()
                try:
                    key,value=attr.split(separator)[-2:]
                except:
                    continue

                if "," in attr:
                    properties[key]=value.split(",")[0]
                else:
                    properties[key]=value
            

            # if the current feature is a gene (does not have a parent)
            if category==featurecategory: # or category==FeatureCategory.Region:
                if "ID" not in properties:
                    print(properties)
                    print(f"tail line to determine the separator: '{tail_line}'")
                    print(f"the detected separator is: '{separator}'")
                    raise RuntimeError(f"no id property found for gene in line: {line}")
            # initialize a new gene entry in the gfffile dict
                newgene=Gene(category,int(start),int(stop),contig,source,strandedness,feature_id=properties["ID"],features={})
                # print(newgene)
                if newgene.contig not in gfffile:
                    gfffile[newgene.contig]=[newgene]
                else:
                    gfffile[newgene.contig].append(newgene)
                count_feature+=1


            # count+=1
            # print(count)
    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tparsing time: {execution_time:.2f} seconds")
        print(f"\tparsed feature type: {featurecategory}")
        print(f"\tnumber of features: {count_feature}")
        print(f"\tnumber of contigs: {len(gfffile)}")
    return(gfffile)



### parse the gff3 files with transcripts and genes both as top-level features
def parse_gff3_general(filepath, verbose = True, only_genes = False):
    if verbose:
        print(f"File that is being parsed: {filepath}")
        if only_genes:
            print(f"------------------------------------------------------\n!! only_genes = True !!  ==> all non-gene features are excluded\n------------------------------------------------------")
        start_time = time.perf_counter()

    gfffile={}
    count = 0
    with open(filepath, "r") as file:
        linelist = file.readlines()

        # determine separator (" " or "=") so that I don't have to test in every line
        tail_line = linelist[-1].split("\t")[-1].split(";")[0]
        separator = " "
        if "=" in tail_line:
            separator = "="
        
        count_mRNA = 0
        count_trans = 0
        count_gene = 0
        for line in tqdm(linelist):

            # Skip empty lines or lines starting with a comment character
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            # parse gff file
            # more info on file format and columns here: https://www.ensembl.org/info/website/upload/gff.html?redirect=no
            contig,source,category_,start,stop,score,strandedness,frame,otherproperties=[c for c in line.split("\t") if len(c)>0]
            # categorize feature with Enum class
            #category=FeatureCategory(category_)#.split()[-1])
            category = categorize_string(category_)


            if only_genes and category != FeatureCategory.Gene:
                # If there's only genes supposed to be included, skip everything that isn't a gene
                continue

            if category_ == "mRNA":
                count_mRNA+=1
            elif category_ == "transcript":
                count_trans +=1
            elif category_ == "gene":
                count_gene +=1


            # parse all the IDs from the last gff column
            properties={}
            for attr in otherproperties.strip().split(";"):
                attr = attr.strip()
                try:
                    key,value=attr.split(separator)[-2:]
                except:
                    continue

                #if "RNA-" in attr: # Sometimes, the parent of an exon is the transcript "mRNA", which itself is a child of the "gene" feature. Since I don't have double-nested fetaures like that I will just remove that level of separation and parse the parent ID from the child feature. The [:-2] is so that the "-m" or "-t" possibilities are both removed
                #    properties[key]=value.split("RNA-")[0][:-2] 
                if "," in attr:
                    properties[key]=value.split(",")[0]
                else:
                    properties[key]=value
            

            # if the current feature is a gene (does not have a parent)
            if category==FeatureCategory.Gene or category==FeatureCategory.Region:
                if "ID" not in properties:
                    raise RuntimeError(f"no id property found for gene in line: {line}")
            # initialize a new gene entry in the gfffile dict
                newgene=Gene(category,start,stop,contig,source,strandedness,feature_id=properties["ID"],features={})
                # print(newgene)
                gfffile[newgene.feature_id]=newgene
                continue # move to the next line in the gff file

            if category==FeatureCategory.Repeat:
                if "ID" not in properties:
                    raise RuntimeError(f"no id property found for gene in line: {line}")
                rep_identification= otherproperties.strip().split(" ")[1]
                try:
                    rep_category, rep_family = rep_identification.strip('"').split(":")[-1].split("_")
                except:
                    rep_category = rep_identification.split(":")[-1]
                    rep_family = "not_assigned" # sometimes no repeat family is assigned for motifs like Motif:(ATT)n or Motif:A-rich
                
                newgene=Repeat(category,start,stop,contig,source,strandedness,feature_id=properties["ID"], repeat_category=rep_category, repeat_family=rep_family)
                # print(newgene)
                gfffile[newgene.feature_id]=newgene
                continue # move to the next line in the gff file

            ### if the current feature is a region then don't include in the annotation dict since we only care about features contained in "gene"
            # if category==FeatureCategory.Region:
            #     if "ID" not in properties:
            #         raise RuntimeError(f"no id property found for gene in line: {line}")
            # # initialize a new gene entry in the gfffile dict
            #     newgene=Gene(category,start,stop,contig,source,strandedness,feature_id=properties["ID"],features={})
            #     gfffile[newgene.feature_id]=newgene
            #      
            #     continue

            if "ID" not in properties:
                raise RuntimeError(f"no id property found for feature in line: {line}")

            if "Parent" not in properties:
                raise RuntimeError(f"feature is not a gene and no parentid property found for feature in line: {line}")
            
            # make a class instance of the new feature

            # if the feature is a transcript, separate it from it's gene-parent so that I can also look for features by the transcript ID
            if category==FeatureCategory.Transcript:
                if "ID" not in properties:
                    raise RuntimeError(f"no id property found for gene in line: {line}")

                # add transcript as new top-level feature
                newtranscript=Gene(category,start,stop,contig,source,strandedness,feature_id=properties["ID"],features={})
                gfffile[newtranscript.feature_id]=newtranscript

                # continue # don't continue to also add the transcript as a child feature of it's parent gene

            newfeature=Feature(category,start,stop,contig,source,strandedness,feature_id=properties["ID"],parent_id=properties["Parent"])

            parent = newfeature.parent_id
            parent_gene = None
            #if the feature is an exon or cds it is possible that the parent is a transcript 
            if parent not in gfffile: 
                parent = find_gene_by_transcript(gfffile, newfeature.parent_id)
                if parent == None:
                    # some annotations don't have transcripts, but RNA 
                    parent = find_gene_by_RNA(gfffile, newfeature.parent_id)
                    if parent == None:
                        raise RuntimeError(f"parent id {newfeature.parent_id} not found in file (at least in the lines parsed so far)")
                        #parent_gene = gfffile[parent]
                    #parent_gene = gfffile[parent]
                parent_gene = gfffile[parent]

            else:
                parent_gene=gfffile[newfeature.parent_id] # if the gene id is longer like in the native annotations
                
            # add current feature to preexisting parent
            parent_gene.addFeature(newfeature)

            # count+=1
            # print(count)
    if verbose:
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"\tparsing time: {execution_time:.2f} seconds")
        if count_mRNA>0:
            print(f"\tmRNA features: {count_mRNA}")
        if count_trans>0:
            print(f"\ttranscript features: {count_trans}")
        print(f"\tgene features: {count_gene}")
        print(f"\ttotal number of features: {len(gfffile)}")
    return(gfffile)


### filter longest isoforms
## returns a dictionary of a parsed gff file and returns one of the same format that only contains the longest transcript per gene if there are multiple present
#TODO this doesn't work yet
def filter_longest_isoforms(gff_dict):
    filt_dict = {}
    after_filtering_one_exon = 0
    for gene_id in gff_dict: # iterates only through the gene names, gff_dict[gene] accesses the actual feature information
        gene = gff_dict[gene_id] # get actual gene feature content

        # get list of all transcripts
        all_transcripts = gene.features[FeatureCategory.Transcript]

        # skip this loop iteration if there is less than two transcripts (only one isoform)
        if len(all_transcripts) < 2:
            filt_dict[gene_id] = gene
            pass
        else:
            # find longest transcript
            transcript_lengths = {}
            for transcript in all_transcripts:
                transcript_lengths[transcript.feature_id]=transcript.length
            longest_transcript_id = max(transcript_lengths, key=lambda k: transcript_lengths[k]) # google for additional info on what the key option in max does and how the lambda function works
            # use filterFeature attribute to keep only features that have longest_transcript_id as their parent
            # update gene.features directly
            gene.features = gene.filterFeature(longest_transcript_id, gene.getFeatures)
            # append to isoform filtered dictionary
            filt_dict[gene_id] = gene
            if len(gene.features[FeatureCategory.Exon])<2:
                after_filtering_one_exon +=1
            #break
    print(after_filtering_one_exon)
    return filt_dict


### get list of single-exon genes
def single_exon_genes(gff_dict):
    single_exon_IDs = []
    for gene_id in gff_dict:
        gene = gff_dict[gene_id]
        if FeatureCategory.Exon in gene.features:
            if len(gene.features[FeatureCategory.Exon])<2:
                single_exon_IDs.append(gene_id)
    return single_exon_IDs


def write_dict_to_file(dict, filename):
    with open(filename, "w") as outfile:
        for key, value_list in dict.items():
            value_str = ','.join(value_list)
            outfile.write(f"{key}:{value_str}\n")
        print("file saved in current working directory as: "+filename)


def split_at_second_occurrence(s, char = "_"): # split the gene string at the second occurence of "_" to get only the species name
    if s.count(char)<2:
        return s
    else:
        first_occurrence = s.find(char)
        second_occurrence = s.find(char, first_occurrence + 1)
        species = s[:second_occurrence]
        return species

def get_single_exon_transcripts(parsed_gff_dict):
    transcripts_single_exon = {}
    transcripts_multi_exon = {}
    for key, value in parsed_gff_dict.items():
        if value.category == FeatureCategory.Transcript:
            exons = value.features[FeatureCategory.Exon]
            if len(exons) == 1:
                transcripts_single_exon[key] = value
            else:
                transcripts_multi_exon[key] = value
        else:
            continue
    return(transcripts_single_exon, transcripts_multi_exon)

def get_transcript_lengths(parsed_gff_dict):
    transcript_lengths = {} # this is a more basic dictionary with { transcript_ID = length}
    for trans_id, attributes in parsed_gff_dict.items():
        if attributes.category == FeatureCategory.Transcript:
            exons = attributes.features[FeatureCategory.Exon]
            length = 0
            for exon in exons:
                exon_length = int(exon.stop) - int(exon.start)
                length += exon_length
            transcript_lengths[trans_id] = length
    return(transcript_lengths)


def print_single_exon_stats(filepath, include_list = True):
    try:
        gff_dict = parse_gff3_general(filepath)
    except:
        raise RuntimeError("parsing gone wrong")
    no_transcripts = {key : value for key, value in gff_dict.items() if value.category == FeatureCategory.Transcript}
    no_genes = {key : value for key, value in gff_dict.items() if value.category == FeatureCategory.Gene}
    print(f"total number of transcripts: {len(no_transcripts)}")
    print(f"total number of top-level features (genes, regions): {len(no_genes)}")

    aobt_single_exon, aobt_multi_exon = get_single_exon_transcripts(gff_dict)
    print(f"no. single exon transcripts: {len(aobt_single_exon)}")
    if include_list:
        single_exon_IDs = ",".join(aobt_single_exon)
        print(f"list of single-exon transcript IDs: {single_exon_IDs}")
    print(f"no. multi exon transcripts: {len(aobt_multi_exon)}")
    if include_list:
        multi_exon_IDs = ",".join(aobt_multi_exon)
        print(f"list of single-exon transcript IDs: {multi_exon_IDs }")

    single_exon_transcript_lengths = get_transcript_lengths(aobt_single_exon)
    mean_len = sum(single_exon_transcript_lengths.values())/len(single_exon_transcript_lengths)
    print(f"{len(single_exon_transcript_lengths)} single-exon transcripts with average length {mean_len:.6}")

    multi_exon_transcript_lengths = get_transcript_lengths(aobt_multi_exon)
    mean_len = sum(multi_exon_transcript_lengths.values())/len(multi_exon_transcript_lengths)
    print(f"{len(multi_exon_transcript_lengths)} multi-exon transcripts with average length {mean_len:.6}")


if __name__ == '__main__':

    # cmac_native_path = "/Users/milena/work/c_maculatus/C_maculatus_native_isoform_filtered.gff"
    cmac_native_path = "/Users/milena/work/c_maculatus/cmac_Lu2024_native_isoform_filtered.gff"
    # cmac_orthoDB_path = "/Users/milena/work/c_maculatus/C_maculatus_braker_isoform_filtered.gtf"
    cmac_orthoDB_path = "/Users/milena/work/c_maculatus/cmac_Lu2024_orthoDB_isoform_filtered_scaffold_name_corrected.gff"

    aobt_native_path = "/Users/milena/work/a_obtectus/a_obtectus_native_isoform_filtered.gff"
    # aobt_native_path = "/Users/milena/work/a_obtectus/a_obtectus_native_isoform_filtered.gff"
    aobt_orthoDB_path = "/Users/milena/work/a_obtectus/a_obtectus_orthodb_isoform_filtered.gff"
    # aobt_orthoDB_path = "/Users/milena/work/a_obtectus/a_obtectus_orthodb_isoform_filtered.gff"
    

    ## testing cases
    if False:
        print()
        start_time = time.perf_counter()
        aobt_native_gff = parse_gff3_general(aobt_native_path)
        aobt_native_gff["augustus_masked-@002915F_arrow_arrow-processed-gene-0.7"].show
        aobt_native_gff["augustus_masked-@002915F_arrow_arrow-processed-gene-0.8-mRNA-1"].show
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        print(f"parsing time for braker annotations: {execution_time:.2f} seconds")


    ## test collapse function for the repeats
    if True:
        Bsil_ori_out = "/Users/miltr339/work/repeatmasking/repeat_gffs/B_siliquastri_masking.ori.out"

    