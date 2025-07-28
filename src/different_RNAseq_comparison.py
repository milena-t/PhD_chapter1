"""
compare annotation statistics for the three versions of the annotation on the superscaffolded Cmac Lome assembly
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FuncFormatter
import numpy as np


def filepaths_SE_stats():
    SE_stats_dict = {
        "Lome_RNA" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/Lome_single_exon_stats.txt",
        "Nigeria_RNA" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/Lu_single_exon_stats.txt",
        "SI_RNA" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/SI_single_exon_stats.txt",
        "Kaufmann_native" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/Kaufmann_nonsuperscaffoleded.txt",
        "Lu_native" : "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/Lu_native.txt",
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

    out_dict = {}
    with open(filepath, "r") as infile:
        lines = infile.readlines()
        for line in lines:
            if transcripts_string in line:
                line_num = line.strip().split(": ")[-1]
                out_dict[transcripts_string] = int(line_num)
            if SE_string in line:
                line_num = line.strip().split(": ")[-1]
                out_dict[SE_string] = int(line_num)
            if ME_String in line:
                line_num = line.strip().split(": ")[-1]
                out_dict[ME_String] = int(line_num)
            if SE_length_String in line:
                line_num = line.strip().split(": ")[-1]
                out_dict[SE_length_String] = float(line_num)
            if ME_length_String in line:
                line_num = line.strip().split(": ")[-1]
                out_dict[ME_length_String] = float(line_num)
        
    return out_dict



def plot_single_exon_no_species_specific_three_annot(numbers, filename = ""):

    print(f" plotting for these {len(numbers)} species: \n{numbers}")

    # X coordinates for the groups
    x = np.arange(len(numbers))
    annotation_names = list(numbers.keys())

    # figure proportions according to the data included (longer or shorter)

    # fontsize scales with the dpi somehow which i have to do extra because i change the aspect ratio manually below
    fs = 45 # 37 originally
    
    
    width = 0.2 # (this is a fraction of the standardized 1 unit of space between axis ticks)
    aspect_ratio = 20 / 12 # nice for a presentation


    height_pixels = 2000  # Height in pixels
    width_pixels = int(height_pixels * aspect_ratio)  # Width in pixels
    fig = plt.figure(figsize=(width_pixels / 100, height_pixels / 100), dpi=100)
    ax = fig.add_subplot(111)

    colors = {
        "RNAseq" : "#9C4C32", # "#21295C",
        "uniform" : "#F2933A",
        "native" : "#b82946",
        "SE_length" : "#94ADC7", # lighter blue
        "ME_length" : "#6C8EB2" # darker blue
    }

    hatch_color = '#ffffff' # '#E2D4CA' #kind of eggshell white
    plt.rcParams['hatch.color'] = hatch_color


    #### plot native annotation ####

    transcript_numbers = [numbers[annot]['total number of transcripts'] for annot in annotation_names]
    SE_numbers = [numbers[annot]['no. single exon transcripts'] for annot in annotation_names]
    ME_numbers = [numbers[annot]['no. multi exon transcripts'] for annot in annotation_names]

    SE_transcript_len = [numbers[annot]['single-exon transcripts with average length'] for annot in annotation_names]
    ME_transcript_len = [numbers[annot]['multi-exon transcripts with average length'] for annot in annotation_names]
    y_max_genes = max(ME_transcript_len)

    color = [colors["native"] if "native" in annot else  colors["RNAseq"] for annot in annotation_names]
    hatching = {
        "yes" : "//", 
        "no" : "" }
    
    print(f"bar width = {width}")

    # total number of genes (with single-exons hatched)
    annotation_base = ax.bar(x - width, SE_numbers, width = width, label='proportion of which are single-exon', color= color, hatch=hatching["yes"])
    annotation_top = ax.bar(x - width, ME_numbers, width = width, bottom=SE_numbers, label='all genes in uniform RNAseq annotation', color=color, hatch=hatching["no"])

    # Make some labels.
    labels = [f"{transcript_number}" for transcript_number in transcript_numbers]

    for i, label in enumerate(labels):
        height = int(label)
        ax.text(
            x[i]-width, height + 5, label, ha="center", va="bottom", fontsize=fs*0.8, color= color[i]
        )

    # single and multi exon average length
    ax2 = ax.twinx()
    ax2.set_ylim(bottom=0, top=y_max_genes*1.3)
    ax2.set_ylabel('average transcript length', color = colors["ME_length"], fontsize = fs)
    ax2.tick_params(axis ='y', labelcolor = colors["ME_length"], labelsize = fs) 
    ax2.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: '' if x < 1 or x>max(ME_transcript_len)*1.1 else f'{x / 1e3:.0f} kbp'))

    narrow=0.5 # make the bars a little more narrow
    SE_transcript_length = ax2.bar(x, SE_transcript_len, width= width*narrow, label='avg. single-exon cds length', color=colors["SE_length"], hatch=hatching["no"])
    ME_transcript_length = ax2.bar(x+width*narrow, ME_transcript_len, width= width*narrow, label='avg. multi-exon cds length', color=colors["ME_length"], hatch=hatching["no"])

    #### set up labels and stuff ####
    
    ax.set_ylabel('Number of genes', fontsize=fs+4)
    ax.set_title('Gene numbers and proportion of single-exon genes of \nsuperscaffolded assembly from the Lome population', fontsize=fs*1.3)
    ax.set_xticks(x)

    ax.set_xlabel('', fontsize=fs+4)
    xtick_labels = [species.replace("_", " ") for species in annotation_names]
    x_rotation=0
    ax.set_xticklabels(xtick_labels, rotation=x_rotation, fontsize=fs)
    ax.set_yticklabels([f'{int(tick)/1e3:.0f}k' for tick in ax.get_yticks()], fontsize=fs)

    # make custom legend patch for the dashed bars
    plt.rcParams.update({'hatch.color': "#3f3832ff"})
    dashed_handle = mpatches.Patch(hatch = "//", alpha = 0.0)
    native_handle = mpatches.Patch(color=colors["native"])
    dashed_label = "proportion of genes with no introns"
    native_label = "all genes in native annotation"

    # Legend with custom order
    handles, labels = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    handles = handles[1:]
    labels = labels[1:]
    handles.append(native_handle)
    labels.append(native_label)
    handles.append(dashed_handle)
    labels.append(dashed_label)

    handles = handles+handles2
    labels = labels+labels2

    # new_order = [1,3,5]
    # handles = [handles[idx] for idx in new_order]
    # labels = [labels[idx] for idx in new_order]
    # new_order = [0,2,1,3,4]
    # handles = [handles[idx] for idx in new_order]
    # labels = [labels[idx] for idx in new_order]

    ax.legend(handles, labels, fontsize=fs, ncol=2, loc='upper center')

    # add space at the top of the plot for the legend
    ymax = max(transcript_numbers)
    ymax = ymax*1.3
    ax.set_ylim(0, int(ymax))
    ax.set_xlim(-0.5, len(xtick_labels)-0.5)

    plt.tight_layout()

    plt.savefig(filename, dpi = 300, transparent = True)
    print("Figure saved in the current working directory directory as: "+filename)



if __name__ == "__main__":
    
    SE_stats_paths_dict = filepaths_SE_stats()
    data = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/Lome_RNA_annot_comparison/"

    numbers = {}
    for annotation, path in SE_stats_paths_dict.items():
        numbers[annotation] = parse_SE_stats(path)

    # print(numbers)
    plot_single_exon_no_species_specific_three_annot(numbers=numbers, filename=f"{data}Cmac_Lome_annot_comparison.png")