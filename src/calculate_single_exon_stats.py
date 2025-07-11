### calculate gff statistics focused on single- and multi-exon genes
## write the results to stdout
"""
This is just a wrapper script around print_single_exon_stats() in parse_gff.py to make automating it a little easier
useage: 
python3 calculate_single_exon_stats.py path/to/annotation.gff [True|False]
"""

import parse_gff as gff
import sys

if True:
    try:
        orthodb_annot = sys.argv[1]
        # native_annot = sys.argv[2]
        arg =  sys.argv[2]
    except Exception as e:
        print(f"error:{e}")
        print(f"\narguments: {sys.argv}")
        print(f"useage: python3 calculate_single_exon_stats.py annotation.gff [True|False]")
        sys.exit(1)

    if arg.lower() == "true":
        boolean_value = True
    elif arg.lower() == "false":
        boolean_value = False
    else:
        print("Invalid third argument. Please pass 'True' or 'False'.")
        sys.exit(1)

    gff.print_single_exon_stats(orthodb_annot, include_list=boolean_value) # include or exclude a list of all the single-exon IDs

    # write to stdout, so if you include the list a long list of gene IDs will be written to stdout
    # parse_gff.print_single_exon_stats(native_annot, include_list=boolean_value)
