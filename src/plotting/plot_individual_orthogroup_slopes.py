"""
for all orthogroups of choice:
    compute the slope of the orthogroup (14 datapoints, one for each species with [GF size, genome size])
plot the slope of all of those together, and highlight OGs with high slopes
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from statistics import mean
import numpy as np
import pandas as pd
import parse_gff as gff
import parse_orthogroups as OGs

