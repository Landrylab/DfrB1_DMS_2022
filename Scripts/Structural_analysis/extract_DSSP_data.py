#!/usr/bin/env python
# coding: utf-8

# ## Extract secondary structure annotations from the DSSP file
# 
# The meaning of the annotations is explained here: http://www.csb.yale.edu/userguides/databases/dssp/dssp_man.html
# 
# - H = alpha-helix
# - B = beta-bridge residue
# - E = extended strand (in beta ladder)
# - G = 3/10-helix
# - I = 5-helix
# - T = H-bonded turn
# - S = bend
# 

# Load libraries
import csv
import numpy as np
import pandas as pd
import re

from matplotlib.pylab import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import matplotlib.colors as mcol
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import OrderedDict

code_dict = {
    'H':'Alpha helix',
    'B':'Beta bridge',
    'E':'Beta ladder', 
    'G':'3/10 helix',
    'I':'5-helix',
    'T':'H-bonded turn',
    'S':'Bend'
}

# A handle to open the file
handle = open('../../Data/Structural_data/DHFR_2rk1_from_file_bio.dssp', 'r')

# Start a dictionary to keep track of the values
sec_struc_annotations = {}

# Initialize all positions in the sequence as NA
for i in range(1, 79):
    sec_struc_annotations[i] = [i, 'Missing', 'Missing']

# A boolean to know when to start looking at the data
bool_data = False

# Loop through the lines
for line in handle:
    
    # Looking for the start of the data
    if not bool_data:
        if line.startswith('  #'):
            # The hashtag marks the header, start reading data from the 
            # following line
            bool_data = True
            continue
        else:
            continue
    else:
        # Start looking at the data from here
    
        # Split by one or more instances of whitespace
        split_line = re.split(pattern = '\s+', string = line)
        
        # Columns 2 and 5 (0-based) contain the position and the secondary structure annotation
        # However, positions with no annotation in column 5 will skip that column,
        # make sure the value in that column is a secondary structure annotation
        
        # If column 2 is "!*", it marks the end of the first subunit. Stop there since 
        # the structure is a symmetric homomer
        if split_line[2] == '!*':
            break
        else:
            position = int(split_line[2])
        sec_struc = split_line[5]
        
        # Save solvent accessibility
        solv_acc = int(line[34:39].strip())
        
        # If this is a secondary structure annotation
        if sec_struc in ('H', 'B', 'E', 'G', 'I', 'T', 'S'):
            sec_struc_annotations[position] = [position, code_dict[sec_struc], solv_acc]
        else:
            sec_struc_annotations[position] = [position, 'none', solv_acc]


handle.close()

sec_struc_annotations

# Save to a pandas dataframe
df_sec_struc = pd.DataFrame.from_dict(sec_struc_annotations, orient='index', columns = ['Position', 'Secondary_structure', 'Solvent_accessibility'])
df_sec_struc

# Save dataframe
df_sec_struc.to_csv(path_or_buf = '../../Data/Structural_data/DHFR_2rk1_DSSP_table_bio.txt', sep = '\t', index = False)

