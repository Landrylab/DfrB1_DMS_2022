#!/usr/bin/env python
# coding: utf-8

# # Calculate Shannon entropy
# 
# This script will look at the DHFR orthologs we got from jackhmmer and look at the variants at each position. It will also filter out some of the lowest quality sequences and run Evol to calculate Shannon entropy at each position.
# 

# Load libraries
import csv
import numpy as np
import pandas as pd
import re

from prody import *
print(prody.__name__, prody.__version__)

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


# Define variables
jackhmmer_tsv_path = '../../Data/Evolutionary_analysis/DHFR_orthologs_jackhmmer.tsv'

# Load the jackhmmer data to filter out some of the low-confidence sequences
jackhmmer_data = pd.read_csv(filepath_or_buffer = jackhmmer_tsv_path, sep = '\t')

jackhmmer_data


# Filter the sequences based on the following criteria:
# - E-score below 1e-6
# - Alignment length (alignment target end - alignment target start) greater than 50 (full sequence would be 78)
# - Sequence length below 100

tmp = jackhmmer_data.loc[jackhmmer_data['Target Ali. End'] - jackhmmer_data['Target Ali. Start'] > 50]
tmp2 = tmp.loc[tmp['Target Length'] < 100]
filtered_dataset = tmp2.loc[tmp2['E-value'] < 1e-6]


filtered_dataset.sort_values(by = ['E-value'])


# ## Look at diversity (Shannon entropy) per position
# 

# Load the alignment
msa_dhfr = parseMSA('../../Data/Evolutionary_analysis/all_sequences_DHFR_orthologs_aln.fasta')

# Use the refine function with the human sequence and different occupancy thresholds
checks = []

for i in range(0,100,5):
    threshold = float(i)/100
    msa_dhfr_refined = refineMSA(msa_dhfr, label = 'DfrB1_Joelle', rowocc=threshold)
    checks.append(len(msa_dhfr_refined))

y_pos = range(0,100,5)

msa_dhfr_refined = refineMSA(msa_dhfr, label = 'DfrB1_Joelle', rowocc=0.80)
occupancy_data = pd.DataFrame(np.column_stack((y_pos, checks)), columns=['Occupancy', 'Sequences'])
occupancy_data.to_csv(path_or_buf=os.path.join('../../Data/Evolutionary_analysis/', 'DHFR_occupancy_filter.txt'), sep = '\t', index = False)

get_ipython().run_line_magic('matplotlib', 'inline')
y_pos = range(0,100,5)
fig, ax = plt.subplots()
plt.bar(y_pos, checks, alpha = 0.5, width = 3, align = 'center')
plt.xticks(y_pos)
plt.ylabel('Number of sequences')
plt.xlabel('Percentage occupancy (%)')
plt.title('Maintained sequences, DHFR')

# Write a function to save the refined alignments
def write_prody_aln(ref_alignment, outfile):
    '''This function receives a refined alignment variable (ProDy's MSA class) and a path to an output file.
    It will save the alignment in the fasta format.
    '''
    # Start the output variable
    records = []
    id_list = ref_alignment.getLabels()
    # Loop through each of the records
    for i in range(len(id_list)):
        new_id = id_list[i]
        full_sequence = ''
        seq_array = ref_alignment.getArray()[i]
        for j in range(len(seq_array)):
            full_sequence = full_sequence + seq_array[j].decode('UTF-8')
        new_record = SeqRecord(Seq(full_sequence), id = new_id)
        records.append(new_record)
    # Write the sequences
    SeqIO.write(records, outfile, 'fasta')

# Save the refined alignment
write_prody_aln(msa_dhfr_refined, '../../Data/Evolutionary_analysis/DataS1_DfrB1_DMS_2022.fasta')

# Calculate Shannon entropy for each position
entropy = calcShannonEntropy(msa_dhfr_refined)
np.savetxt(os.path.join('../../Data/Evolutionary_analysis/', 'DHFR_entropy.txt'), entropy)

# Show the figure for Shannon entropy
showShannonEntropy(entropy)

positions = range(1, 79)
figure(figsize=(10, 6), dpi=80)
plot(positions, entropy)
plt.ylabel('Shannon entropy')
plt.xlabel('Position')

