#!/usr/bin/env python
# coding: utf-8

# ## AlphaFold2 pLDDT values
# 
# This script will read the pLDDT confidence values from the AF2 model of DfrB1 and save them to a simpler table to have uniform figures in R.

# Load libraries
import re
import os
from collections import OrderedDict
import math
from Bio.PDB import *
from Bio import SeqIO
from Bio.Seq import *
from Bio.SeqRecord import *
import csv
import numpy as np
import pandas as pd

# Read the PDB file
parser = PDBParser()
structure = parser.get_structure("DfrB1", 
                                 "../../Data/Structural_data/Dfrb1_6ca29_unrelaxed_rank_1_model_1.pdb")

structure

# Define a dictionary to go from three-letter residues to one-letter
aa_three2one = {'ALA': 'A', 
               'ARG': 'R',
                'ASN':'N',
                'ASP':'D',
                'CYS':'C',
                'GLU':'E',
                'GLN':'Q',
                'GLY':'G',
                'HIS':'H',
                'ILE':'I',
                'LEU':'L',
                'LYS':'K',
                'MET':'M',
                'PHE':'F',
                'PRO':'P',
                'SER':'S',
                'THR':'T',
                'TRP':'W',
                'TYR':'Y',
                'VAL':'V'
               }


residue_list = []
pLDDT_list = []

curr_pos = 1

# Loop through the residues in the structure
for model in structure:
    for chain in model:
        for residue in chain:
            new_residue = aa_three2one[residue.get_resname()] + str(curr_pos)
            residue_list.append(new_residue)
            curr_pos = curr_pos + 1
            for atom in residue:
                # The pLDDT is saved as the B-factor in the file and it is the same for
                # all atoms of the same residue
                new_pLDDT = atom.get_bfactor()
                pLDDT_list.append(new_pLDDT)
                break
        break
        print('----')


print(residue_list)
print(pLDDT_list)

# Organize the data in a dataframe
data_dict = {
    'Residue': residue_list,
    'pLDDT': pLDDT_list
}
new_df = pd.DataFrame(data_dict)
new_df

# Save the dataframe
new_df.to_csv('../../Data/Structural_data/AF2model1_pLDDT.txt',
              sep = '\t', index = False)

