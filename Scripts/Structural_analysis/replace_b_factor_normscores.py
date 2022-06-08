#!/usr/bin/env python
# coding: utf-8

# ## Replacing B factor with normalized scores
# 
# This script will take an input PDB and replace its b-factor column with values from a table for visualization in PyMOL or ChimeraX

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

# Define helper functions
def parse_pdb_line(pdb_line):
    '''This function will receive a line from a PDB file and parse it as a list. It will do so based on the
    PDB format explanation from this site:

    https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html.
    '''
    atom = pdb_line[0:4].strip(' ')
    atom_num = pdb_line[6:11].strip(' ')
    atom_name = pdb_line[12:16].strip(' ')
    resname = pdb_line[17:20].strip(' ')
    chain = pdb_line[21]
    res_num = pdb_line[22:26].strip(' ')
    x = pdb_line[30:38].strip(' ')
    y = pdb_line[38:46].strip(' ')
    z = pdb_line[46:54].strip(' ')

    return [atom, atom_num, atom_name, resname, chain, res_num, x, y, z]

##################################

def replace_b_factor(pdb_infile, in_dict, outfile):
    '''This function uses an input PDB file and replaces b-factors with values from
    a dictionary.'''
    in_handle = open(pdb_infile, 'r')
    out_handle = open(outfile, 'w')

    for line in in_handle:
        if line.startswith('ATOM'):

            # Parse the line
            parsed_line = parse_pdb_line(line)

            chain = parsed_line[4]
            resid = int(parsed_line[5])

            # Replace the b-factor
            dict_sasa = str(in_dict[chain][resid]) +'0'
            final_sasa = (6-len(dict_sasa))*' ' + dict_sasa

            # final_line = line.replace(line[60:66], final_sasa)
            final_line = line[0:60] + final_sasa + line[66:]
            out_handle.write(final_line)

    # Close the outfile
    out_handle.close()

##################################

def replace_admin(pdb_file, outfile, score_file):
    '''This function receives the files and feeds the information to the replace_b_factor function.'''
    ## Read the scores
    diffNormScores = pd.read_csv(score_file, sep = '\t', index_col = None)
    
    print('Maximum:', max(diffNormScores[diffNormScores.columns[1]]))
    print('Minimum:', min(diffNormScores[diffNormScores.columns[1]]))
    print('-------------')
    
    ## Convert the pandas dataframe to a dictionary compatible with the previous function
    diffNormScores_dict = {'A':{}, 
                          'B':{},
                          'C':{},
                          'D':{}}

    # Assign a zero to position one (no data for it)
    diffNormScores_dict['A'][1] = round(0, 2)
    diffNormScores_dict['B'][1] = round(0, 2)
    diffNormScores_dict['C'][1] = round(0, 2)
    diffNormScores_dict['D'][1] = round(0, 2)

    for index, row in diffNormScores.iterrows():

        ### For the new file
        ## Add values to dictionary
        diffNormScores_dict['A'][int(row['Position'])] = round(row[1], 2)
        diffNormScores_dict['B'][int(row['Position'])] = round(row[1], 2)
        diffNormScores_dict['C'][int(row['Position'])] = round(row[1], 2)
        diffNormScores_dict['D'][int(row['Position'])] = round(row[1], 2)

    # Replace the b-factors and write the output file
    replace_b_factor(pdb_file, diffNormScores_dict, outfile)
    


# ## Use the admin function to save structures for maximum, minimum, and mean deltaS at each position (s_weak - s_opt)
## Work with the files for the minimum effects
# This PDB file contains the AF2 model translated to the coordinates of the biological assembly
# of 2RK1 based on a structural alignment.
pdb_file = '../../Data/Structural_data/DfrB1_alphafold_rank1_2RK1Coords.pdb'

## Max delta s
replace_admin(pdb_file = pdb_file, 
             outfile = '../../Figures/Chimerax_figures/2rk1_max_deltaS_ara0.01_ara0.2_af2.pdb',
             score_file = '../../Figures/Chimerax_figures/max_deltaS_ara0.2_ara0.01.txt')

## Min delta s
replace_admin(pdb_file = pdb_file, 
             outfile = '../../Figures/Chimerax_figures/2rk1_min_deltaS_ara0.01_ara0.2_af2.pdb',
             score_file = '../../Figures/Chimerax_figures/min_deltaS_ara0.2_ara0.01.txt')

## Mean delta s
replace_admin(pdb_file = pdb_file, 
             outfile = '../../Figures/Chimerax_figures/2rk1_mean_deltaS_ara0.01_ara0.2_af2.pdb',
             score_file = '../../Figures/Chimerax_figures/mean_deltaS_ara0.2_ara0.01.txt')

