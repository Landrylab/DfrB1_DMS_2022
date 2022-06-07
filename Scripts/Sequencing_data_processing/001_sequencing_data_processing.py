#!/usr/bin/env python
# coding: utf-8

# # 001_sequencing_data_processing
# 
# This script demultiplexes the MiSeq data and NovaSeq data based on the metadata file in supplementary table S2 (sample description). It then merges the reads and aligns to the reference to count the number of reads that map to each mutant.
# 
# This script assumes the following folder layout
# - Home folder
#     - Scripts/Sequencing_data_processing folder with this script
#     - Data/Sequencing_data folder with the raw sequencing data from the SRA
#     - Data/uni_amplicon_sequences with the indexes needed to demultiplex the data
#     - Data/Analysis_MiSeq for the results of the demultiplexing for the MiSeq data
#     - Data/Analysis_NovaSeq for the results of the demultiplexing for the NovaSeq data
#     
# This script assumes the following programs are installed (paths to the executables need to be specified below):
# - FastQC
# - Trimmomatic
# - Bowtie 1
# - Pandaseq
# - vsearch
# - Needle


import os
import subprocess

import pandas as pd
print(pd.__name__, pd.__version__)

import numpy as np
print(np.__name__, np.__version__)

import matplotlib.pyplot as plt
import matplotlib
print(matplotlib.__name__, matplotlib.__version__)

import scipy.stats as stats
import scipy
print(scipy.__name__, scipy.__version__)

import re
print(re.__name__, re.__version__)

from collections import Counter

import sys
print(sys.version)

import seaborn as sns
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec

from Bio import SeqIO

## Define paths of executables
fastqc_path = '/path/to/fastqc'
trimmomatic_path = '/path/to/Trimmomatic'
bowtie_build_path = '/path/to/bowtie-build'
bowtie_path = '/path/to/bowtie'

## Define some helper functions
def reverse_complement(dna):
    """ function that reverse complements DNA
    dna: input dna sequence
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def get_dict_of_seq(fasta_file):
    """ function that converts a fasta file to a dictionnary of sequences
    fasta_file: the input fasta file
    """
    
    file_fasta_dict = {}
    # output dict of imported seqs
    
    with open(fasta_file, 'r') as fasta:    
        for line in fasta:
            # loops through the file

            if line.startswith('>') == True:
                seq_info = line.strip('>').strip('\n').split('\t')[0]
                file_fasta_dict[seq_info] = ''
                # checks if seq header

            else:
                file_fasta_dict[seq_info] += line.strip('\n')
                # If not, append nt to current seq
                
    return file_fasta_dict

universal_seqs_fasta = '../../Data/uni_amplicon_sequences/universal_amplicon.fa'
universal_seqs = get_dict_of_seq(universal_seqs_fasta)
print(universal_seqs)

indexes_fasta = '../../Data/uni_amplicon_sequences/indexes.fa'
indexes = get_dict_of_seq(indexes_fasta)
print(indexes)

degen_seq = 'AAAAA'


# ## Work with the MiSeq data

# ### Check sequencing quality and select for seqs w/ correct length

experiment_name = 'R67_DMS'

path_to_R1_file = '../../Data/Sequencing_data/libs-R67_S1_L001_R1_001.fastq.gz'
path_to_R2_file = '../../Data/Sequencing_data/libs-R67_S1_L001_R2_001.fastq.gz'


## FastQC quality control
subprocess.check_output(fastqc_path + ' ' + path_to_R1_file, shell=True)
subprocess.check_output(fastqc_path + ' ' +path_to_R2_file, shell=True)


# Initialize intermediate folders for the analysis
amplicon_sequences_path = '../../Data/Analysis_MiSeq/amplicon_sequences'
os.makedirs(amplicon_sequences_path, exist_ok = True)

bowtie_indexes_path = '../../Data/Analysis_MiSeq/bowtie_indexes'
os.makedirs(bowtie_indexes_path, exist_ok = True)

temp_path = '../../Data/Analysis_MiSeq/temp'
os.makedirs(temp_path, exist_ok = True)

## Trim MiSeq reads with Trimmomatic
trim_call = 'java -jar ' + trimmomatic_path + ' PE -threads 6 -trimlog log.txt '
trim_call += path_to_R1_file+' '+path_to_R2_file+' '
trim_call += '-baseout ../../Data/Analysis_MiSeq/temp/minlen299.fastq MINLEN:299 CROP:270'

trim_call

subprocess.check_output(trim_call, shell=True)

## FastQC quality control on the trimmed reads
subprocess.check_output('fastqc '+'../../Data/Analysis_MiSeq/temp/minlen299_1P.fastq', shell=True)
subprocess.check_output('fastqc '+'../../Data/Analysis_MiSeq/temp/minlen299_2P.fastq', shell=True)

# ### check abundance of libraries

# Output files
plate_pcr_amplicon_for = os.path.join(amplicon_sequences_path, 'plate_pcr_for.fa')
plate_pcr_amplicon_rev = os.path.join(amplicon_sequences_path, 'plate_pcr_rev.fa')

for_plate_indexes = [1,2,3]

rev_plate_indexes = [1,2,3]

print(for_plate_indexes)
print(rev_plate_indexes)

# Make a bowtie index of the forward plate indices
with open(plate_pcr_amplicon_for, 'w') as fasta_out:
    
    for index in for_plate_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = degen_seq+index_seq+universal_seqs['plate_fivep_sticky']+'\n'
        fasta_out.write(amplicon_frag)
        
        
for_plate_index = os.path.join(bowtie_indexes_path, experiment_name+'_for_plate')

for_plate_bowtie_index_call = bowtie_build_path + ' -f -r -o 4 ' +plate_pcr_amplicon_for+' '+for_plate_index
for_plate_indexing_log = subprocess.check_output(for_plate_bowtie_index_call, shell=True)
        
# Make a bowtie index of the forward plate indices    
with open(plate_pcr_amplicon_rev, 'w') as fasta_out:
    
    for index in rev_plate_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = degen_seq + index_seq+ reverse_complement(universal_seqs['plate_threep_sticky'])+'\n'
        fasta_out.write(amplicon_frag)    
    
rev_plate_index = os.path.join(bowtie_indexes_path, experiment_name+'_rev_plate')
rev_plate_bowtie_index_call = bowtie_build_path + ' -f -r -o 4 '+plate_pcr_amplicon_rev+' '+rev_plate_index
rev_plate_indexing_log = subprocess.check_output(rev_plate_bowtie_index_call, shell=True)


# ## Align for barcode

## Bowtie arguments:
# -t print time
# -v maximum number of mismatches tolerated
# -p number of threads
# -k number of valid alignments to report
# --trim3 trim x nt from 3' end
# --trim5 trim x nt from 5' end
# --norc do not align to reverse complement
# --chunkmbs Mbs of memory for each thread

test_align_call = bowtie_path + ' -t -v 3 -p 6 -k 1 --trim3 250 --trim5 5 --norc --chunkmbs 256 '

test_align_call += for_plate_index+' '
test_align_call += os.path.join(temp_path, 'minlen299_1P.fastq') + ' '
test_align_call += os.path.join(temp_path, 'plate_for_align.txt')

print (test_align_call)

subprocess.check_output(test_align_call, shell = True)


test_align_call = bowtie_path + ' -t -v 3 -p 6 -k 1 --trim3 250 --trim5 5 --norc --chunkmbs 256 '

test_align_call += rev_plate_index+' '
test_align_call += os.path.join(temp_path, 'minlen299_2P.fastq') + ' '
test_align_call += os.path.join(temp_path, 'plate_rev_align.txt')
print (test_align_call)
subprocess.check_output(test_align_call, shell = True)

plate_align_for_output = os.path.join(temp_path, 'plate_for_align.txt')
plate_align_rev_output = os.path.join(temp_path, 'plate_rev_align.txt')


plate_index_found_for = {}
plate_index_found_rev = {}
# empty containers that will hold    read_name: index    data pairs for the forward and reverse reads that
# had valid alignments

with open(plate_align_for_output, 'r') as for_index_positives:
    for line in for_index_positives:
        # opens P1 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info, not the illumina annotations in the read 
        # name. See http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/
        #           Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm
        # For info on illumina read name convention
        plate_index_found_for[read_name] = line[2]
        # add read P1 index to dict

with open(plate_align_rev_output, 'r') as rev_index_positives:
    for line in rev_index_positives:
        # opens P2 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info
        plate_index_found_rev[read_name] = line[2]
        # add read P2 index to dict


plate_both = list(set(plate_index_found_for.keys()) & set(plate_index_found_rev.keys()))

print(len(plate_both))


# This block looks for the different indexes in the data

RC_pcr_amplicon_for = os.path.join(amplicon_sequences_path, 'RC_pcr_for.fa')
RC_pcr_amplicon_rev = os.path.join(amplicon_sequences_path, 'RC_pcr_rev.fa')

RC_for_indexes = range(1,16)

RC_rev_indexes = range(1,16)

print(RC_for_indexes)
print(RC_rev_indexes)


with open(RC_pcr_amplicon_for, 'w') as fasta_out:
    
    for index in RC_for_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = universal_seqs['plate_fivep_sticky']+index_seq+universal_seqs['RC_fivep_sticky']+'\n'
        fasta_out.write(amplicon_frag)
        
        
for_RC_index = os.path.join(bowtie_indexes_path, experiment_name+'_for_RC')
for_RC_bowtie_index_call = bowtie_build_path + ' -f -r -o 4 '+RC_pcr_amplicon_for+' '+for_RC_index
for_RC_indexing_log = subprocess.check_output(for_RC_bowtie_index_call, shell=True)
        
        
with open(RC_pcr_amplicon_rev, 'w') as fasta_out:
    
    for index in RC_rev_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = reverse_complement(universal_seqs['plate_threep_sticky']) + index_seq
        amplicon_frag += reverse_complement(universal_seqs['RC_threep_sticky'])+'\n'
        fasta_out.write(amplicon_frag)    
    
rev_RC_index = os.path.join(bowtie_indexes_path, experiment_name+'_rev_RC')
rev_RC_bowtie_index_call = bowtie_build_path + ' -f -r -o 4 '+RC_pcr_amplicon_rev+' '+rev_RC_index
rev_RC_indexing_log = subprocess.check_output(rev_RC_bowtie_index_call, shell=True)

test_align_call = bowtie_path + ' -t -v 3 -p 6 -k 1 --trim3 215 --trim5 28 --norc --chunkmbs 256 '
# adjusted for 50pb trimming

test_align_call += for_RC_index+' '
test_align_call += os.path.join(temp_path, 'minlen299_1P.fastq') + ' '
test_align_call += os.path.join(temp_path, 'RC_for_align.txt')

subprocess.check_output(test_align_call, shell = True)

test_align_call = bowtie_path + ' -t -v 3 -p 6 -k 1 --trim3 217 --trim5 26 --norc --chunkmbs 256 '
# adjusted for 50pb trimming

test_align_call += rev_RC_index+' '
test_align_call += os.path.join(temp_path, 'minlen299_2P.fastq') + ' '
test_align_call += os.path.join(temp_path, 'RC_rev_align.txt')

subprocess.check_output(test_align_call, shell = True)

RC_align_for_output = os.path.join(temp_path, 'RC_for_align.txt')
RC_align_rev_output = os.path.join(temp_path, 'RC_rev_align.txt')


index_found_for_RC = {}
index_found_rev_RC = {}
# empty containers that will hold    read_name: index    data pairs for the forward and reverse reads that
# had valid alignments

with open(RC_align_for_output, 'r') as for_index_positives:
    for line in for_index_positives:
        # opens P1 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info, not the illumina annotations in the read 
        # name. See http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/
        #           Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm
        # For info on illumina read name convention
        index_found_for_RC[read_name] = line[2]
        # add read P1 index to dict


with open(RC_align_rev_output, 'r') as rev_index_positives:
    for line in rev_index_positives:
        # opens P2 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info
        index_found_rev_RC[read_name] = line[2]
        # add read P2 index to dict


both_indexes_RC = list(set(index_found_for_RC.keys()) & set(index_found_rev_RC.keys()))

print(len(both_indexes_RC))

# Valid plates for the MiSeq data
valid_plates = [tuple((3,1))]

def plot_reads_in_plate(plate_both, RC_both, plate_index_found_for, plate_index_found_rev, index_found_for_RC, index_found_rev_RC, plate,
                       silent=True):
    
    plate_rows = sorted(set(index_found_for_RC.values()))
    plate_columns = sorted(set(index_found_rev_RC.values()))
    
    read_plate_RC_dict = {}
    grid_dict = {}
    
    for row in plate_rows:
        
        grid_dict[int(row)] = {}
                
        for column in plate_columns:
                        
            grid_dict[int(row)][int(column)] = 0
            
    if silent == False:
        
        print(grid_dict)
    
    plate_and_rc_both = list((set(plate_both)&set(both_indexes_RC)))
    
    for read in plate_and_rc_both:
        
        for_index_plate = int(plate_index_found_for[read])
        rev_index_plate = int(plate_index_found_rev[read])
        
        if tuple((int(for_index_plate), int(rev_index_plate))) == plate:
                        
            for_index_RC = int(index_found_for_RC[read])
            rev_index_RC = int(index_found_rev_RC[read])
            
            grid_dict[for_index_RC][rev_index_RC] += 1
            
            read_plate_RC_dict[read] = [int(for_index_plate), int(rev_index_plate), for_index_RC, rev_index_RC]
            
            
            
    plate_df = pd.DataFrame.from_dict(grid_dict, orient='index')
    return plate_df, read_plate_RC_dict
            
            
# Keep in mind I am only using the first one of the valid plates array
# since it's the only one I need for this analysis
F_1 = plot_reads_in_plate(plate_both, both_indexes_RC, plate_index_found_for, plate_index_found_rev, index_found_for_RC, 
                                index_found_rev_RC, valid_plates[0], silent = False)

F_1[0].head(10)

col_order = [1,2,3,4,5,6,7,8,9,10]
plt.figure(figsize=(14/1.5, 10/1.5))


sns.heatmap(F_1[0][col_order], vmax=6400, cmap='viridis',)
plt.xlabel('Column Barcode', fontsize = 16)
plt.ylabel('Row Barcode', fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.axhline(2, color='white')
plt.axhline(4, color='white')
plt.axhline(6, color='white')
plt.axhline(8, color='white')
plt.title('Fragment 1')

all_plate_index_dicts = []

for plate in valid_plates:
    
    mapping = plot_reads_in_plate(plate_both, both_indexes_RC, plate_index_found_for, plate_index_found_rev, index_found_for_RC, 
                                index_found_rev_RC, plate, silent=True)[1]
    
    all_plate_index_dicts.append(mapping)

demux_reads_path = '../../Data/Analysis_MiSeq/demultiplexed_reads'
os.makedirs(demux_reads_path, exist_ok = True)


# ## Look for the row column barcodes
plate_RC_combi_to_dict = {}

# The plate indexes are always 3 and 1
plate_for = 3
plate_rev = 1

# There are 16 sample IDs in my data (0 through 15)
# Combinations of RCs are given in Isa's table
for i in range(16):
   
    # Takes values from 1 to 8, then restarts from 1 to 8
    RC_for = (i % 8) + 1
    
    # This RC is 1 for the first 8 samples and 2 for the last 8 samples
    # I will use the ceiling function
    RC_rev = int(np.ceil((i + 1) / 8))
    
    new_key = 'pool_' + str(i)
    
    plate_RC_combi_to_dict[new_key] = [plate_for, plate_rev, RC_for, RC_rev]

    
plate_RC_combi_to_dict


# ## Keep in mind that this block appends to the files in the demultiplexed reads folder

plate_RC_combi_to_dict = {}

out_fastq_path_dict_for_df = {}

out_fastq_path_dict_demult = {}

dir_path = demux_reads_path

# The plate indexes are always 3 and 1
plate_for = 3
plate_rev = 1

# There are 16 sample IDs in my data (0 through 15)
# Combinations of RCs are given in Isa's table
for i in range(16):
   
    # Takes values from 1 to 8, then restarts from 1 to 8
    RC_for = (i % 8) + 1
    
    # This RC is 1 for the first 8 samples and 2 for the last 8 samples
    # I will use the ceiling function
    RC_rev = int(np.ceil((i + 1) / 8))
    
    sample = 'pool_' + str(i)
    
    plate_RC_combi_to_dict[sample] = [plate_for, plate_rev, RC_for, RC_rev]

    out_path = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)

    out_path_1 = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)+'_for'+'.fastq'
    out_path_2 = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)+'_rev'+'.fastq'

    out_fastq_path_dict_for_df[sample] = out_path

    # Save the paths to the files so that they can be used in the following function
    out_fastq_path_dict_demult[tuple([plate_for, plate_rev, RC_for, RC_rev])] = out_path

    # Creates empty files to append to them
    create_file_for = 'touch '+ out_path_1
    create_file_rev = 'touch '+ out_path_2

    subprocess.check_output(create_file_for, shell=True)
    subprocess.check_output(create_file_rev, shell=True)
    


for plate_dict in all_plate_index_dicts:
    
    print(len(plate_dict.keys()))

    for read in list(plate_dict.keys()):

        index_combi = plate_dict[read]

        if index_combi not in plate_RC_combi_to_dict.values():

            plate_dict.pop(read, None)
        
    print(len(plate_dict.keys()))

## This function will add the reads to the files we created above
def split_reads_from_fastq(file, file_type, plate_mapping_list):
    
    to_write = 'none'
    read_name = 'none'
    
    
    with open(file, 'r') as source_fastq:
        
        for line in source_fastq:

            # If this line indicates the name of the read
            if line.startswith('@M') == True:
                
                if to_write == 'none':
                    
                    read_name = line.split(' ')[0].strip('@')
                    to_write = line

                # Otherwise, we are looking at the information about the read
                else:
                    
                    # Check in which plate it is
                    for plate_read_dict in plate_mapping_list:
                        
                                                
                        if read_name in plate_read_dict.keys():

                            indexes = tuple(plate_read_dict[read_name])
                            
                            # This is where we distinguish between the files for forward and reverse reads
                            if file_type == 'forward':
                                
                                file_path = out_fastq_path_dict_demult[indexes]+'_for.fastq'
                                
                            elif file_type == 'reverse':
                                
                                file_path = out_fastq_path_dict_demult[indexes]+'_rev.fastq'

                            
                            
                            with open(file_path, 'a') as dest:
                                
                                   dest.write(to_write)

                            break





                    read_name = line.split(' ')[0].strip('@')

                    to_write = line
                

            else:

                to_write += line

                
# Separate the forward from the reverse reads
split_reads_from_fastq(os.path.join(temp_path, 'minlen299_1P.fastq'), 'forward', all_plate_index_dicts)
split_reads_from_fastq(os.path.join(temp_path, 'minlen299_2P.fastq'), 'reverse', all_plate_index_dicts)


# ## Merge the reads
os.makedirs('./../../Data/Analysis_MiSeq/merged_reads', exist_ok = True)

merged_file_list = []

# Loop through the files
for filepath in out_fastq_path_dict_demult.values():

    filepath_for = filepath+'_for.fastq'
    filepath_rev = filepath+'_rev.fastq'

    sample_read_count = 0

    with open(filepath_for, 'r') as source:

        for line in source:

            if line.startswith('@M'):

                sample_read_count += 1

    print(sample_read_count)

    if sample_read_count >= 1:

        # Derive the name of the output file based on the input file
        filepath_out_prefix = './../../Data/Analysis_MiSeq/merged_reads/'
        
        # Use regular expressions to get the name of the oputput file
        pool_number = re.search('pool_\d+', filepath).group(0)
        
        filepath_out = os.path.join(filepath_out_prefix, pool_number + '_merged.fasta')
        print(filepath_out)

        # Add to the list of merged files so we can keep working with it
        merged_file_list.append(filepath_out)
        
        ## Pandaseq arguments:
        # -f input file with forward reads
        # -r input file with reverse reads
        # -L maximum length for a sequence
        # -O maximum overlap for a sequence
        # -k kmers
        # -B allow unbarcoded sequences
        # -N eliminate all sequences with unknown nucleotides
        # -t threshold (minimum probability a sequence must have to assemble)
        # -T threads
        # -w output file in fasta.bz2 format
        panda_seq_call = 'pandaseq -f '+filepath_for+' -r '+filepath_rev+ ' -L 550 -O 400 -k 4 -B -N -t 0.5 -T 6 -w '+ filepath_out

        print (panda_seq_call)

        subprocess.check_output(panda_seq_call, shell=True)
        
        print('--------')


# ## Aggregate with vsearch
outfile_trimmed_merged_prefix = './../../Data/Analysis_MiSeq/trimmed_merged_reads/'

os.makedirs(outfile_trimmed_merged_prefix, exist_ok = True)

trimmed_merged_list = []

## Loop through the merged files
for merged_file in merged_file_list:
    
    
    pool_number = re.search('pool_\d+', merged_file).group(0)
    filename_out = pool_number + '_trimmed_merged.fa'
    
    outfile_trimmed_merged = os.path.join(outfile_trimmed_merged_prefix, filename_out)
    
    print(outfile_trimmed_merged)
    
    trimmed_merged_list.append(outfile_trimmed_merged)
    
    ## vsearch parameters:
    # --fastx_filter trim and filter sequences in this file
    # --fastq_stripleft delete given number of bases from the 5' end
    # --fastq_stripright delete given number of bases from the 3' end
    # --fastaout output file
    vsearch_trim_call = 'vsearch --fastx_filter ' + merged_file + ' --fastq_stripleft 76 --fastq_stripright 76 --fastaout '
    vsearch_trim_call += outfile_trimmed_merged

    print (vsearch_trim_call)

    subprocess.check_output(vsearch_trim_call, shell=True)
    
    print('--------')

aggregate_list = []

aggregate_prefix = './../../Data/Analysis_MiSeq/aggregate_reads/'

## Loop through the trimmed merged files and aggregate
for trimmed_merged_file in trimmed_merged_list:
    
    pool_number = re.search('pool_\d+', trimmed_merged_file).group(0)
    filename_out = pool_number + '_aggregate.fa'
    
    outfile_aggregate = os.path.join(aggregate_prefix, filename_out)
    
    print(outfile_aggregate)
    
    aggregate_list.append(outfile_aggregate)
        
    ## vsearch parameters:
    # --derep_fulllength dereplicate sequences in the given fasta file
    # --relabel add a given prefix
    # --output output file
    # -- sizeout include abundance information

    vsearch_aggregate_call = 'vsearch --derep_fulllength ' + trimmed_merged_file + ' --relabel seq --sizeout --output '
    vsearch_aggregate_call += outfile_aggregate
    print (vsearch_aggregate_call)
    subprocess.check_output(vsearch_aggregate_call, shell=True)
    print('-----')


# ## Define the reference amplicon
os.makedirs('./../../Data/Analysis_MiSeq/amplicon_sequences/', exist_ok = True)

path_to_amplicons = './../../Data/uni_amplicon_sequences/miseq_r67_bacteria.fasta'

amplicon_info_dict = get_dict_of_seq(path_to_amplicons)

amplicon_length_dict = {}

amplicon_dict = {}

for amplicon in amplicon_info_dict:
    
    amplicon_name = amplicon.split('|')[0]

    amplicon_dict[amplicon_name] = amplicon_info_dict[amplicon]
    
    amplicon_fasta_path = './../../Data/Analysis_MiSeq/amplicon_sequences/'+amplicon_name+'.fasta'
    
    with open(amplicon_fasta_path, 'w') as dest:
        
        seq_ID = '>'+amplicon_name+'\n'
        seq = amplicon_dict[amplicon_name]+'\n'
        
        dest.write(seq_ID)
        dest.write(seq)
    

print (amplicon_info_dict.keys())
print (amplicon_dict)


# ## Align to the reference amplicon
os.mkdir('./../../Data/Analysis_MiSeq/amplicon_align')

## Define a function to do the alignments with needle
def needle_align_on_ref(ref_orf, filepath):
               
    ref_seq = amplicon_dict[ref_orf]
    
    ref_fasta_path = './../../Data/Analysis_MiSeq/amplicon_sequences/'+ref_orf+'.fasta '
    
    # Extract the name of the pool
    pool_number = re.search('pool_\d+', filepath).group(0)
    
    ## Needle arguments:
    # -auto Turn of prompts
    # -gapopen penalty for opening a gap
    # -asequence one of the two input sequences to align
    # -bsequence the second input sequence to align
    # -aformat output format for the alignment (markx10 is an easily parseable format http://emboss.sourceforge.net/docs/themes/alnformats/align.markx10 )
    # -outfile path to the output file
    
    # For the analysis_wt folder
    needle_out = './../../Data/Analysis_MiSeq/amplicon_align/'+ ref_orf + '_' + pool_number +'.needle'
    
    needle_call = 'needle -auto -gapopen 50 -asequence '+ ref_fasta_path
    
    needle_call += '-bsequence '+ filepath +' -aformat3 markx10 -outfile '+needle_out
    
    print(needle_call)
    
    subprocess.check_output(needle_call, shell = True)
          
    return

## Loop through the files of aggregate reads
for filepath in aggregate_list:
    needle_align_on_ref('MiSeq_R67_bacteria', filepath)


# ## Parse needle output and generate the dataframes with read counts

def parse_needle_output(needle_align_path):
    
    n_aligns = 0
    align_seqs_dict = {}       
   
    with open(needle_align_path, 'r') as source:

        current_align = ''
        current_qseq = ''
        current_sseq = ''
        qseq_done = 0

        for line in source:

            if line.startswith('>>>') == True:

                n_aligns +=1
                align_name = line.strip('>>>')

                if n_aligns != 1:

                    align_seqs_dict[current_align] = [current_qseq, current_sseq]
                    current_align = align_name
                    current_qseq = ''
                    current_sseq = ''
                    qseq_done = 0

                else:

                    current_align = align_name

            elif line.startswith(';') == False and line.startswith('>') == False and line.startswith('\n') == False and line.startswith('#') == False:

                if qseq_done == 1:
                    current_sseq += line.strip('\n').upper()

                else:
                    current_qseq += line.strip('\n').upper()

            elif line.startswith('#--') == True:

                align_seqs_dict[align_name] = [current_qseq, current_sseq]

            else:

                if qseq_done == 0 and current_qseq != '':
                    qseq_done =1            
                             
    return align_seqs_dict, n_aligns

def find_mutations(path, ref_orf):
    
    allele_dict = {}

    align_dict, align_count = parse_needle_output(path)

    for entry in list(align_dict.keys()):

        read_var_list = []

        query_seq = align_dict[entry][1]
        # aligned prot sequence of the strain

        align_ref = align_dict[entry][0]
        # aligned prot sequence of the reference

        gap_adjust = 0
        # value used to adjust the protein sequence index for the presence of insertions in the strain sequence vs the 
        # reference strain

        backtrack_adjust = 0

        temp_var = None
        # temporary variable to hold the sequence of an insertion or deletion as a string. When the gap ends, annotation 
        # will be added to strain_var_list

        indel_start = 0
        # position start of the indel annotation in the reference sequence, with adjustment for gap presence

        ref_seq_no_gaps = align_ref.replace('-','')
        # Make a copy of the reference sequence without gaps 
        
        # align_start = (amplicon_dict[ref_orf].index(ref_seq_no_gaps))+1
        align_start = (amplicon_dict[ref_orf].upper().index(ref_seq_no_gaps))+1
        # Look for the starting position of the gapless part of the reference sequence
        # This helps remove gaps at the start and the end of the alignment
        
        query_seq_no_gaps = len(query_seq.replace('-',''))
        # Make a copy of the query sequence without gaps
        
        for nt in range(0, len(align_ref)):
            # iterates through the entire alignment of the strain prot sequence

            if query_seq[nt] == '-':
                # detect a deletion variant

                # logic for indel detection/annotation:
                #
                # suppose we have this alignment  
                #
                # 1 2 3 4 5 6 7 8 9
                # A T - - A A A T G    strain variant: del gaps are indexed because the aa index is based on reference
                # A T K P A - - T G
                # 1 2 3 4 5     6 7    reference: insert gaps not indexed because aa positions do (actually don't?) exist in reference
                #
                # following this logic, every time an insertion is detected and annotated, the gap_adjust value is 
                # incremented by the length of the gap and used to adjust the variant mapping to make it match the 
                # reference index values. The indel aa postion is the first residue detected as part of the indel


                if indel_start == 0:
                    # checks if the character is the start or the continuation of a gap in the alignment

                    temp_var = 'del'+ align_ref[nt]
                    indel_start = (nt+1-gap_adjust)
                    # if it is, starts a new annotation entry with a start position compensated for previous insertions
                    # (if any)

                    backtrack_adjust += 1

                else:

                    temp_var += align_ref[nt]
                    # if it is not, adds the following aa to the deletion annotation

                    backtrack_adjust += 1


            elif align_ref[nt] == '-':
                # detects an insertion variant

                if indel_start == 0:
                    # checks if the character is the start or the continuation of a gap in the alignment

                    temp_var = 'ins'+ query_seq[nt]

                    indel_start = (nt+1-gap_adjust)
                    # if it is, starts a new annotation entry with a start position compensated for previous insertions
                    # (if any)

                    gap_adjust += 1
                    # increments the gap adjust for the this added aa in the strain sequence                   

                else:

                    temp_var += query_seq[nt]
                    # if it is not, adds the following aa to the insertion annotation

                    gap_adjust += 1
                    # increments the gap adjust for the this added aa in the strain sequence

            elif query_seq[nt] != align_ref[nt]:
                # detects a mismatch between the strain sequence and the reference

                variant = align_ref[nt]+'|'+str((nt+1-gap_adjust))+'|'+query_seq[nt]
                read_var_list.append(variant)
                # creates an annotation for the strain-reference aa mismatch and appends it to the list of 
                # annotations

            else:

                 if indel_start != 0:
                    # detects if there is currently an open gap entry. If there is, then the detected mismatch means 
                    # that it has now concluded

                    read_var_list.append(str((indel_start))+temp_var)
                    temp_var = None
                    indel_start = 0
                    # adds the indel annotation to the strain variant list and resets temporary variables for the next 
                    # indel entry

        if len(read_var_list)<25:  
            allele_dict[entry] = read_var_list, align_start
                           
    return allele_dict  

def get_variant_count_1(mutation_set, ref_seq, frag_start, codon_start, n_aa):
    
    variant_abundance_dict ={}
    variants = list(mutation_set.keys())
    codon_groups = {}
    codon = 1
    wt_count =0
    valid_seq=0
    
    # Calculate the end of the fragment with the frag start and the number of residues
    frag_end = frag_start + ((n_aa -1 )* 3) - 1
    
    for nt in range(0, (n_aa - 1)*3):
        pos = nt + frag_start
        
        if nt % 3 == 0:
            codon += 1
            
        codon_groups[pos] = codon
        
        variant_abundance_dict[codon] = {}       
        
    wt_codons = {}
        
    ref = amplicon_dict[ref_seq].upper()
    
    # Set the abundance of WT codons to not a number (nan) since they would be overrepresented 
    for aa in range(0, n_aa - 1):

        offset = frag_start - 1
        start = offset+(aa*3)
        wt_codon=ref[start:(start+3)]
        wt_codons[(aa+codon_start)] = wt_codon
        variant_abundance_dict[aa+codon_start][wt_codon]=np.nan 
    
    for variant in variants:
        
        var_info = variant.split(',')
        var_count =int(var_info[1].split(';')[1].strip('size='))      
        mut_list = mutation_set[variant][0]
        filtered_list = []
        
        # Only keep variants that appeared in more than 20 reads
        if var_count>=20:
        
            for mutation in mut_list:

                if 'del' in mutation:
                    mut_info = mutation.split('del')
                    mut_pos = int(mut_info[0])
                    
                    if mut_pos == 1:
                        mut_list.remove(mutation)
            
            # If not an indel
            if 'ins' not in str(mut_list) and 'del' not in str(mut_list):
                
                if len(mut_list) ==0:
                    wt_count += var_count
                    
                else:
                    mut_nt_list = []
                    out_list = []
                    
                    for mutation in mut_list:
                        
                        mut_pos = int(mutation.split('|')[1])
                        
                        if mut_pos >= frag_start and mut_pos <= frag_end:
                        
                            # mut_nt_list contains the list of mutated positions inside the coding sequence
                            mut_nt_list.append(codon_groups[mut_pos])
                            
                        else:
                            # out_list contains the list of mutated positions outside the coding sequence
                            out_list.append(mut_pos)
                        
                    if len(set(mut_nt_list)) == 1:
                        valid_seq+=var_count
                        
                        codon = int(list(set(mut_nt_list))[0])
                        
                        wt_seq = wt_codons[codon]
                                                                      
                        new_seq = [x for x in wt_seq]
                        
                        for mutation in mut_list:
                            mut_pos = int(mutation.split('|')[1])
                            mutation = mutation.split('|')[2]
                            
                            codon_pos = (mut_pos-1)%3
                            
                            new_seq[codon_pos] = mutation
                            
                        new_codon = ''.join(new_seq)
                        
                        if new_codon in list(variant_abundance_dict[codon].keys()):

                            variant_abundance_dict[codon][new_codon]+=var_count
                            
                        else:
                            variant_abundance_dict[codon][new_codon]=var_count
  
                    # This is for mutations outside of the coding sequence
                    elif len(set(mut_nt_list)) == 0 and len(out_list)>=1:
                        wt_count+=var_count
                    
        ## For low abundance variants            
        elif var_count < 20:

            for mutation in mut_list:

                if 'del' in mutation:
                    mut_info = mutation.split('del')
                    mut_pos = int(mut_info[0])

                    if mut_pos == 1:
                        mut_list.remove(mutation)
                  
            if len(mut_list) <=3 and 'ins' not in str(mut_list) and 'del' not in str(mut_list):

                if len(mut_list) ==0:
                    wt_count += var_count

                else:
                    #print(mut_list)
                    mut_nt_list = []

                    out_list = []

                    for mutation in mut_list:

                        mut_pos = int(mutation.split('|')[1])

                        if mut_pos >= frag_start and mut_pos <= frag_end:

                            mut_nt_list.append(codon_groups[mut_pos])

                        else:
                            out_list.append(mut_pos)

                    if len(set(mut_nt_list)) == 1:
                        valid_seq+=var_count

                        codon = int(list(set(mut_nt_list))[0])

                        wt_seq = wt_codons[codon]

                        new_seq = [x for x in wt_seq]

                        for mutation in mut_list:
                            mut_pos = int(mutation.split('|')[1])
                            mutation = mutation.split('|')[2]

                            codon_pos = (mut_pos-1)%3

                            new_seq[codon_pos] = mutation

                        new_codon = ''.join(new_seq)
                        
                        if new_codon in list(variant_abundance_dict[codon].keys()):
                            
                            variant_abundance_dict[codon][new_codon]+=var_count
  
                        else:
                            variant_abundance_dict[codon][new_codon]=var_count
                    
    return variant_abundance_dict, wt_count, wt_codons

codontable_standard = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

def convert_dict_to_aa(codon_dict):
    
    aa_dict_of_dicts = {}
    
    for pos in list(codon_dict.keys()):
        
        aa_dict_of_dicts[pos] = {}
        
        for codon in list(codon_dict[pos].keys()):            
            
            if np.isnan(codon_dict[pos][codon]):
                continue
            
            aa = codontable_standard[codon]
            
            if aa in list(aa_dict_of_dicts[pos].keys()):           
                aa_dict_of_dicts[pos][aa] += codon_dict[pos][codon]
                
            else:
                aa_dict_of_dicts[pos][aa] = codon_dict[pos][codon]
                
    return aa_dict_of_dicts

os.makedirs('./../../Data/Analysis_MiSeq/aggregate_dataframes/', exist_ok = True)
os.makedirs('./../../Data/Analysis_MiSeq/aggregate_dataframes/Codons/', exist_ok = True)
os.makedirs('./../../Data/Analysis_MiSeq/aggregate_dataframes/Residues/', exist_ok = True)

## Need to loop through the needle files and call the following functions on each
# Loop through the files
for i in range(0,16):
    needle_file = './../../Data/Analysis/amplicon_align/MiSeq_R67_bacteria_pool_' + str(i) + '.needle'
    
    # Initialize dictionaries
    align_count_dict = {}
    mut_count_dict = {}
    
    # Find mutations (this function parses the needle file)
    pool_1_muts = find_mutations(needle_file, 'MiSeq_R67_bacteria')
    
    ## Get variant count

    frag1_dict = get_variant_count_1(pool_1_muts, 'MiSeq_R67_bacteria', 145, 2, 78)
    
    wt_count = frag1_dict[1]
    
    # Extract the dictionary with the WT codons
    wt_codons = frag1_dict[2]
    
    # Add a loop to fill in the WT codons
    for position, wt_codon in wt_codons.items():
        frag1_dict[0][position][wt_codon] = wt_count
    
    ## Convert to the dataframes we need
    # Codons
    F1_KI_df = pd.DataFrame(frag1_dict[0])
    
    # Arrange codons
    codon_order = ['TAA', 'TAG', 'TGA', # *
              'GGA', 'GGC', 'GGG', 'GGT', # G
              'GCA', 'GCC', 'GCG', 'GCT', # A
              'GTA', 'GTC', 'GTG', 'GTT', # v
              'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG', # L
              'ATA', 'ATC', 'ATT', # I
              'ATG', # M
              'TGC', 'TGT', # C
              'CCA', 'CCC', 'CCG', 'CCT', # P
              'TGG', # W
              'TTC', 'TTT', # F
              'TAC', 'TAT', # Y
              'AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT',  # S
              'ACA', 'ACC', 'ACG', 'ACT', # T
              'AAC', 'AAT', # N
              'CAA', 'CAG', # Q
              'CAC', 'CAT', # H
              'AAA', 'AAG', # K
              'AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT', # R
              'GAC', 'GAT', # D
              'GAA', 'GAG' # E
              ]
    
    F1_KI_df.index = pd.CategoricalIndex(F1_KI_df.index, categories= codon_order)
    F1_KI_df.sort_index(level=0, inplace=True)
    
    ## Added this line to repace NAs with zeroes
    F1_KI_df = F1_KI_df.fillna(0)
    
    # Residues
    aa_frag_1 = convert_dict_to_aa(frag1_dict[0])
    F1_KI_aa_df = pd.DataFrame(aa_frag_1)
    
    # Save dataframes

    F1_KI_df.to_csv('./../../Data/Analysis_MiSeq/aggregate_dataframes/Codons/MiSeq_R67_bacteria_pool_' + str(i) + '_codons.txt', sep = '\t')
    F1_KI_aa_df.to_csv('./../../Data/Analysis_MiSeq/aggregate_dataframes/Residues/MiSeq_R67_bacteria_pool_' + str(i) + '_residues.txt', sep = '\t')
 
os.makedirs('./../../Data/Analysis_MiSeq/read_abundances/', exist_ok = True)
os.makedirs('./../../Data/Analysis_MiSeq/read_abundances/Codons/', exist_ok = True)
os.makedirs('./../../Data/Analysis_MiSeq/read_abundances/Residues/', exist_ok = True)


def get_timepoint_fraction_df(needle_file, Sample_ID):
    
    # Find mutations (this function parses the needle file)
    pool_1_muts = find_mutations(needle_file, 'MiSeq_R67_bacteria')
    
    frag_1_dict = get_variant_count_1(pool_1_muts, 'MiSeq_R67_bacteria', 145, 2, 78)
    mut_dict = frag_1_dict[0]
    wt_count = frag_1_dict[1]
        
    variant_df = pd.DataFrame(mut_dict)
    
    array_size=len(variant_df.to_numpy().flatten())
    
    # Calculate the total number of reads before adding the WT count to the matrix
    read_total = variant_df.sum().sum()+wt_count+array_size
    print(Sample_ID, read_total, wt_count, array_size)
    
    # Extract the dictionary with the WT codons
    wt_codons = frag1_dict[2]
    
    # Add a loop to fill in the WT codons
    for position, wt_codon in wt_codons.items():
        mut_dict[position][wt_codon] = wt_count
    
    # Update variant_df
    variant_df = pd.DataFrame(mut_dict)
    variant_df_no_NaN = variant_df.fillna(0)
    
    variant_df_no_NaN = variant_df_no_NaN + 1
        
    read_fraction_df = variant_df_no_NaN/(read_total)
    print(read_total, 1/read_total)
    
    read_fraction_df.rename_axis(read_total, inplace=True)
    
    df_out_path='./../../Data/Analysis_MiSeq/read_abundances/Codons/'+str(Sample_ID)+'_read_frac.csv'
    
    read_fraction_df.to_csv(df_out_path, sep=',')
    
    # Convert dataframe to the residue level
    aa_frag_1 = convert_dict_to_aa(mut_dict)
    variant_df_aa = pd.DataFrame(aa_frag_1)
    
    # Get proportion of reads at the residue level
    array_size=len(variant_df_aa.to_numpy().flatten())
    read_total = variant_df_aa.sum().sum() + wt_count + array_size
    
    print(Sample_ID, read_total, wt_count, array_size)
    
    variant_df_aa_no_NaN = variant_df_aa.fillna(0)
    variant_df_aa_no_NaN = variant_df_aa_no_NaN + 1
    
    # Normalize by the total of reads
    read_fraction_df_aa = variant_df_aa_no_NaN/(read_total)
    print(read_total, 1/read_total)
    
    read_fraction_df_aa.rename_axis(read_total, inplace=True)
    
    # Save file
    df_out_path='./../../Data/Analysis_MiSeq/read_abundances/Residues/'+str(Sample_ID)+'_read_frac.csv'
    read_fraction_df_aa.to_csv(df_out_path, sep=',')
    
    # variant_df has the raw counts for each codon
    # read_fraction_df has counts normalized by the total of reads
    return read_fraction_df, variant_df

# The loop to calculate read proportions in each sample
for i in range(0,16):
    needle_file = './../Analysis_wt/amplicon_align/MiSeq_R67_bacteria_pool_' + str(i) + '.needle'
    
    read_fraction_df, variant_df = get_timepoint_fraction_df(needle_file, i)
    
    print('----------------')


# ## Work with the NovaSeq data

# ### Check sequencing quality and select for seqs w/ correct length
experiment_name = 'R67_DMS'

path_to_R1_file = '../../Data/Sequencing_data/NovaSeq_Feb2022_DfrB1_1P.fastq.gz'
path_to_R2_file = '../../Data/Sequencing_data/NovaSeq_Feb2022_DfrB1_2P.fastq.gz'

## FastQC quality control
subprocess.check_output(fastqc_path + ' ' + path_to_R1_file, shell=True)
subprocess.check_output(fastqc_path + ' ' +path_to_R2_file, shell=True)

# Initialize intermediate folders for the analysis
amplicon_sequences_path = '../../Data/Analysis_NovaSeq/amplicon_sequences'
os.makedirs(amplicon_sequences_path, exist_ok = True)

bowtie_indexes_path = '../../Data/Analysis_NovaSeq/bowtie_indexes'
os.makedirs(bowtie_indexes_path, exist_ok = True)

temp_path = '../../Data/Analysis_NovaSeq/temp'
os.makedirs(temp_path, exist_ok = True)

## Trim MiSeq reads with Trimmomatic
trim_call = 'java -jar ' + trimmomatic_path + ' PE -threads 6 -trimlog log.txt '
trim_call += path_to_R1_file+' '+path_to_R2_file+' '
trim_call += '-baseout ../../Data/Analysis_NovaSeq/temp/minlen250.fastq MINLEN:250 CROP:225'

trim_call

subprocess.check_output(trim_call, shell=True)

## FastQC quality control on the trimmed reads
subprocess.check_output('fastqc '+'../../Data/Analysis_NovaSeq/temp/minlen250_1P.fastq', shell=True)
subprocess.check_output('fastqc '+'../../Data/Analysis_NovaSeq/temp/minlen250_2P.fastq', shell=True)


# ### check abundance of libraries
# Output files
plate_pcr_amplicon_for = os.path.join(amplicon_sequences_path, 'plate_pcr_for.fa')
plate_pcr_amplicon_rev = os.path.join(amplicon_sequences_path, 'plate_pcr_rev.fa')

# Load the metadata to extract the plate indexes for the NovaSeq data
sheet = '../../Data/TableS2_DMS_sample_description.tsv'
sample_sheet = pd.read_csv(sheet, header=0, delimiter = '\t')

novaseq_sheet = sample_sheet[sample_sheet['Sequencing platform'] == 'NovaSeq']

for_plate_indexes = novaseq_sheet['PE1_index'].unique()
rev_plate_indexes = novaseq_sheet['PE2_index'].unique()

# Valid plates for the NovaSeq data 
valid_plates = list(zip(novaseq_sheet.PE1_index, novaseq_sheet.PE2_index))

print(for_plate_indexes)
print(rev_plate_indexes)
print(valid_plates)

# Make a bowtie index of the forward plate indices
with open(plate_pcr_amplicon_for, 'w') as fasta_out:
    
    for index in for_plate_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = degen_seq+index_seq+universal_seqs['plate_fivep_sticky']+'\n'
        fasta_out.write(amplicon_frag)
        
        
for_plate_index = os.path.join(bowtie_indexes_path, experiment_name+'_for_plate')

for_plate_bowtie_index_call = bowtie_build_path + ' -f -r -o 4 ' +plate_pcr_amplicon_for+' '+for_plate_index
for_plate_indexing_log = subprocess.check_output(for_plate_bowtie_index_call, shell=True)
        
# Make a bowtie index of the forward plate indices    
with open(plate_pcr_amplicon_rev, 'w') as fasta_out:
    
    for index in rev_plate_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = degen_seq + index_seq+ reverse_complement(universal_seqs['plate_threep_sticky'])+'\n'
        fasta_out.write(amplicon_frag)    
    
rev_plate_index = os.path.join(bowtie_indexes_path, experiment_name+'_rev_plate')
rev_plate_bowtie_index_call = bowtie_build_path + ' -f -r -o 4 '+plate_pcr_amplicon_rev+' '+rev_plate_index
rev_plate_indexing_log = subprocess.check_output(rev_plate_bowtie_index_call, shell=True)


# ## Align for barcode

## Bowtie arguments:
# -t print time
# -v maximum number of mismatches tolerated
# -p number of threads
# -k number of valid alignments to report
# --trim3 trim x nt from 3' end
# --trim5 trim x nt from 5' end
# --norc do not align to reverse complement
# --chunkmbs Mbs of memory for each thread

test_align_call = bowtie_path + ' -t -v 3 -p 6 -k 1 --trim3 250 --trim5 5 --norc --chunkmbs 256 '

test_align_call += for_plate_index+' '
test_align_call += os.path.join(temp_path, 'minlen250_1P.fastq') + ' '
test_align_call += os.path.join(temp_path, 'plate_for_align.txt')

print (test_align_call)

subprocess.check_output(test_align_call, shell = True)

test_align_call = bowtie_path + ' -t -v 3 -p 6 -k 1 --trim3 250 --trim5 5 --norc --chunkmbs 256 '

test_align_call += rev_plate_index+' '
test_align_call += os.path.join(temp_path, 'minlen250_2P.fastq') + ' '
test_align_call += os.path.join(temp_path, 'plate_rev_align.txt')
print (test_align_call)
subprocess.check_output(test_align_call, shell = True)

plate_align_for_output = os.path.join(temp_path, 'plate_for_align.txt')
plate_align_rev_output = os.path.join(temp_path, 'plate_rev_align.txt')


plate_index_found_for = {}
plate_index_found_rev = {}
# empty containers that will hold    read_name: index    data pairs for the forward and reverse reads that
# had valid alignments

with open(plate_align_for_output, 'r') as for_index_positives:
    for line in for_index_positives:
        # opens P1 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info, not the illumina annotations in the read 
        # name. See http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/
        #           Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm
        # For info on illumina read name convention
        plate_index_found_for[read_name] = line[2]
        # add read P1 index to dict

with open(plate_align_rev_output, 'r') as rev_index_positives:
    for line in rev_index_positives:
        # opens P2 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info
        plate_index_found_rev[read_name] = line[2]
        # add read P2 index to dict


plate_both = list(set(plate_index_found_for.keys()) & set(plate_index_found_rev.keys()))

print(len(plate_both))


# This block looks for the different indexes in the data

RC_pcr_amplicon_for = os.path.join(amplicon_sequences_path, 'RC_pcr_for.fa')
RC_pcr_amplicon_rev = os.path.join(amplicon_sequences_path, 'RC_pcr_rev.fa')

RC_for_indexes = range(1,16)

RC_rev_indexes = range(1,16)

print(RC_for_indexes)
print(RC_rev_indexes)


with open(RC_pcr_amplicon_for, 'w') as fasta_out:
    
    for index in RC_for_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = universal_seqs['plate_fivep_sticky']+index_seq+universal_seqs['RC_fivep_sticky']+'\n'
        fasta_out.write(amplicon_frag)
        
        
for_RC_index = os.path.join(bowtie_indexes_path, experiment_name+'_for_RC')
for_RC_bowtie_index_call = bowtie_build_path + ' -f -r -o 4 '+RC_pcr_amplicon_for+' '+for_RC_index
for_RC_indexing_log = subprocess.check_output(for_RC_bowtie_index_call, shell=True)
        
        
with open(RC_pcr_amplicon_rev, 'w') as fasta_out:
    
    for index in RC_rev_indexes:

        index_seq = indexes[str(index)]
        
        header = '>'+str(index)+'\n'
        fasta_out.write(header)
        
        amplicon_frag = reverse_complement(universal_seqs['plate_threep_sticky']) + index_seq
        amplicon_frag += reverse_complement(universal_seqs['RC_threep_sticky'])+'\n'
        fasta_out.write(amplicon_frag)    
    
rev_RC_index = os.path.join(bowtie_indexes_path, experiment_name+'_rev_RC')
rev_RC_bowtie_index_call = bowtie_build_path + ' -f -r -o 4 '+RC_pcr_amplicon_rev+' '+rev_RC_index
rev_RC_indexing_log = subprocess.check_output(rev_RC_bowtie_index_call, shell=True)

test_align_call = bowtie_path + ' -t -v 3 -p 6 -k 1 --trim3 215 --trim5 28 --norc --chunkmbs 256 '

test_align_call += for_RC_index+' '
test_align_call += os.path.join(temp_path, 'minlen250_1P.fastq') + ' '
test_align_call += os.path.join(temp_path, 'RC_for_align.txt')

subprocess.check_output(test_align_call, shell = True)

test_align_call = bowtie_path + ' -t -v 3 -p 6 -k 1 --trim3 217 --trim5 26 --norc --chunkmbs 256 '

test_align_call += rev_RC_index+' '
test_align_call += os.path.join(temp_path, 'minlen250_2P.fastq') + ' '
test_align_call += os.path.join(temp_path, 'RC_rev_align.txt')

subprocess.check_output(test_align_call, shell = True)

RC_align_for_output = os.path.join(temp_path, 'RC_for_align.txt')
RC_align_rev_output = os.path.join(temp_path, 'RC_rev_align.txt')


index_found_for_RC = {}
index_found_rev_RC = {}
# empty containers that will hold    read_name: index    data pairs for the forward and reverse reads that
# had valid alignments


with open(RC_align_for_output, 'r') as for_index_positives:
    for line in for_index_positives:
        # opens P1 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info, not the illumina annotations in the read 
        # name. See http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/
        #           Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm
        # For info on illumina read name convention
        index_found_for_RC[read_name] = line[2]
        # add read P1 index to dict


with open(RC_align_rev_output, 'r') as rev_index_positives:
    for line in rev_index_positives:
        # opens P2 index alignment pack file and loops through lines

        line = line.split('\t')
        # split line to extract data
        read_name = line[0].split(' ')[0]
        # Only keep the cluster coordinates and run associated info
        index_found_rev_RC[read_name] = line[2]
        # add read P2 index to dict


both_indexes_RC = list(set(index_found_for_RC.keys()) & set(index_found_rev_RC.keys()))

print(len(both_indexes_RC))

def plot_reads_in_plate(plate_both, RC_both, plate_index_found_for, plate_index_found_rev, index_found_for_RC, index_found_rev_RC, plate,
                       silent=True):
    
    plate_rows = sorted(set(index_found_for_RC.values()))
    plate_columns = sorted(set(index_found_rev_RC.values()))
    
    read_plate_RC_dict = {}
    grid_dict = {}
    
    for row in plate_rows:
        
        grid_dict[int(row)] = {}
        
        for column in plate_columns:
            
            grid_dict[int(row)][int(column)] = 0
            
    if silent == False:
        
        print(grid_dict)
    
    plate_and_rc_both = list((set(plate_both)&set(both_indexes_RC)))
    
    for read in plate_and_rc_both:
        
        for_index_plate = int(plate_index_found_for[read])
        rev_index_plate = int(plate_index_found_rev[read])
        
        if tuple((int(for_index_plate), int(rev_index_plate))) == plate:
                        
            for_index_RC = int(index_found_for_RC[read])
            rev_index_RC = int(index_found_rev_RC[read])
            
            grid_dict[for_index_RC][rev_index_RC] += 1
            
            read_plate_RC_dict[read] = [int(for_index_plate), int(rev_index_plate), for_index_RC, rev_index_RC]
            
            
            
    plate_df = pd.DataFrame.from_dict(grid_dict, orient='index')
    return plate_df, read_plate_RC_dict
            
            
# Keep in mind I am only using the first one of the valid plates array
# since it's the only one I need for this analysis
F_1 = plot_reads_in_plate(plate_both, both_indexes_RC, plate_index_found_for, plate_index_found_rev, index_found_for_RC, 
                                index_found_rev_RC, valid_plates[0], silent = False)

F_1[0].head(10)


col_order = [1,2,3,4,5,6,7,8,9,10]
plt.figure(figsize=(14/1.5, 10/1.5))


sns.heatmap(F_1[0][col_order], vmax=6400, cmap='viridis',)
plt.xlabel('Column Barcode', fontsize = 16)
plt.ylabel('Row Barcode', fontsize = 16)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.axhline(2, color='white')
plt.axhline(4, color='white')
plt.axhline(6, color='white')
plt.axhline(8, color='white')
plt.title('Fragment 1')


all_plate_index_dicts = []

for plate in valid_plates:
    
    mapping = plot_reads_in_plate(plate_both, both_indexes_RC, plate_index_found_for, plate_index_found_rev, index_found_for_RC, 
                                index_found_rev_RC, plate, silent=True)[1]
    
    all_plate_index_dicts.append(mapping)

demux_reads_path = '../../Data/Analysis_NovaSeq/demultiplexed_reads'
os.makedirs(demux_reads_path, exist_ok = True)


# ## Look for the row column barcodes
# Organize the combinations of plate and row column barcodes in a dictionary
plate_RC_combi_to_dict = {}

out_fastq_path_dict_for_df = {}

out_fastq_path_dict_demult = {}

dir_path = demux_reads_path

plate_RC_combi = list(zip(novaseq_sheet.Sample,
    novaseq_sheet.PE1_index, novaseq_sheet.PE2_index,
    novaseq_sheet.RC_for_index, novaseq_sheet.RC_rev_index
))

for sample_data in plate_RC_combi:
    sample = 'pool_' + str(sample_data[0])
    
    plate_for = sample_data[1]
    plate_rev = sample_data[2]
    RC_for = sample_data[3]
    RC_rev = sample_data[4]
    
    plate_RC_combi_to_dict[new_key] = [plate_for, plate_rev, RC_for, RC_rev]

    out_path = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)

    out_path_1 = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)+'_for'+'.fastq'
    out_path_2 = dir_path+str(sample)+'_indexes_'+str(plate_for)+'_'+str(plate_rev)+'_'+str(RC_for)+'_'+str(RC_rev)+'_rev'+'.fastq'

    out_fastq_path_dict_for_df[sample] = out_path

    # Saves the paths to the files so that they can be used in the following function
    out_fastq_path_dict_demult[tuple([plate_for, plate_rev, RC_for, RC_rev])] = out_path

    # Create empty files to append
    create_file_for = 'touch '+ out_path_1
    create_file_rev = 'touch '+ out_path_2

    subprocess.check_output(create_file_for, shell=True)
    subprocess.check_output(create_file_rev, shell=True)
    

for plate_dict in all_plate_index_dicts:
    
    print(len(plate_dict.keys()))

    for read in list(plate_dict.keys()):

        index_combi = plate_dict[read]

        if index_combi not in plate_RC_combi_to_dict.values():

            plate_dict.pop(read, None)
        
    print(len(plate_dict.keys()))

## This function will add the reads to the files we created above
def split_reads_from_fastq(file, file_type, plate_mapping_list):
    
    to_write = 'none'
    read_name = 'none'
    
    
    with open(file, 'r') as source_fastq:
        
        for line in source_fastq:

            # If this line indicates the name of the read
            if line.startswith('@M') == True:
                
                if to_write == 'none':
                    
                    read_name = line.split(' ')[0].strip('@')
                    to_write = line

                # Otherwise, we are looking at the information about the read
                else:
                    
                    # Check in which plate it is
                    for plate_read_dict in plate_mapping_list:
                        
                                                
                        if read_name in plate_read_dict.keys():

                            indexes = tuple(plate_read_dict[read_name])
                            
                            # This is where we distinguish between the files for forward and reverse reads
                            if file_type == 'forward':
                                
                                file_path = out_fastq_path_dict_demult[indexes]+'_for.fastq'
                                
                            elif file_type == 'reverse':
                                
                                file_path = out_fastq_path_dict_demult[indexes]+'_rev.fastq'

                            
                            
                            with open(file_path, 'a') as dest:
                                
                                   dest.write(to_write)

                            break





                    read_name = line.split(' ')[0].strip('@')

                    to_write = line
                

            else:

                to_write += line

                
# Separate the forward from the reverse reads
split_reads_from_fastq(os.path.join(temp_path, 'minlen250_1P.fastq'), 'forward', all_plate_index_dicts)
split_reads_from_fastq(os.path.join(temp_path, 'minlen250_2P.fastq'), 'reverse', all_plate_index_dicts)


# ## Load the reference amplicon
path_to_amplicons = './../../Data/uni_amplicon_sequences/miseq_r67_bacteria.fasta'

amplicon_info_dict = get_dict_of_seq(path_to_amplicons)

amplicon_length_dict = {}

amplicon_dict = {}

for amplicon in amplicon_info_dict:
    
    amplicon_name = amplicon.split('|')[0]

    amplicon_dict[amplicon_name] = amplicon_info_dict[amplicon]
    
    amplicon_fasta_path = './../../Data/Analysis_NovaSeq/amplicon_sequences/'+amplicon_name+'.fasta'
    
    with open(amplicon_fasta_path, 'w') as dest:
        
        seq_ID = '>'+amplicon_name+'\n'
        seq = amplicon_dict[amplicon_name]+'\n'
        
        dest.write(seq_ID)
        dest.write(seq)
    

print (amplicon_info_dict.keys())
print (amplicon_dict)


# ## Align to the reference amplicon
os.mkdir('./../../Data/Analysis_NovaSeq/amplicon_align')

## Define a function to do the alignments with needle
def needle_align_on_ref(ref_orf, filepath):
               
    ref_seq = amplicon_dict[ref_orf]
    
    ref_fasta_path = './../../Data/Analysis_NovaSeq/amplicon_sequences/'+ref_orf+'.fasta '
    
    # Extract the name of the pool
    pool_number = re.search('pool_\d+', filepath).group(0)
    
    ## Needle arguments:
    # -auto Turn of prompts
    # -gapopen penalty for opening a gap
    # -asequence one of the two input sequences to align
    # -bsequence the second input sequence to align
    # -aformat output format for the alignment (markx10 is an easily parseable format http://emboss.sourceforge.net/docs/themes/alnformats/align.markx10 )
    # -outfile path to the output file
    
    # For the analysis_wt folder
    needle_out = './../../Data/Analysis_NovaSeq/amplicon_align/'+ ref_orf + '_' + pool_number +'.needle'
    
    needle_call = 'needle -auto -gapopen 50 -asequence '+ ref_fasta_path
    
    needle_call += '-bsequence '+ filepath +' -aformat3 markx10 -outfile '+needle_out
    
    print(needle_call)
    
    subprocess.check_output(needle_call, shell = True)
          
    return

aggregate_list = glob.glob('../../Data/Analysis_NovaSeq/aggregate_reads/*fasta')

## Loop through the files of aggregate reads
for filepath in aggregate_list:
    needle_align_on_ref('MiSeq_R67_bacteria', filepath)


# ## Parse needle output and generate the dataframes with read counts
def parse_needle_output(needle_align_path):
    
    n_aligns = 0
    align_seqs_dict = {}       
   
    with open(needle_align_path, 'r') as source:

        current_align = ''
        current_qseq = ''
        current_sseq = ''
        qseq_done = 0

        for line in source:

            if line.startswith('>>>') == True:

                n_aligns +=1
                align_name = line.strip('>>>')

                if n_aligns != 1:

                    align_seqs_dict[current_align] = [current_qseq, current_sseq]
                    current_align = align_name
                    current_qseq = ''
                    current_sseq = ''
                    qseq_done = 0

                else:

                    current_align = align_name

            elif line.startswith(';') == False and line.startswith('>') == False and line.startswith('\n') == False and line.startswith('#') == False:

                if qseq_done == 1:
                    current_sseq += line.strip('\n').upper()

                else:
                    current_qseq += line.strip('\n').upper()

            elif line.startswith('#--') == True:

                align_seqs_dict[align_name] = [current_qseq, current_sseq]

            else:

                if qseq_done == 0 and current_qseq != '':
                    qseq_done =1            
                             
    return align_seqs_dict, n_aligns

def find_mutations(path, ref_orf):
    
    allele_dict = {}

    align_dict, align_count = parse_needle_output(path)

    for entry in list(align_dict.keys()):

        read_var_list = []

        query_seq = align_dict[entry][1]
        # aligned prot sequence of the strain

        align_ref = align_dict[entry][0]
        # aligned prot sequence of the reference

        gap_adjust = 0
        # value used to adjust the protein sequence index for the presence of insertions in the strain sequence vs the 
        # reference strain

        backtrack_adjust = 0

        temp_var = None
        # temporary variable to hold the sequence of an insertion or deletion as a string. When the gap ends, annotation 
        # will be added to strain_var_list

        indel_start = 0
        # position start of the indel annotation in the reference sequence, with adjustment for gap presence

        ref_seq_no_gaps = align_ref.replace('-','')
        # Make a copy of the reference sequence without gaps 
        
        # align_start = (amplicon_dict[ref_orf].index(ref_seq_no_gaps))+1
        align_start = (amplicon_dict[ref_orf].upper().index(ref_seq_no_gaps))+1
        # Look for the starting position of the gapless part of the reference sequence
        # This helps remove gaps at the start and the end of the alignment
        
        query_seq_no_gaps = len(query_seq.replace('-',''))
        # Make a copy of the query sequence without gaps
        
        for nt in range(0, len(align_ref)):
            # iterates through the entire alignment of the strain prot sequence

            if query_seq[nt] == '-':
                # detect a deletion variant

                # logic for indel detection/annotation:
                #
                # suppose we have this alignment  
                #
                # 1 2 3 4 5 6 7 8 9
                # A T - - A A A T G    strain variant: del gaps are indexed because the aa index is based on reference
                # A T K P A - - T G
                # 1 2 3 4 5     6 7    reference: insert gaps not indexed because aa positions do (actually don't?) exist in reference
                #
                # following this logic, every time an insertion is detected and annotated, the gap_adjust value is 
                # incremented by the length of the gap and used to adjust the variant mapping to make it match the 
                # reference index values. The indel aa postion is the first residue detected as part of the indel


                if indel_start == 0:
                    # checks if the character is the start or the continuation of a gap in the alignment

                    temp_var = 'del'+ align_ref[nt]
                    indel_start = (nt+1-gap_adjust)
                    # if it is, starts a new annotation entry with a start position compensated for previous insertions
                    # (if any)

                    backtrack_adjust += 1

                else:

                    temp_var += align_ref[nt]
                    # if it is not, adds the following aa to the deletion annotation

                    backtrack_adjust += 1


            elif align_ref[nt] == '-':
                # detects an insertion variant

                if indel_start == 0:
                    # checks if the character is the start or the continuation of a gap in the alignment

                    temp_var = 'ins'+ query_seq[nt]

                    indel_start = (nt+1-gap_adjust)
                    # if it is, starts a new annotation entry with a start position compensated for previous insertions
                    # (if any)

                    gap_adjust += 1
                    # increments the gap adjust for the this added aa in the strain sequence                   

                else:

                    temp_var += query_seq[nt]
                    # if it is not, adds the following aa to the insertion annotation

                    gap_adjust += 1
                    # increments the gap adjust for the this added aa in the strain sequence

            elif query_seq[nt] != align_ref[nt]:
                # detects a mismatch between the strain sequence and the reference

                variant = align_ref[nt]+'|'+str((nt+1-gap_adjust))+'|'+query_seq[nt]
                read_var_list.append(variant)
                # creates an annotation for the strain-reference aa mismatch and appends it to the list of 
                # annotations

            else:

                 if indel_start != 0:
                    # detects if there is currently an open gap entry. If there is, then the detected mismatch means 
                    # that it has now concluded

                    read_var_list.append(str((indel_start))+temp_var)
                    temp_var = None
                    indel_start = 0
                    # adds the indel annotation to the strain variant list and resets temporary variables for the next 
                    # indel entry

        if len(read_var_list)<25:  
            allele_dict[entry] = read_var_list, align_start
                           
    return allele_dict  

def get_variant_count_1(mutation_set, ref_seq, frag_start, codon_start, n_aa):
    
    variant_abundance_dict ={}
    variants = list(mutation_set.keys())
    codon_groups = {}
    codon = 1
    wt_count =0
    valid_seq=0
    
    # Calculate the end of the fragment with the frag start and the number of residues
    frag_end = frag_start + ((n_aa -1 )* 3) - 1
    
    for nt in range(0, (n_aa - 1)*3):
        pos = nt + frag_start
        
        if nt % 3 == 0:
            codon += 1
            
        codon_groups[pos] = codon
        
        variant_abundance_dict[codon] = {}       
        
    wt_codons = {}
        
    ref = amplicon_dict[ref_seq].upper()
    
    # Set the abundance of WT codons to not a number (nan) since they would be overrepresented 
    for aa in range(0, n_aa - 1):

        offset = frag_start - 1
        start = offset+(aa*3)
        wt_codon=ref[start:(start+3)]
        wt_codons[(aa+codon_start)] = wt_codon
        variant_abundance_dict[aa+codon_start][wt_codon]=np.nan 
    
    for variant in variants:
        
        var_info = variant.split(',')
        var_count =int(var_info[1].split(';')[1].strip('size='))      
        mut_list = mutation_set[variant][0]
        filtered_list = []
        
        # Only keep variants that appeared in more than 20 reads
        if var_count>=20:
        
            for mutation in mut_list:

                if 'del' in mutation:
                    mut_info = mutation.split('del')
                    mut_pos = int(mut_info[0])
                    
                    if mut_pos == 1:
                        mut_list.remove(mutation)
            
            # If not an indel
            if 'ins' not in str(mut_list) and 'del' not in str(mut_list):
                
                if len(mut_list) ==0:
                    wt_count += var_count
                    
                else:
                    mut_nt_list = []
                    out_list = []
                    
                    for mutation in mut_list:
                        
                        mut_pos = int(mutation.split('|')[1])
                        
                        if mut_pos >= frag_start and mut_pos <= frag_end:
                        
                            # mut_nt_list contains the list of mutated positions inside the coding sequence
                            mut_nt_list.append(codon_groups[mut_pos])
                            
                        else:
                            # out_list contains the list of mutated positions outside the coding sequence
                            out_list.append(mut_pos)
                        
                    if len(set(mut_nt_list)) == 1:
                        valid_seq+=var_count
                        
                        codon = int(list(set(mut_nt_list))[0])
                        
                        wt_seq = wt_codons[codon]
                                                                      
                        new_seq = [x for x in wt_seq]
                        
                        for mutation in mut_list:
                            mut_pos = int(mutation.split('|')[1])
                            mutation = mutation.split('|')[2]
                            
                            codon_pos = (mut_pos-1)%3
                            
                            new_seq[codon_pos] = mutation
                            
                        new_codon = ''.join(new_seq)
                        
                        if new_codon in list(variant_abundance_dict[codon].keys()):

                            variant_abundance_dict[codon][new_codon]+=var_count
                            
                        else:
                            variant_abundance_dict[codon][new_codon]=var_count
  
                    # This is for mutations outside of the coding sequence
                    elif len(set(mut_nt_list)) == 0 and len(out_list)>=1:
                        wt_count+=var_count
                    
        ## For low abundance variants            
        elif var_count < 20:

            for mutation in mut_list:

                if 'del' in mutation:
                    mut_info = mutation.split('del')
                    mut_pos = int(mut_info[0])

                    if mut_pos == 1:
                        mut_list.remove(mutation)
                  
            if len(mut_list) <=3 and 'ins' not in str(mut_list) and 'del' not in str(mut_list):

                if len(mut_list) ==0:
                    wt_count += var_count

                else:
                    #print(mut_list)
                    mut_nt_list = []

                    out_list = []

                    for mutation in mut_list:

                        mut_pos = int(mutation.split('|')[1])

                        if mut_pos >= frag_start and mut_pos <= frag_end:

                            mut_nt_list.append(codon_groups[mut_pos])

                        else:
                            out_list.append(mut_pos)

                    if len(set(mut_nt_list)) == 1:
                        valid_seq+=var_count

                        codon = int(list(set(mut_nt_list))[0])

                        wt_seq = wt_codons[codon]

                        new_seq = [x for x in wt_seq]

                        for mutation in mut_list:
                            mut_pos = int(mutation.split('|')[1])
                            mutation = mutation.split('|')[2]

                            codon_pos = (mut_pos-1)%3

                            new_seq[codon_pos] = mutation

                        new_codon = ''.join(new_seq)
                        
                        if new_codon in list(variant_abundance_dict[codon].keys()):
                            
                            variant_abundance_dict[codon][new_codon]+=var_count
  
                        else:
                            variant_abundance_dict[codon][new_codon]=var_count
                    
    return variant_abundance_dict, wt_count, wt_codons

codontable_standard = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

def convert_dict_to_aa(codon_dict):
    
    aa_dict_of_dicts = {}
    
    for pos in list(codon_dict.keys()):
        
        aa_dict_of_dicts[pos] = {}
        
        for codon in list(codon_dict[pos].keys()):            
            
            if np.isnan(codon_dict[pos][codon]):
                continue
            
            aa = codontable_standard[codon]
            
            if aa in list(aa_dict_of_dicts[pos].keys()):           
                aa_dict_of_dicts[pos][aa] += codon_dict[pos][codon]
                
            else:
                aa_dict_of_dicts[pos][aa] = codon_dict[pos][codon]
                
    return aa_dict_of_dicts

os.makedirs('./../../Data/Analysis_NovaSeq/aggregate_dataframes/', exist_ok = True)
os.makedirs('./../../Data/Analysis_NovaSeq/aggregate_dataframes/Codons/', exist_ok = True)
os.makedirs('./../../Data/Analysis_NovaSeq/aggregate_dataframes/Residues/', exist_ok = True)

## Need to loop through the needle files and call the following functions on each
list_files = glob.glob('./../../Data/Analysis_NovaSeq/amplicon_align/*')

for needle_file in list_files:

    # Get the sample number
    i = os.path.basename(needle_file)[0:-7].split('_')[-1]
    
    # Initialize dictionaries
    align_count_dict = {}
    mut_count_dict = {}
    
    # Find mutations (this function parses the needle file)
    pool_1_muts = find_mutations(needle_file, 'MiSeq_R67_bacteria')
    
    ## Get variant count
    frag1_dict = get_variant_count_1(pool_1_muts, 'MiSeq_R67_bacteria', 145, 2, 78)
    
    wt_count = frag1_dict[1]
    
    # Extract the dictionary with the WT codons
    wt_codons = frag1_dict[2]
    
    # Add a loop to fill in the WT codons
    for position, wt_codon in wt_codons.items():
        frag1_dict[0][position][wt_codon] = wt_count
    
    ## Convert to the dataframes we need
    # Codons
    F1_KI_df = pd.DataFrame(frag1_dict[0])
    
    # Arrange codons
    codon_order = ['TAA', 'TAG', 'TGA', # *
              'GGA', 'GGC', 'GGG', 'GGT', # G
              'GCA', 'GCC', 'GCG', 'GCT', # A
              'GTA', 'GTC', 'GTG', 'GTT', # v
              'CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG', # L
              'ATA', 'ATC', 'ATT', # I
              'ATG', # M
              'TGC', 'TGT', # C
              'CCA', 'CCC', 'CCG', 'CCT', # P
              'TGG', # W
              'TTC', 'TTT', # F
              'TAC', 'TAT', # Y
              'AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT',  # S
              'ACA', 'ACC', 'ACG', 'ACT', # T
              'AAC', 'AAT', # N
              'CAA', 'CAG', # Q
              'CAC', 'CAT', # H
              'AAA', 'AAG', # K
              'AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT', # R
              'GAC', 'GAT', # D
              'GAA', 'GAG' # E
              ]
    
    F1_KI_df.index = pd.CategoricalIndex(F1_KI_df.index, categories= codon_order)
    F1_KI_df.sort_index(level=0, inplace=True)
    
    ## Added this line to repace NAs with zeroes
    F1_KI_df = F1_KI_df.fillna(0)
    
    # Residues
    aa_frag_1 = convert_dict_to_aa(frag1_dict[0])
    F1_KI_aa_df = pd.DataFrame(aa_frag_1)
    
    # Save dataframes
    F1_KI_df.to_csv('./../../Data/Analysis_NovaSeq/aggregate_dataframes/Codons/MiSeq_R67_bacteria_' + str(i) + '_codons.txt', sep = '\t')
    F1_KI_aa_df.to_csv('./../../Data/Analysis_NovaSeq/aggregate_dataframes/Residues/MiSeq_R67_bacteria_' + str(i) + '_residues.txt', sep = '\t')

def get_timepoint_fraction_df(needle_file, Sample_ID):
    
    # Find mutations (this function parses the needle file)
    pool_1_muts = find_mutations(needle_file, 'MiSeq_R67_bacteria')
    
    frag_1_dict = get_variant_count_1(pool_1_muts, 'MiSeq_R67_bacteria', 145, 2, 78)
    mut_dict = frag_1_dict[0]
    wt_count = frag_1_dict[1]
        
    variant_df = pd.DataFrame(mut_dict)
    
    array_size=len(variant_df.to_numpy().flatten())
    
    # Calculate the total number of reads before adding the WT count to the matrix
    read_total = variant_df.sum().sum()+wt_count+array_size
    print(Sample_ID, read_total, wt_count, array_size)
    
    # Extract the dictionary with the WT codons
    wt_codons = frag1_dict[2]
    
    # Add a loop to fill in the WT codons
    for position, wt_codon in wt_codons.items():
        mut_dict[position][wt_codon] = wt_count
    
    # Update variant_df
    variant_df = pd.DataFrame(mut_dict)
    variant_df_no_NaN = variant_df.fillna(0)
    
    variant_df_no_NaN = variant_df_no_NaN + 1
    
    #wt_fraction=wt_count/read_total
    
    read_fraction_df = variant_df_no_NaN/(read_total)
    print(read_total, 1/read_total)
    
    read_fraction_df.rename_axis(read_total, inplace=True)
    
    df_out_path='./../../Data/Analysis_NovaSeq/read_abundances/Codons/'+str(Sample_ID)+'_read_frac.csv'
    
    read_fraction_df.to_csv(df_out_path, sep=',')
    
    
    # Convert dataframe to the residue level
    aa_frag_1 = convert_dict_to_aa(mut_dict)
    variant_df_aa = pd.DataFrame(aa_frag_1)
    
    # Get proportion of reads at the residue level
    array_size=len(variant_df_aa.to_numpy().flatten())
    read_total = variant_df_aa.sum().sum() + wt_count + array_size
    
    print(Sample_ID, read_total, wt_count, array_size)
    
    variant_df_aa_no_NaN = variant_df_aa.fillna(0)
    variant_df_aa_no_NaN = variant_df_aa_no_NaN + 1
    
    # Normalize by the total of reads
    read_fraction_df_aa = variant_df_aa_no_NaN/(read_total)
    print(read_total, 1/read_total)
    
    read_fraction_df_aa.rename_axis(read_total, inplace=True)
    
    # Save file
    df_out_path='./../../Data/Analysis_NovaSeq/read_abundances/Residues/'+str(Sample_ID)+'_read_frac.csv'
    read_fraction_df_aa.to_csv(df_out_path, sep=',')
    
    # variant_df has the raw counts for each codon
    # read_fraction_df has counts normalized by the total of reads
    return read_fraction_df, variant_df

# Create the folders I will need
os.makedirs('./../../Data/Analysis_NovaSeq/read_abundances/', exist_ok = True)
os.makedirs('./../../Data/Analysis_NovaSeq/read_abundances/Codons/', exist_ok = True)
os.makedirs('./../../Data/Analysis_NovaSeq/read_abundances/Residues/', exist_ok = True)

## Loop through the files to get the read fractions
file_list = glob.glob('../../Data/Analysis_NovaSeq/amplicon_align/*')
for needle_file in file_list:
    ## Extract the sample name (remove the ".needle" extension)
    sample_name = os.path.basename(needle_file)[0:-7]
    
    read_fraction_df, variant_df = get_timepoint_fraction_df(needle_file, sample_name)
    
    print('----------------')

