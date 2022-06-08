This folder contains the scripts used to preprocess the raw sequencing data.

- 001_sequencing_data_processing.py: Does the quality control and trimming of the sequencing data before merging paired end reads and
aligning to the reference to aggregate all reads that map to the same mutant.

- 002_complete_datasets_MiSeq_NovaSeq.R: Takes the preprocessed data and the data from the other analyses to organize it in complete
datasets to be used in the scripts that generate the main and supplementary figures for the paper.
