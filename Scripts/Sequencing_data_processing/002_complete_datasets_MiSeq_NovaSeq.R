############################################################
####        Complete data MiSeq_NovaSeq                 ####
#### This script will generate a complete dataset with  ####
#### both the MiSeq and NovaSeq data and all the        ####
#### selection coefficients.                            ####
############################################################

# Load libraries
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(ggpubr)
library(xlsx)
library(BiocManager)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

theme_set(theme_cowplot())

# Set path to home folder (DfrB1_DMS_2022)
setwd('/path/to/DfrB1_DMS_2022/')

#### Load the NovaSeq data ####

metadata <- read.xlsx(
  'Data/SuppTables_ManuscriptDfrB1.xlsx',
                      sheetName = 'TableS2_DMS_sample_description', rowIndex = 1:89)

#### Load the NovaSeq data at the codon level ####

### The 001_sequencing_data_processing must have been run before to load the data
codon_file_list <- list.files(
  'Data/Analysis_NovaSeq/read_abundances/Codons/',
  include.dirs = F, full.names = T)

# Define the genetic code
codons <- c(
  'ATA', 'ATC', 'ATT', 'ATG',
  'ACA', 'ACC', 'ACG', 'ACT',
  'AAC', 'AAT', 'AAA', 'AAG',
  'AGC', 'AGT', 'AGA', 'AGG',
  'CTA', 'CTC', 'CTG', 'CTT',
  'CCA', 'CCC', 'CCG', 'CCT',
  'CAC', 'CAT', 'CAA', 'CAG',
  'CGA', 'CGC', 'CGG', 'CGT',
  'GTA', 'GTC', 'GTG', 'GTT',
  'GCA', 'GCC', 'GCG', 'GCT',
  'GAC', 'GAT', 'GAA', 'GAG',
  'GGA', 'GGC', 'GGG', 'GGT',
  'TCA', 'TCC', 'TCG', 'TCT',
  'TTC', 'TTT', 'TTA', 'TTG',
  'TAC', 'TAT', 'TAA', 'TAG',
  'TGC', 'TGT', 'TGA', 'TGG'
)

residues <- c(
  'I', 'I', 'I', 'M',
  'T', 'T', 'T', 'T',
  'N', 'N', 'K', 'K',
  'S', 'S', 'R', 'R',
  'L', 'L', 'L', 'L',
  'P', 'P', 'P', 'P',
  'H', 'H', 'Q', 'Q',
  'R', 'R', 'R', 'R',
  'V', 'V', 'V', 'V',
  'A', 'A', 'A', 'A',
  'D', 'D', 'E', 'E',
  'G', 'G', 'G', 'G',
  'S', 'S', 'S', 'S',
  'F', 'F', 'L', 'L',
  'Y', 'Y', '*', '*',
  'C', 'C', '*', 'W'
)

genetic_code <- data.frame(Codons = codons, Encoded_residues = residues)

all_codon_data <- c()

for(infile in codon_file_list){
  ## Extract the sample ID from the name
  sample_id <- str_split(string = basename(infile), pattern = '_')[[1]][4]
  
  # Read the file
  codon_df <- read_delim(delim = ',', col_names = T, file = infile)
  colnames(codon_df)[1] <- 'Codon'
  
  # Transform into a tidy formatted df and add the pool number
  new_codon_df <- codon_df %>% gather(-Codon, key = Position, value = read_abundance)
  new_codon_df %<>% mutate(Sample = sample_id)
  
  all_codon_data <- rbind(all_codon_data, new_codon_df)
  
}

all_codon_data$Position <- as.numeric(all_codon_data$Position)
all_codon_data$Sample <- as.numeric(all_codon_data$Sample)

# Add the arabinose, TMP, timepoint from the metadata
codon_read_abundances <- inner_join(x = all_codon_data %>% mutate(Sequencer = 'NovaSeq'), 
                                    y = metadata %>% 
                                      select(Sample, ID, Experiment, Timepoint, Arabinose, TMP, Biological.replicate,
                                             Technical.replicate, Sequencing.platform),
                                    by = c('Sample' = 'Sample', 'Sequencer' = 'Sequencing.platform'))

## Add the data for the WT codons
wt_matrix <- read_delim(file = 'Data/WT_sequence_table.txt', delim = '\t')

# Add the WT residues
codon_read_abundances <- left_join(x = codon_read_abundances, 
                                   y = wt_matrix, 
                                   by = c('Position' = 'Position'))

# Add the encoded residue
codon_read_abundances <- left_join(x = codon_read_abundances, 
                                   y = genetic_code, by = c('Codon' = 'Codons'))

## Remove UAG codon and group by encoded residue to aggregate read abundances
novaseq_read_abundances_res <- codon_read_abundances %>% ungroup() %>%
  filter(Codon != 'TAG') %>%
  group_by(Position, Sample, Sequencer, Timepoint, Arabinose, TMP, Biological.replicate, 
           Technical.replicate, Experiment, ID, WT_Residue, Encoded_residues) %>%
  summarise(read_abundance = sum(read_abundance))

# Get the median of WT for each sample
wt_medians <- codon_read_abundances %>% rowwise() %>%
  filter(Codon == WT_Codon) %>% 
  group_by(Sample, ID, Timepoint, Arabinose, TMP, Biological.replicate) %>%
  summarise(median_wt = median(read_abundance))

# Add these medians to the dataframe of read abundances
codon_read_abundances <- left_join(x = codon_read_abundances, 
                                   y = wt_medians %>% ungroup() %>% select(Sample, median_wt),
                                   by = c('Sample' = 'Sample'))

# Calculate Rmut / Rwt for all the data
codon_read_abundances %<>% mutate(mut_wt_ratio = read_abundance / median_wt)

# Separate the t0 samples and match them to their respective t10
codon_read_abundances_t0 <- codon_read_abundances %>% filter(Timepoint == 0)
codon_read_abundances_t10 <- codon_read_abundances %>% filter(Timepoint != 0)

# Average technical replicates
codon_read_abundances_t0_techavg <- codon_read_abundances_t0 %>% ungroup() %>%
  group_by(WT_Codon, Experiment, Codon, Position, Timepoint, Arabinose, TMP, Sequencer, 
           Biological.replicate, WT_Residue, Encoded_residues) %>%
  summarise(mean_mut_wt_ratio_t0 = mean(mut_wt_ratio))

codon_read_abundances_t10_techavg <- codon_read_abundances_t10 %>% ungroup() %>%
  group_by(WT_Codon, Experiment, Codon, Position, Timepoint, Arabinose, TMP, Sequencer, 
           Biological.replicate, WT_Residue, Encoded_residues) %>%
  summarise(mean_mut_wt_ratio = mean(mut_wt_ratio))

codon_read_abundances_matched <- left_join(x = codon_read_abundances_t10_techavg %>% ungroup(), 
                                           y = codon_read_abundances_t0_techavg %>% ungroup() %>%
                                             select(WT_Codon, Codon, Experiment, Position, mean_mut_wt_ratio_t0, TMP, 
                                                    Arabinose), 
                                           by = c('Codon' = 'Codon', 'Position' = 'Position', 'TMP' = 'TMP',
                                                  'Arabinose' = 'Arabinose', 'WT_Codon' = 'WT_Codon', 
                                                  'Experiment' = 'Experiment'))

# Calculate the selection coefficient
codon_read_abundances_matched %<>% mutate(Timepoint = as.numeric(Timepoint)) %>%
  mutate(sel_coeff = log2(mean_mut_wt_ratio / mean_mut_wt_ratio_t0) / Timepoint)

# Add the IDs now that the technical replicates have been averaged and selection
# coefficients have been calculated
codon_read_abundances_matched %<>% rowwise() %>%
  mutate(temp_TMP = ifelse(TMP == 10, 'T', 'NT'), 
         temp_sequencer = ifelse(Sequencer == 'NovaSeq', 'N', 'M'), 
         temp_Arabinose = toString(format(round(as.numeric(Arabinose), 3), nsmall = 3))
  ) %>%
  mutate(ID = str_c('E', toString(Experiment), 
                    '.BR', toString(Biological.replicate),
                    '.', substr(temp_Arabinose, start = 3,
                                stop = nchar(temp_Arabinose)
                    ),
                    '.', temp_TMP,
                    '.S', temp_sequencer,
                    sep = ''
  )
  )


## Remove unneeded columns
all_novaseq_data <- codon_read_abundances_matched %>%
  select(Position, WT_Codon, Codon, WT_Residue, Encoded_residues, ID, Timepoint, Arabinose, TMP, 
         sel_coeff, Sequencer)

#### Load the MiSeq data with TMP ####

### The 001_sequencing_data_processing must have been run before to load the data
codon_file_list <- list.files('Data/Analysis_MiSeq/read_abundances/Codons/', include.dirs = F,
                              full.names = T)

all_codon_data <- c()

for(infile in codon_file_list){
  # Extract the sample ID from the name
  sample_id <- str_split(string = basename(infile), pattern = '_')[[1]][1]
  
  # Read the file
  codon_df <- read_delim(delim = ',', col_names = T, file = infile)
  colnames(codon_df)[1] <- 'Codon'
  
  # Transform into a tidy formatted df and add the pool number
  new_codon_df <- codon_df %>% gather(-Codon, key = Position, value = read_abundance)
  new_codon_df %<>% mutate(Sample = sample_id)
  
  all_codon_data <- rbind(all_codon_data, new_codon_df)
  
}

all_codon_data$Position <- as.numeric(all_codon_data$Position)
all_codon_data$Sample <- as.numeric(all_codon_data$Sample)

# Add the arabinose, TMP, timepoint from the metadata
codon_read_abundances_miseq <- inner_join(x = all_codon_data %>% mutate(Sequencer = 'MiSeq'), 
                                          y = metadata %>% 
                                            select(Sample, Timepoint, Arabinose, TMP, Biological.replicate, 
                                                   Technical.replicate, Sequencing.platform, Experiment, ID),
                                          by = c('Sample' = 'Sample', 'Sequencer' = 'Sequencing.platform'))

## Calculate selection coefficients

# Add the WT residues
codon_read_abundances_miseq <- left_join(x = codon_read_abundances_miseq, 
                                         y = wt_matrix, 
                                         by = c('Position' = 'Position'))

# Add the encoded residue
codon_read_abundances_miseq <- left_join(x = codon_read_abundances_miseq, 
                                         y = genetic_code, by = c('Codon' = 'Codons'))

## Remove UAG codon and group by encoded residue to aggregate read abundances
miseq_read_abundances_res <- codon_read_abundances_miseq %>% ungroup() %>%
  filter(Codon != 'TAG') %>%
  group_by(Position, Sample, Sequencer, Timepoint, Arabinose, TMP, Biological.replicate, 
           Technical.replicate, Experiment, ID, WT_Residue, Encoded_residues) %>%
  summarise(read_abundance = sum(read_abundance))

# Get the median of WT for each sample
wt_medians_miseq <- codon_read_abundances_miseq %>% rowwise() %>%
  filter(Codon == WT_Codon) %>% 
  group_by(Sample, ID, Timepoint, Arabinose, TMP, Biological.replicate) %>%
  summarise(median_wt = median(read_abundance))

# Add these medians to the dataframe of read abundances
codon_read_abundances_miseq <- left_join(x = codon_read_abundances_miseq, 
                                         y = wt_medians_miseq,
                                         by = c('Sample' = 'Sample', 'Arabinose' = 'Arabinose', 
                                                'Timepoint' = 'Timepoint', 'TMP' = 'TMP', 
                                                'Biological.replicate' = 'Biological.replicate', 
                                                'ID' = 'ID'))

# Calculate Rmut / Rwt for all the data
codon_read_abundances_miseq %<>% mutate(mut_wt_ratio = read_abundance / median_wt)

# Separate the t0 samples and match them to their respective t10
codon_read_abundances_miseq_t0 <- codon_read_abundances_miseq %>% filter(Timepoint == 0)
codon_read_abundances_miseq_t10 <- codon_read_abundances_miseq %>% filter(Timepoint != 0)

codon_read_abundances_matched_miseq <- left_join(x = codon_read_abundances_miseq_t10, 
                                                 y = codon_read_abundances_miseq_t0 %>% 
                                                   mutate(mut_wt_ratio_t0 = mut_wt_ratio) %>%
                                                   select(Codon, Position, Experiment, mut_wt_ratio_t0, TMP, Arabinose), 
                                                 by = c('Codon' = 'Codon', 'Position' = 'Position',
                                                        'TMP' = 'TMP', 'Arabinose' = 'Arabinose', 
                                                        'Experiment' = 'Experiment'))

# Calculate the selection coefficient
codon_read_abundances_matched_miseq %<>% mutate(Timepoint = as.numeric(Timepoint)) %>%
  mutate(sel_coeff = log2(mut_wt_ratio / mut_wt_ratio_t0) / Timepoint)

# Add the IDs now that the technical replicates have been averaged and selection
# coefficients have been calculated
codon_read_abundances_matched_miseq %<>% rowwise() %>%
  mutate(temp_TMP = ifelse(TMP == 10, 'T', 'NT'), 
         temp_sequencer = ifelse(Sequencer == 'NovaSeq', 'N', 'M'), 
         temp_Arabinose = toString(format(round(as.numeric(Arabinose), 3), nsmall = 3))
  ) %>%
  mutate(ID = str_c('E', toString(Experiment), 
                    '.BR', toString(Biological.replicate),
                    '.', substr(temp_Arabinose, start = 3,
                                stop = nchar(temp_Arabinose)
                    ),
                    '.', temp_TMP,
                    '.S', temp_sequencer,
                    sep = ''
  )
  )

all_miseq_data <- bind_rows(codon_read_abundances_matched_miseq %>% mutate(Sequencer = 'MiSeq') %>%
                              select(Position, WT_Codon, Codon, WT_Residue, Encoded_residues,
                                     ID, Timepoint, Arabinose, TMP,
                                     sel_coeff, Sequencer)) 

#### Put the MiSeq and NovaSeq data together (codons) ####

colnames(all_novaseq_data)
colnames(all_miseq_data)

# Average selection coefficients of all replicates
all_data_all_reps <- bind_rows(all_novaseq_data %>% 
                                 mutate(Timepoint = as.numeric(Timepoint)),
                               all_miseq_data %>%
                                 mutate(Timepoint = as.numeric(Timepoint))
) %>% ungroup() %>%
  group_by(Position, WT_Codon, Codon, WT_Residue, Encoded_residues,
           Timepoint, Arabinose, TMP) %>%
  summarise(mean_sel_coeff = mean(sel_coeff),
            sd_sel_coeff = sd(sel_coeff),
            sem_sel_coeff = sd(sel_coeff) / sqrt(n()), 
            num_samples = n())


## Concatenate all the data points
all_data_all_reps_notavg <- bind_rows(all_novaseq_data %>%
                                        mutate(Timepoint = as.numeric(Timepoint)),
                                      all_miseq_data %>%
                                        mutate(Timepoint = as.numeric(Timepoint))
)

# Arrange the tables by position
all_data_all_reps %<>% mutate(Position = as.numeric(Position)) %>%
  arrange(Position)
all_data_all_reps_notavg %<>% mutate(Position = as.numeric(Position)) %>%
  arrange(Position)


#### Save the datasets ####
write.table(x = all_data_all_reps, quote = F, sep = '\t', row.names = F, col.names = T, append = F, 
            file = 'Data/Complete_datasets/complete_dataset_avgBothSequencers_Codons.txt')
write.table(x = all_data_all_reps_notavg, quote = F, sep = '\t', row.names = F, col.names = T, append = F, 
            file = 'Data/Complete_datasets/all_data_all_reps_bothSequencers_Codons.txt')

#### Continue working with the data at the residue level ####

#### NovaSeq residue level ####

# Get the median of WT for each sample
wt_medians_novaseq <- novaseq_read_abundances_res  %>% ungroup() %>% rowwise() %>%
  filter(Encoded_residues == WT_Residue) %>% 
  group_by(Sample, ID, Timepoint, Arabinose, TMP, Biological.replicate) %>%
  summarise(median_wt = median(read_abundance))

# Add these medians to the dataframe of read abundances
residue_read_abundances_novaseq <- left_join(x = novaseq_read_abundances_res, 
                                           y = wt_medians,
                                           by = c('Sample' = 'Sample', 'Arabinose' = 'Arabinose', 
                                                  'Timepoint' = 'Timepoint', 'TMP' = 'TMP', 
                                                  'Biological.replicate' = 'Biological.replicate', 
                                                  'ID' = 'ID'))

# Calculate Rmut / Rwt for all the data
residue_read_abundances_novaseq %<>% mutate(mut_wt_ratio = read_abundance / median_wt)

# Separate the t0 samples and match them to their respective t10
residue_read_abundances_novaseq_t0 <- residue_read_abundances_novaseq %>% filter(Timepoint == 0)
residue_read_abundances_novaseq_t10 <- residue_read_abundances_novaseq %>% filter(Timepoint != 0)

# Average technical replicates
residue_read_abundances_novaseq_t0_techavg <- residue_read_abundances_novaseq_t0 %>% ungroup() %>%
  group_by(WT_Residue, Experiment, Encoded_residues, Position, Timepoint, Arabinose, TMP, Sequencer, 
           Biological.replicate) %>%
  summarise(mean_mut_wt_ratio_t0 = mean(mut_wt_ratio))

residue_read_abundances_novaseq_t10_techavg <- residue_read_abundances_novaseq_t10 %>% ungroup() %>%
  group_by(WT_Residue, Experiment, Encoded_residues, Position, Timepoint, Arabinose, TMP, Sequencer, 
           Biological.replicate) %>%
  summarise(mean_mut_wt_ratio = mean(mut_wt_ratio))


residue_read_abundances_novaseq_matched <-
  left_join(x = residue_read_abundances_novaseq_t10_techavg %>% ungroup(),
            y = residue_read_abundances_novaseq_t0_techavg %>% ungroup() %>%
              select(WT_Residue, Encoded_residues, Experiment, Position,
                     mean_mut_wt_ratio_t0, TMP,Arabinose), 
            by = c('Encoded_residues' = 'Encoded_residues', 'Position' = 'Position', 'TMP' = 'TMP',
                   'Arabinose' = 'Arabinose', 'WT_Residue' = 'WT_Residue', 
                   'Experiment' = 'Experiment'))

# Calculate the selection coefficient
residue_read_abundances_novaseq_matched %<>% mutate(Timepoint = as.numeric(Timepoint)) %>%
  mutate(sel_coeff = log2(mean_mut_wt_ratio / mean_mut_wt_ratio_t0) / Timepoint)

# Add the IDs now that the technical replicates have been averaged and selection
# coefficients have been calculated
residue_read_abundances_novaseq_matched %<>% rowwise() %>%
  mutate(temp_TMP = ifelse(TMP == 10, 'T', 'NT'), 
         temp_sequencer = ifelse(Sequencer == 'NovaSeq', 'N', 'M'), 
         temp_Arabinose = toString(format(round(as.numeric(Arabinose), 3), nsmall = 3))
  ) %>%
  mutate(ID = str_c('E', toString(Experiment), 
                    '.BR', toString(Biological.replicate),
                    '.', substr(temp_Arabinose, start = 3,
                                stop = nchar(temp_Arabinose)
                    ),
                    '.', temp_TMP,
                    '.S', temp_sequencer,
                    sep = ''
  )
  )


## Remove unneeded columns
all_novaseq_data <- residue_read_abundances_novaseq_matched %>%
  select(Position, WT_Residue, Encoded_residues, ID, Timepoint, Arabinose, TMP, 
         sel_coeff, Sequencer)


#### MiSeq residue level ####
# Get the median of WT for each sample
wt_medians_miseq <- miseq_read_abundances_res  %>% ungroup() %>% rowwise() %>%
  filter(Encoded_residues == WT_Residue) %>% 
  group_by(Sample, ID, Timepoint, Arabinose, TMP, Biological.replicate) %>%
  summarise(median_wt = median(read_abundance))

# Add these medians to the dataframe of read abundances
residue_read_abundances_miseq <- left_join(x = miseq_read_abundances_res, 
                                           y = wt_medians_miseq,
                                           by = c('Sample' = 'Sample', 'Arabinose' = 'Arabinose', 
                                                  'Timepoint' = 'Timepoint', 'TMP' = 'TMP', 
                                                  'Biological.replicate' = 'Biological.replicate', 
                                                  'ID' = 'ID'))

# Calculate Rmut / Rwt for all the data
residue_read_abundances_miseq %<>% mutate(mut_wt_ratio = read_abundance / median_wt)

# Separate the t0 samples and match them to their respective t10
residue_read_abundances_miseq_t0 <- residue_read_abundances_miseq %>% filter(Timepoint == 0)
residue_read_abundances_miseq_t10 <- residue_read_abundances_miseq %>% filter(Timepoint != 0)

residue_read_abundances_matched_miseq <- left_join(x = residue_read_abundances_miseq_t10 %>% ungroup(), 
                                                   y = residue_read_abundances_miseq_t0 %>% ungroup() %>%
                                                     mutate(mut_wt_ratio_t0 = mut_wt_ratio) %>%
                                                     select(Encoded_residues, Position, Experiment, mut_wt_ratio_t0, TMP, Arabinose), 
                                                   by = c('Encoded_residues' = 'Encoded_residues',
                                                          'Position' = 'Position',
                                                          'TMP' = 'TMP', 'Arabinose' = 'Arabinose', 
                                                          'Experiment' = 'Experiment'))

# Calculate the selection coefficient
residue_read_abundances_matched_miseq %<>% mutate(Timepoint = as.numeric(Timepoint)) %>%
  mutate(sel_coeff = log2(mut_wt_ratio / mut_wt_ratio_t0) / Timepoint)

# Add the IDs now that the technical replicates have been averaged and selection
# coefficients have been calculated
residue_read_abundances_matched_miseq %<>% rowwise() %>%
  mutate(temp_TMP = ifelse(TMP == 10, 'T', 'NT'), 
         temp_sequencer = ifelse(Sequencer == 'NovaSeq', 'N', 'M'), 
         temp_Arabinose = toString(format(round(as.numeric(Arabinose), 3), nsmall = 3))
  ) %>%
  mutate(ID = str_c('E', toString(Experiment), 
                    '.BR', toString(Biological.replicate),
                    '.', substr(temp_Arabinose, start = 3,
                                stop = nchar(temp_Arabinose)
                    ),
                    '.', temp_TMP,
                    '.S', temp_sequencer,
                    sep = ''
  )
  )


all_miseq_data <- residue_read_abundances_matched_miseq %>%
  select(Position, WT_Residue, Encoded_residues, ID, Timepoint, Arabinose, TMP, 
         sel_coeff, Sequencer)

colnames(all_novaseq_data)
colnames(all_miseq_data)

# Average selection coefficients of all replicates
all_data_all_reps <- bind_rows(all_novaseq_data %>% 
                                 mutate(Timepoint = as.numeric(Timepoint)),
                               all_miseq_data %>%
                                 mutate(Timepoint = as.numeric(Timepoint))
) %>% ungroup() %>%
  mutate(Residue = Encoded_residues) %>%
  group_by(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP) %>%
  summarise(mean_sel_coeff = mean(sel_coeff),
            sd_sel_coeff = sd(sel_coeff),
            sem_sel_coeff = sd(sel_coeff) / sqrt(n()), 
            num_samples = n())


## Concatenate all the data points
all_data_all_reps_notavg <- bind_rows(all_novaseq_data %>%
                                        mutate(Timepoint = as.numeric(Timepoint)),
                                      all_miseq_data %>%
                                        mutate(Timepoint = as.numeric(Timepoint))
) %>% mutate(Residue = Encoded_residues) %>% select(-Encoded_residues)

#### Add the rest of the columns from the other data sources ####
# Prepare a table of amino acid names
aa_three2one <- data.frame(cbind(c('A', 'R', 'D', 'N', 'C', 
                                   'E', 'Q', 'G', 'H', 'I', 
                                   'L', 'K', 'M', 'F', 'P',
                                   'S', 'T', 'W', 'Y', 'V'),
                                 c('ALA', 'ARG', 'ASN', 'ASP', 'CYS',
                                   'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
                                   'LEU', 'LYS', 'MET', 'PHE', 'PRO',
                                   'SER', 'THR', 'TRP', 'TYR', 'VAL')))
colnames(aa_three2one) <- c('One-letter', 'Three-letter')

## Stability ddGs
HET_stab<-read_delim('Data/In_silico_mutagenesis_results/006_gathered_data/2rk1_monomer_HM/stability_monomer_ddGs_tabs.txt', 
                     delim='\t', col_names = c('Chain', 'WT_res', 'Mut_res', 'Mean_ddG_stab_HET', 'Std_dev_ddG_stab_HET', 'Min_ddG_stab_HET', 'Max_ddG_stab_HET'))

## Binding energy ddGs (only HMs because this experiment only works with one gene copy)
HM_int_A_B<-read_delim('Data/In_silico_mutagenesis_results/006_gathered_data/2rk1_HM/A-B/interface_ddGs_tabs.txt', 
                       delim='\t', col_names = c('Chain', 'WT_res', 'Mut_res', 'Mean_ddG_int_HM', 'Std_dev_ddG_int_HM', 'Min_ddG_int_HM', 'Max_ddG_int_HM'))

HM_int_A_C<-read_delim('Data/In_silico_mutagenesis_results/006_gathered_data/2rk1_HM/A-C/interface_ddGs_tabs.txt', 
                       delim='\t', col_names = c('Chain', 'WT_res', 'Mut_res', 'Mean_ddG_int_HM', 'Std_dev_ddG_int_HM', 'Min_ddG_int_HM', 'Max_ddG_int_HM'))

HM_int_A_D<-read_delim('Data/In_silico_mutagenesis_results/006_gathered_data/2rk1_HM/A-D/interface_ddGs_tabs.txt', 
                       delim='\t', col_names = c('Chain', 'WT_res', 'Mut_res', 'Mean_ddG_int_HM', 'Std_dev_ddG_int_HM', 'Min_ddG_int_HM', 'Max_ddG_int_HM'))

# For the dataframes on mutational effects, separate the WT residue from the position
HET_stab %<>% separate(col = WT_res, into = c('WT', 'Position'), sep = 1) %>%
  mutate(Position = as.numeric(Position))
HM_int_A_B %<>% separate(col = WT_res, into = c('WT', 'Position'), sep = 1) %>%
  mutate(Position = as.numeric(Position))
HM_int_A_C %<>% separate(col = WT_res, into = c('WT', 'Position'), sep = 1) %>%
  mutate(Position = as.numeric(Position))
HM_int_A_D %<>% separate(col = WT_res, into = c('WT', 'Position'), sep = 1) %>%
  mutate(Position = as.numeric(Position))

# Start joining the data
all_data_tmp1 <- left_join(x = all_data_all_reps, 
                           y = HET_stab %>% select(Position, Mut_res, Mean_ddG_stab_HET),
                           by = c('Residue' = 'Mut_res', 'Position' = 'Position'))

all_data_tmp2 <- left_join(x = all_data_tmp1,
                           y = HM_int_A_B %>% mutate(Mean_ddG_int_HM_A_B = Mean_ddG_int_HM) %>%
                             select(Position, Mut_res, Mean_ddG_int_HM_A_B),
                           by = c('Residue' = 'Mut_res', 'Position' = 'Position'))

all_data_tmp3 <- left_join(x = all_data_tmp2,
                           y = HM_int_A_C %>% mutate(Mean_ddG_int_HM_A_C = Mean_ddG_int_HM) %>%
                             select(Position, Mut_res, Mean_ddG_int_HM_A_C),
                           by = c('Residue' = 'Mut_res', 'Position' = 'Position'))

all_data_mutEffects <- left_join(x = all_data_tmp3,
                                 y = HM_int_A_D %>% mutate(Mean_ddG_int_HM_A_D = Mean_ddG_int_HM) %>%
                                   select(Position, Mut_res, Mean_ddG_int_HM_A_D),
                                 by = c('Residue' = 'Mut_res', 'Position' = 'Position'))

#### Add data on Shannon entropy ####

# Read the Shannon entropy data
entropy <- read_delim(file = 'Data/Evolutionary_analysis/DHFR_entropy.txt', 
                      delim = '\t', col_names = c('Entropy')) %>%
  mutate(Position = row_number())

# Add the data on Shannon entropy
all_data_mutEffects_Entropy <- left_join(x = all_data_mutEffects, y = entropy,
                                         by = c('Position' = 'Position'))

#### Read the secondary structure annotation data ####

sec_struc <- read_delim('Data/Structural_data/DHFR_2rk1_DSSP_table_bio_rSASA.txt', delim = '\t')

all_data_mutEffects_Entropy_SASA <- left_join(x = all_data_mutEffects_Entropy, y = sec_struc,
                                              by = c('Position' = 'Position'))

# Load matrix of aminoacid indices
aa_indices <- read_delim(file = 'Data/aminoacid_inidices/002_all_data/all_indices_final_table_propensity.txt', delim = '\t')

## Left joins to add amino acid indices for WT and mutant codons (two separate tables)
# Left joins allow us to keep mutations to stop codons
all_codon_means_wt_indices <- left_join(x = all_data_mutEffects_Entropy_SASA,
                                        y = aa_indices, 
                                        by = c('WT_Residue' = 'Aminoacid.1.letter'))

all_codon_means_mut_indices <- left_join(x = all_data_mutEffects_Entropy_SASA,
                                         y = aa_indices, 
                                         by = c('Residue' = 'Aminoacid.1.letter'))


# Get the subtraction to get the final matrix
change_matrix <- all_codon_means_mut_indices[,19:ncol(all_codon_means_mut_indices)] - all_codon_means_wt_indices[,19:ncol(all_codon_means_wt_indices)]
# all_codon_means_change_indices <- cbind(all_codon_means_wt_indices[,1:5], change_matrix)
all_data_complete <- all_codon_means_wt_indices
all_data_complete[, 19:ncol(all_data_complete)] <- change_matrix

# Repeat the join for the data that were not averaged
all_data_complete_notavg <- left_join(x = all_data_all_reps_notavg %>% ungroup(), 
                                      y = all_data_complete %>% ungroup() %>%
                                        filter(Arabinose == 0.01, TMP == 0, Timepoint == 5) %>% # To avoid repeating values
                                        select(-Timepoint, -Arabinose, -TMP, -mean_sel_coeff, -sd_sel_coeff,
                                               -sem_sel_coeff, -num_samples), 
                                      by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                                             'Residue' = 'Residue')
)

#### Save the datasets ####
write.table(x = all_data_complete, quote = F, sep = '\t', row.names = F, col.names = T, append = F, 
            file = 'Data/Complete_datasets/complete_dataset_avgBothSequencers.txt')
write.table(x = all_data_complete_notavg, quote = F, sep = '\t', row.names = F, col.names = T, append = F, 
            file = 'Data/Complete_datasets/all_data_all_reps_bothSequencers.txt')
