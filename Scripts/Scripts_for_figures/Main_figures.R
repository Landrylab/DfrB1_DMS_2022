###########################################################
####        DMS_figures_paper_2022-05-03               ####
#### This script will organize the data and prepare    ####
#### the figures for the paper. This new script will   ####
#### get the updated figures with Moira's demux data.  ####
#### This script includes the ara 0.4 data             ####
###########################################################

# Load libraries
library(tidyverse)
library(magrittr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(ggpubr)

library(BiocManager)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(wesanderson)
library(xlsx)
library(growthcurver)

library(car)
library(ggrepel)
library(nlme)
library(ARTool)
library(lme4)
library(rlme)
library(FDRestimation)
library(gridGraphics)
library(shadowtext)
library(viridis)
# library(nls)

theme_set(theme_cowplot() + 
            theme(panel.background = element_rect(fill = 'white'),
                  plot.background = element_rect(fill = 'white'))
          ) 

## Set working directory
setwd('/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/Angel/PhD_projects/R67_DMS_December2020/')

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

# Load the complete dataset, calculate selection coefficients from normalized scores
all_data_complete <- read_delim(
  '../R67_DMS_February2022/Complete_datasets_MiSeq_NovaSeq/complete_dataset_avgBothSequencers.txt', delim = '\t')
  
#### Figure 1: Study system ####

#### Figure 1C: Effect of arabinose on fitness ####

annotation_df <- read_delim(
  'Data/Data_cytometry/May_2022/Stats_run2_05_05_22_matched.csv',
  locale = locale(decimal_mark = ','), delim = ',')

# Use a loop to load the rest of the data
# path_files <- 'Data/Data_cytometry/All_raw_files'
path_files <- 'Data/Data_cytometry/May_2022/All_raw_files/'
list_files <- list.files(path_files)

all_data_expression <- c()

for(infile in list_files){
  # Load the data
  new_data <- read_delim(file.path(path_files, infile), delim = ',', col_names = T) %>%
    select(TIME, 'GRN-B-HLin', 'FSC-HLin', 'SSC-HLin', 'FSC-HLog', 'SSC-HLog')
  
  
  # Add the name of the file as an ID (remove the csv)
  new_data %<>% mutate(ID = substr(x = infile, start = 1, stop = (nchar(infile) - 4)))
  
  all_data_expression <- bind_rows(all_data_expression, new_data)
}

# Add the annotation data to identify each sample
all_data_expression %<>% separate(col = ID, into = c('tmp', 'Well'), sep = 'am.')

all_data_expression <- inner_join(x = all_data_expression, 
                       y = annotation_df,
                       by = c('Well' = 'Well'))


# Add the logarithms of the green fluorescence
all_data_expression %<>% mutate(GRNBHLog = ifelse(log10(`GRN-B-HLin`) < 0, 0, log10(`GRN-B-HLin`)),
                     `FSC-HLog` = ifelse(log10(`FSC-HLin`) < 0, 0, log10(`FSC-HLin`)),
                     `SSC-HLog` = ifelse(log10(`SSC-HLin`) < 0, 0, log10(`SSC-HLin`))
)

all_data_processed_expression <- all_data_expression %>% rowwise() %>%
  mutate(circle_test = (`FSC-HLog` - 1.25)^2 + (`SSC-HLog` - 1.25)^2) %>%
  filter(circle_test < 1) %>%
  # mutate(Time = factor(Time, levels = c('Stationary phase', 'After 5 generations', 'After 10 generations'))) %>%
  mutate(
    # Sample = str_c('Sample', toString(Sample), sep = ' '),
    Arabinose = str_c(toString(Arabinose), '% arabinose', sep = '') #,
    # Replicate = ifelse(Sample <= 6, 1, 
    #                    ifelse(Sample <= 24, floor((Sample - 1) / 6),
    #                           ifelse(Sample > 24, floor((Sample -1) / 6) - 3, NA)
    #                    )
  )

# all_data_processed %<>% arrange(Sample)

all_data_processed_expression %<>% 
  mutate(Arabinose = factor(Arabinose, levels = c('0% arabinose', 
                                                  '0.01% arabinose', 
                                                  '0.025% arabinose', 
                                                  '0.05% arabinose',
                                                  '0.2% arabinose',
                                                  '0.4% arabinose')
  )
  ) 

all_data_processed_summary_expression <- all_data_processed_expression %>%
  ungroup() %>%
  # group_by(Arabinose, Time, Replicate) %>%
  group_by(Arabinose, Well) %>%
  summarise(
    # med_fluo = median(GRNBHLog)
    med_fluo = mean(GRNBHLog),
    sem_fluo = sd(GRNBHLog) / sqrt(n()),
    num_cells = n()
    )

all_data_processed_summary_expression %<>%
  mutate(exp_level = ifelse(Arabinose == '0.2% arabinose', 'Optimal',
                            ifelse(Arabinose == '0.05% arabinose', 'Near-optimal',
                                   ifelse(Arabinose == '0.025% arabinose', 'Suboptimal',
                                          ifelse(Arabinose == '0.01% arabinose', 'Weak',
                                                 ifelse(Arabinose == '0.4% arabinose', 'Overexpressed',
                                                        'Absent')))))
  )



#### Figures for growth curve validations ####
plate.ind <- 'Data/GrowthCurves_February2022/Fitness_TMP_ID_IGA_23_02_22.xlsx'
file.od <- 'Data/GrowthCurves_February2022/Data_23_02_2022_bact.xlsx'

## Define function to read plate data
# function to process plates ----------------------------------------------
read.my.gc <- function(file, plate.index){
  pl <- read.xlsx(file,sheetIndex = 1, rowIndex = 5:100, stringsAsfators = FALSE,
                  header = F)
  ind <- read.xlsx(plate.index, sheetIndex = 1, rowIndex = 1:97, header = T) 
  
  time <- seq(0,0.30*(ncol(pl)-2), 0.30)
  
  colnames(pl)[1] <- "Well"
  colnames(pl)[2:ncol(pl)] <- time
  pl %<>% select(1:(ncol(pl)-2)) 
  
  data.pl <- gather(pl, key = "time", value = "OD",2:ncol(pl), convert = F)
  data.pl$time <- as.numeric(data.pl$time)
  data.pl$OD <- as.numeric(data.pl$OD)
  data.pl %<>% left_join(ind, by = "Well")
  
  d <-  select(data.pl, 1:3)
  d %<>% spread(key = "Well", value = "OD", convert = F) 
  colnames(d)[1] <- "time"
  
  ## Use Growthcurver
  ## For the WT controls I will trim at t = 13.5
  gc_out <- SummarizeGrowthByPlate(d, t_trim =13.5)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
}

data.od1 <- read.my.gc(file.od, plate.ind)

# Subtract the blank for this experiment and multiply by 5 to make it OD / mL
# (experiment was carried out in 0.2 mL)
data.od1 %<>% mutate(OD = (OD - 0.088) * 5)

data.od1.summary <- data.od1 %>%
  group_by(TMP, Arabinose, Replicate, Mutant) %>%
  summarise(auc = mean(auc_e))

data.od1.summary.wide <- data.od1.summary %>% ungroup() %>%
  filter(Mutant == 'WT') %>%
  pivot_wider(names_from = TMP, values_from = auc) %>%
  mutate(auc_diff = `10` - `0`)

# Merge all the data
fig1C_data <- left_join(x = data.od1.summary.wide, #  %>% filter(Arabinose != 0.4), 
                        y = all_data_processed_summary_expression %>% # filter(Time == 'After 10 generations') %>%
                          separate(Arabinose, into = c('Arabinose', 'tmp'), sep = '%') %>%
                          mutate(Arabinose = as.numeric(Arabinose)) %>% ungroup() %>%
                          group_by(Arabinose) %>% 
                          summarise(mean_fluo = round(mean(med_fluo), 2)),
                        by = c('Arabinose' = 'Arabinose'))

# Add a column for the maximum difference (0% arabinose)
max_diff <- fig1C_data %>% filter(Arabinose == 0) %>% group_by(Arabinose) %>%
  summarise(max_diff = median(auc_diff))

# Calculate the percentage of recovered growth as 100*(1 - (curr_diff / max_diff))
fig1C_data %<>% mutate(max_auc_diff = max_diff$max_diff[1]) %>%
  mutate(pct_recovery = 100*(1 - (auc_diff / max_auc_diff)))

fig1C_data %<>%
  mutate(exp_level = ifelse(Arabinose == 0.2, 'Optimal', 
                            ifelse(Arabinose == 0.05, 'Near-optimal',
                                   ifelse(Arabinose == 0.025, 'Suboptimal', 
                                          ifelse(Arabinose == 0.01, 'Weak', 
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        'Absent')))))
  )

p_fig1C <- fig1C_data %>%
  mutate(Arabinose = as.factor(Arabinose), 
         exp_level = factor(exp_level, 
                            levels = c('Absent', 'Weak', 'Suboptimal',
                                       'Near-optimal', 'Optimal', 'Overexpressed'))) %>%
  ggplot(aes(x = Arabinose, y = pct_recovery, colour = exp_level)) +
  geom_jitter(width = 0.4, size = 3) +
  xlab('Arabinose concentration (% m/v)') +
  ylab('Recovered growth (%)') +
  labs(colour = 'Expression level') +
  scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = 'top', 
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank())
p_fig1C + theme(legend.position = 'none')

## Let's see Fig 1C with boxplots
p_fig1C <- fig1C_data %>%
  mutate(Arabinose = as.factor(Arabinose), 
         exp_level = factor(exp_level, 
                            levels = c('Absent', 'Weak', 'Suboptimal',
                                       'Near-optimal', 'Optimal', 'Overexpressed'))) %>%
  ggplot(aes(x = Arabinose, y = pct_recovery, colour = exp_level, fill = exp_level)) +
  geom_boxplot(width = 0.4, alpha = 0.4) +
  geom_jitter(size = 3) +
  xlab('Arabinose concentration (% m/v)') +
  ylab('Recovered growth (%)') +
  labs(colour = 'Expression level') +
  scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  scale_fill_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = 'top', 
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank())
p_fig1C + theme(legend.position = 'none')

## Let's try with a big dot for the mean with errorbars
fig1C_data_sum <- fig1C_data %>% ungroup() %>% 
  group_by(Arabinose, exp_level) %>%
  summarise(mean_recovery = mean(pct_recovery), 
            sem_recovery = sd(pct_recovery) / sqrt(n()), 
            num_samples = n())

# Add an annotation box for the legend of figure 1C
annotation_box <- roundrectGrob(x = unit(0.7, 'npc'), y = unit(0.3, 'npc'), 
                                width = unit(0.3, 'npc'), height = unit(0.4, 'npc'), 
                                gp = gpar(lwd = 2, col = 'black',
                                          # fill = '#999999', 
                                          fill = '#737373',
                                          lty = 2))

p_fig1C <- fig1C_data_sum %>% ungroup() %>%
  mutate(Arabinose = as.factor(Arabinose), 
         exp_level = factor(exp_level, 
                            levels = c('Absent', 'Weak', 'Suboptimal',
                                       'Near-optimal', 'Optimal', 'Overexpressed'))
         ) %>%
  ggplot(aes(x = Arabinose, y = mean_recovery,
             ymax = mean_recovery + 2.5*sem_recovery, 
             ymin = mean_recovery - 2.5*sem_recovery, 
             colour = exp_level)
         ) +
  # geom_line(aes(group = -1), size = 1.5, show.legend = F) +
  geom_point(size = 3) +
  # geom_jitter(data = fig1C_data, alpha = 0.4, inherit.aes = F, width = 0, size = 3,
  #             aes(x = as.factor(Arabinose), y = pct_recovery, colour = exp_level)) +
  geom_errorbar(width = 0.2, size = 0.7) +
  xlab('Arabinose concentration (% m/v)') +
  ylab('Recovered growth (%)') +
  labs(colour = 'Expression level') +
  scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20, face = 'bold'),
        legend.text = element_text(size = 18),
        # legend.position = 'top', 
        legend.position = c(0.7, 0.3),
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  annotation_custom(annotation_box)
  
p_fig1C

#### Figure 1D: Showing the medians of the distribution of GFP fluorescence ####

## Load the main table with the sample IDs
# annotation_df <- read.xlsx(
#   'Data/Data_cytometry/Cytometry_file_names_2021_06_02_matched.xlsx',
#   sheetIndex = 1, rowIndex = 1:43, header = T)

# 
# p_fig1C <- all_data_processed_summary %>% ungroup() %>% 
#   # separate(Arabinose, into = c('Arabinose_num', 'temp'), sep = ' ') %>%
#   separate(Arabinose, into = c('Arabinose_num', 'temp'), sep = '% ') %>%
#   mutate(Arabinose_num = factor(Arabinose_num, 
#                                 # levels = c('0%', '0.01%', '0.025%', '0.05%', '0.2%', '0.4%'))) %>%
#                                 levels = c('0', '0.01', '0.025', '0.05', '0.2', '0.4')), 
#          conf_interval = 1.96*sem_fluo, 
#          exp_level = factor(exp_level, levels = c('Absent', 'Weak', 'Suboptimal',
#                                                   'Near-optimal', 'Optimal', 'Overexpressed'))) %>%
#   # filter(Time != 'Stationary phase') %>%
#   # filter(Time == 'After 10 generations') %>%
#   ggplot(aes(x = Arabinose_num, y = med_fluo, 
#              colour = exp_level)) + 
#   # facet_wrap(~Time) +
#   # geom_jitter(width = 0.2, size = 3) +
#   # geom_errorbar(aes(ymax = med_fluo + sem_fluo, ymin = med_fluo - sem_fluo), width = 0.1) +
#   geom_pointrange(aes(ymax = med_fluo + conf_interval, ymin = med_fluo - conf_interval), 
#                  position = position_jitter(width = 0.2)) +
#   scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
#   # geom_violin(fill = 'grey', alpha = 0.7) +
#   # geom_boxplot(fill = 'grey', alpha = 0.5, outlier.shape = NA) +
#   xlab('Arabinose concentration (% m/v)') + ylab('GFP fluorescence (log10)') +
#   theme(axis.text = element_text(size = 18),
#         # strip.background = element_rect(fill = 'white'), strip.text = element_text(face = 'bold', size = 16), 
#         axis.title = element_text(size = 20, face = 'bold'), 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'none')
# p_fig1C
# ggsave(p_fig1C, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-06_Figures_paper/Fig1C_ara0.4_medians.pdf')

# p_fig1D <- all_data_processed_expression %>% ungroup() %>% 
#   # filter(Arabinose != '0% arabinose') %>%
#   mutate(
#     exp_level = ifelse(Arabinose == '0.2% arabinose', 'Optimal', 
#                        ifelse(Arabinose == '0.05% arabinose', 'Near-optimal',
#                               ifelse(Arabinose == '0.025% arabinose', 'Suboptimal', 
#                                      ifelse(Arabinose == '0.01% arabinose', 'Weak', 
#                                             ifelse(Arabinose == '0.4% arabinose', 'Overexpressed',
#                                                    'Absent')))))
#   ) %>%
#   separate(Arabinose, into = c('Arabinose_num', 'temp'), sep = '% ') %>%
#   mutate(Arabinose_num = factor(Arabinose_num, 
#                                 levels = c('0', '0.01', 
#                                   '0.025', '0.05', '0.2', '0.4'))) %>%
#   mutate(exp_level = factor(exp_level, levels = c('Absent', 'Weak', 'Suboptimal',
#                                                   'Near-optimal', 'Optimal', 'Overexpressed'))) %>%
#   ggplot(aes(x = Arabinose_num, y = GRNBHLog, 
#              colour = exp_level, group = Well)) + 
#   geom_violin(aes(fill = exp_level), alpha = 0.6, scale = 'width',
#               position = 'dodge', width = 0.8) +
#   geom_boxplot(aes(fill = exp_level), alpha = 0.4, position = 'dodge', width = 0.8) +
#   scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
#   scale_fill_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
#   xlab('Arabinose concentration (% m/v)') + ylab('GFP fluorescence (log10)') +
#   theme(axis.text = element_text(size = 18),
#         axis.title = element_text(size = 20, face = 'bold'), 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'none')
# p_fig1D

#### New version with expression level on the x-axis ####

p_fig1D <- all_data_processed_expression %>% ungroup() %>% 
  # filter(Arabinose != '0% arabinose') %>%
  mutate(
    exp_level = ifelse(Arabinose == '0.2% arabinose', 'Optimal', 
                       ifelse(Arabinose == '0.05% arabinose', 'Near-optimal',
                              ifelse(Arabinose == '0.025% arabinose', 'Suboptimal', 
                                     ifelse(Arabinose == '0.01% arabinose', 'Weak', 
                                            ifelse(Arabinose == '0.4% arabinose', 'Overexpressed',
                                                   'Absent')))))
  ) %>%
  # separate(Arabinose, into = c('Arabinose_num', 'temp'), sep = '% ') %>%
  # mutate(Arabinose_num = factor(Arabinose_num, 
  #                               levels = c('0', '0.01', 
  #                                          '0.025', '0.05', '0.2', '0.4'))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Absent', 'Weak', 'Suboptimal',
                                                  'Near-optimal', 'Optimal', 'Overexpressed'))) %>%
  ggplot(aes(x = exp_level, y = GRNBHLog, 
             colour = exp_level, group = Well)) + 
  geom_violin(aes(fill = exp_level), alpha = 0.6, scale = 'width',
              position = 'dodge', width = 0.8) +
  geom_boxplot(aes(fill = exp_level), alpha = 0.4, position = 'dodge', width = 0.8) +
  scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  scale_fill_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  xlab('Expression level') + ylab('GFP fluorescence (log10)') +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20, face = 'bold'), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'none')
p_fig1D


#### Panels E and F ####

expression_level <- fig1C_data %>% ungroup() %>%
  group_by(Arabinose, exp_level) %>%
  summarise(mean_fluo = mean(mean_fluo)) %>%
  mutate(exp_level = factor(exp_level, levels = c('Absent', 'Weak', 'Suboptimal', 'Near-optimal', 'Optimal', 'Overexpressed')))

# Show the average growth recovery for each class in figure 1D
fig1C_data %>% ungroup() %>% 
  group_by(Arabinose) %>%
  summarise(pct_recovery = mean(pct_recovery))

#### Figure 1E: Scatterplot of fitness effects at different concentrations vs 0.01% arabinose ####

data_fig_1e <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose)

# Separate the data for optimal expression from the rest
data_part_1 <- data_fig_1e %>%
  filter(Arabinose == 0.2)

data_part_2 <- data_fig_1e %>%
  filter(Arabinose != 0.2)

# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "sem_sel_coeff_2", "Arabinose_2")

# Join
data_fig_1e_final <- inner_join(x = data_part_1, y = data_part_2, 
                               by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
)

data_fig_1e_final_exp <- left_join(x = data_fig_1e_final, y = expression_level, 
                                  by = c('Arabinose_2' = 'Arabinose'))

# Make sure WT appears only once
data_fig_1e_final_exp_nowt <- data_fig_1e_final_exp %>% filter(WT_Residue != Residue)
data_fig_1e_final_exp_wt <- data_fig_1e_final_exp %>% filter(WT_Residue == Residue, Position == 2)

data_fig_1e_final_exp_plot <- bind_rows(data_fig_1e_final_exp_nowt, data_fig_1e_final_exp_wt)

# Define a function for npc annotations
## Old annotation box for width = 17
# annotation_box <- roundrectGrob(x = unit(0.225, 'npc'), y = unit(0.755, 'npc'), 
#                            width = unit(0.3, 'npc'), height = unit(0.3, 'npc'), 
#                            gp = gpar(lwd = 2, col = 'black',
#                                      # fill = '#999999', 
#                                      fill = '#737373',
#                                      lty = 2))

## New annotation box for width = 20
annotation_box <- roundrectGrob(x = unit(0.205, 'npc'), y = unit(0.755, 'npc'), 
                                width = unit(0.26, 'npc'), height = unit(0.3, 'npc'), 
                                gp = gpar(lwd = 2, col = 'black',
                                          # fill = '#999999', 
                                          fill = '#737373',
                                          lty = 2))

p_fig1e <- data_fig_1e_final_exp_plot %>%
  ggplot(aes(x = mean_sel_coeff, y = mean_sel_coeff_2, 
             colour = as.factor(exp_level))
  ) +
  geom_point(alpha = 1, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'black') +
  # xlab('s (optimal expression)') + ylab('s (non-optimal expression)') + 
  # labs(y = expression(bold(italic(s) + '(non-optimal expression)'))) +
  labs(colour = 'Expression level in y-axis', 
       x = expression(paste(bolditalic('s '), bold('(optimal expression)'), sep = ' ')), 
       y = expression(paste(bolditalic('s '), bold('(non-optimal expression)'), sep = ' '))
       ) +
  # scale_colour_manual(values = c('#fd8d3c', '#e6002e', '#80001a')) +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', 'black')) +
  # annotate_npc(geom = geom_rect, xmin = 0.03, xmin = 0.3,
  #              ymax = 0.9, ymax = 0.6, fill = 'grey') +
  # geom_rect(xmin = -0.85, xmax = -0.55,
  #           ymin = -0.15, ymax = 0.3,
  #           fill = 'grey', ) +
  # geom_tile(x = 0.05, y = 0.9, width = 0.15, height = 0.15, fill = 'grey', 
  #           inherit.aes = FALSE) +
  annotation_custom(annotation_box) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 6,
           label.x.npc = 0.05, label.y.npc = 0.9,
           method = 'spearman', show.legend = F, cor.coef.name = 'rho'
  ) +
  geom_smooth(method = 'loess', show.legend = F) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'top', 
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()
  ) + xlim(-1, 0.5) + ylim(-1, 0.5) +
  guides(fill = 'none')
p_fig1e  

#### Figure 1F: Paired fitness effects per mutants ####
all_data_all_reps <- read_delim('../R67_DMS_February2022/Complete_datasets_MiSeq_NovaSeq/all_data_all_reps_bothSequencers.txt', 
                                delim = '\t')

all_data_all_reps_TMP10 <- all_data_all_reps %>% rowwise() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  # Remove WT and stop codons (keep only missense mutants)
  filter(WT_Residue != Residue, Residue != '*') %>%
  mutate(Genotype = str_c(WT_Residue, Residue, Position)) %>%
  select(ID, Genotype, Timepoint, TMP, Arabinose, sel_coeff, Sequencer)

genotypes <- unique(all_data_all_reps_TMP10$Genotype)
all_anovas <- c()

## Do an ANOVA of each mutant
for(genotype in genotypes){
  ## Compile a table with each mutant and its p-value from the ANOVA
  # Select data from one genotype to do the ANOVA
  data_genotype <- all_data_all_reps_TMP10 %>% filter(Genotype == genotype) %>% select(-Timepoint, -TMP) %>%
    mutate(Arabinose = as.factor(Arabinose), 
           ID = as.factor(ID)) %>% group_by(Arabinose) %>%
    mutate(Replicate = as.factor(row_number()))
  
  #### ANOVA on ranks ####
  
  ## Do the aligned rank transform and the ANOVA
  # m <- art(sel_coeff ~ Arabinose + Error(Sample), data=data_genotype)
  # summ_m <- summary(m)
  
  # Make sure that the column sums of aligned responses are close to zero
  # if(summ_m$aligned.col.sums[[1]] > 0.001){
  #   print(str_c('Aligned column sums for ', genotype, ' = ', toString(summ_m$aligned.col.sums[[1]])))
  # }
  
  # anova_m <- anova(m)
  # p_value <- anova_m$`Pr(>F)`
  
  #### Regular ANOVA ####
  m <- aov(sel_coeff ~ Arabinose, data=data_genotype)
  anova_test <- anova(m)
  p_value <- anova_test$`Pr(>F)`[1]
  
  new_row <- data.frame(Genotype = c(genotype), p_val = c(p_value))
  all_anovas <- bind_rows(all_anovas, new_row)
}

## Apply the Benjamini-Hochberg correction for multiple hypotheses
fdr_test <-  p.fdr(pvalues = all_anovas$p_val, adjust.method = 'BH', threshold = 0.05)
all_anovas$p.adj <- fdr_test$fdrs

summary(all_anovas$p.adj)
sum(all_anovas$p.adj < 0.05)

# Prepare the data
temp_data <- data_fig_1e %>% rowwise() %>%
  mutate(ID = str_c(Position, Residue, sep = ''))

data_fig1f <- data_fig_1e_final %>% ungroup() %>%
 mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff,
        Arabinose_2 = str_c(Arabinose_2, '% arabinose', sep = ''))

# Relevel the data set to make sure the boxplots appear in the proper order
data_fig1f %<>% mutate(Arabinose_2 = factor(Arabinose_2, 
                                           levels = c('0.01% arabinose', '0.025% arabinose',
                                                      '0.05% arabinose', '0.4% arabinose'))
)

data_fig1f_exp <- left_join(x = data_fig1f %>% 
                             separate(col = Arabinose_2, into = c('Arabinose_num', 'tmp'), sep = '% arabinose') %>%
                             mutate(Arabinose_num = as.numeric(Arabinose_num)),
                           y = expression_level, 
                           by = c('Arabinose_num' = 'Arabinose')) %>%
 mutate(exp_level = as.factor(exp_level))

# Let's add a column with the difference in scores at low and high expression
score_diff <- data_fig1f %>% 
  filter(Arabinose_2 == '0.01% arabinose') %>%
  mutate(mut_check_diff = ifelse(abs(diffNormScore) > 0.3, TRUE, FALSE))

score_significant <- left_join(x = data_fig1f_exp %>% filter(Arabinose_num == 0.01), 
                               y = all_anovas %>% rowwise() %>%
                                 separate(Genotype, into = c('WT_Residue', 'Residue', 'Position'), 
                                          sep = c(1, 2)) %>%
                                 mutate(Position = as.numeric(Position)), 
                               by = c('WT_Residue' = 'WT_Residue', 'Position' = 'Position', 
                                      'Residue' = 'Residue')) %>%
  mutate(mut_check = p.adj < 0.05)

table(score_diff$mut_check_diff)
table(score_significant$mut_check)

# Join the data to add these marks
data_fig1f <- left_join(x = temp_data, 
                        y= score_significant %>% select(Position, Residue, p.adj, mut_check),
                        by = c('Position' = 'Position', 'Residue' = 'Residue'))


# Add the expression level
data_fig1f_exp <- left_join(x = data_fig1f, y = expression_level, 
                            by = c('Arabinose' = 'Arabinose')) %>%
  mutate(exp_level = as.factor(exp_level))

# Add the data about the minimum change in expression
data_fig1f_exp <- left_join(x = data_fig1f_exp, 
                            y = score_diff %>% select(Position, Residue, diffNormScore, mut_check_diff), 
                            by = c('Position' = 'Position', 'Residue' = 'Residue')) %>%
  rowwise() %>%
  mutate(mut_check_final = and(mut_check, mut_check_diff))

## Figure for all mutants
p_fig1f <- data_fig1f_exp %>% 
  filter(Arabinose %in% c(0.01, 0.025, 0.05, 0.2, 0.4)) %>%
  mutate(Arabinose = factor(Arabinose, levels = c(0.01, 0.025, 0.05, 0.2, 0.4))) %>%
  ggplot(aes(x = exp_level, y = mean_sel_coeff, fill = exp_level)) +
  geom_violin(alpha = 0.4) +
  stat_summary(fun="mean", geom="point", size = 5, colour = 'black') +
  geom_point() +
  geom_line(aes(group = ID, colour = mut_check, alpha = mut_check_final)) +
  scale_alpha_manual(values = c(0.1, 0.4)) +
  scale_colour_manual(values = c('grey', 'black')) +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'none', 
        axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.justification = 0.5) +
  xlab('Expression level') +
  ylab('s')
p_fig1f

## Draw the panels A and B
# p_fig1A <- ggdraw() + draw_image('Figures/2022-04-06_Figures_paper/Fig1A_DfrB1.png')
# p_fig1B <- ggdraw() + draw_image('Figures/2022-04-06_Figures_paper/Fig1B_DfrB1.png')

# p_fig1A <- ggdraw() + draw_image('Figures/2022-05-06_Figures_paper/Fig1A_DfrB1.png')
# p_fig1B <- ggdraw() + draw_image('Figures/2022-05-06_Figures_paper/Fig1B_DfrB1.png')
# 
# # Draw figure 1
# p_fig1AB <- plot_grid(p_fig1A, p_fig1B, ncol = 2, 
#                       labels = c('A', 'B'), label_fontface = 'bold', 
#                       label_size = 20)
# legend_fig1CD <- get_legend(p_fig1C + 
#                               theme(legend.title = element_text(size = 18), 
#                                     legend.text = element_text(size = 18)))
# p_fig1CD <- plot_grid(p_fig1C + theme(legend.position = 'none'), p_fig1D, ncol = 2,
#                       labels = c('C', 'D'), label_fontface = 'bold', 
#                       label_size = 20)
# 
# p_fig1EF <- plot_grid(p_fig1e + theme(legend.position = 'none'),
#                       p_fig1f, labels = c('E', 'F'), label_fontface = 'bold', 
#                       label_size = 20)
# 
# p_fig1 <- plot_grid(p_fig1AB, legend_fig1CD, p_fig1CD, p_fig1EF, nrow = 4, 
#                     rel_heights = c(0.5, 0.1, 1, 1))
# p_fig1
# ggsave(p_fig1, device = cairo_pdf, width = 17, height = 14, dpi = 300, 
#        filename = 'Figures/2022-05-06_Figures_paper/1.Fig1_ara0.4_new.pdf')
# ggsave(p_fig1, device = 'png', width = 17, height = 14, dpi = 300, 
#        filename = 'Figures/2022-05-06_Figures_paper/1.Fig1_ara0.4_new.png')

#### Try using k-means to cluster the data from Figure 1F ####

data_fig1f_wide <- data_fig1f_exp %>% ungroup() %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, exp_level) %>%
  group_by(Position, WT_Residue, Residue) %>%
  pivot_wider(names_from = exp_level, values_from = mean_sel_coeff)

## Make sure WT appears only once
wt_row <- data_fig1f_wide %>% filter(WT_Residue == Residue, Position == 2)

data_kmeans <- bind_rows(data_fig1f_wide %>% filter(WT_Residue != Residue), 
                         wt_row)

set.seed(100)
kmeans_ss <- c()
for(k in 1:10){
  print(k)
  kmeans_fitness <- kmeans(data_kmeans %>% ungroup() %>% rowwise() %>% 
                             filter(WT_Residue != Residue) %>%
                             # mutate(ID = str_c(WT_Residue, Position, Residue)) %>%
                             select(-WT_Residue, -Position, -Residue), 
                           centers = k, iter.max = 10, nstart = 25
                           )
  print('kmeans complete')
  new_row <- c(k, kmeans_fitness$tot.withinss)
  print('Row complete')
  
  kmeans_ss <- rbind(kmeans_ss, new_row)
}

kmeans_ss %<>% as.data.frame()
colnames(kmeans_ss) <- c('k', 'tot_within_ss')

p <- kmeans_ss %>% 
  ggplot(aes(x = k, y = tot_within_ss)) +
  geom_point() + 
  geom_line() +
  scale_x_continuous(breaks = 1:10) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab('Number of clusters (k)') + ylab('Sum of squared errors')
p

## Based on the figure above, let's focus on k = 2 and k = 3 (k = 3 is probably the best)
k = 4
kmeans_fitness <- kmeans(data_kmeans %>% ungroup() %>% rowwise() %>% 
                           # mutate(ID = str_c(WT_Residue, Position, Residue)) %>%
                           select(-WT_Residue, -Position, -Residue), 
                         centers = k, iter.max = 10, nstart = 25
)

kmeans_fitness$centers

# Rearrange the clusters automatically
temp_clusters <- kmeans_fitness$centers
temp_clusters %<>% as.data.frame() %>% mutate(old_cluster = row_number())

cluster_relabel <- temp_clusters %>% arrange(Weak) %>% mutate(new_cluster = row_number())

# Run a principal component analysis
pca_fitness <- princomp(x = data_kmeans %>% ungroup() %>%
                          select(-WT_Residue, -Position, -Residue))

# Save the centers and the loadings in variables
pca_loadings <- pca_fitness$loadings[1:5, 1:2]
pca_centers <- pca_fitness$center
print(pca_fitness$sdev)

pca_coords_points <- cbind(pca_fitness$scores, kmeans_fitness$cluster)
pca_coords_points <- as.data.frame(pca_coords_points)
colnames(pca_coords_points) <- c('PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'Cluster')
# Use the new numbering for the clusters
pca_coords_points <- left_join(x = pca_coords_points, 
                               y = cluster_relabel %>% ungroup() %>%
                                 select(old_cluster, new_cluster), 
                               by = c('Cluster' = 'old_cluster')
                               )

pca_coords_points %<>% mutate(Cluster = new_cluster) %>%
  select(-new_cluster)

# table(pca_coords_points$Cluster)

# Need to apply the transformation to the kmeans centers

new_pca_centers <- c()
for(i in 1:k){
  new_pca_centers <- rbind(new_pca_centers, pca_centers)
}

kmeans_centers_new <- kmeans_fitness$centers - new_pca_centers
kmeans_centers_projected <- kmeans_centers_new %*% pca_loadings
kmeans_centers_projected <- as.data.frame(kmeans_centers_projected) %>%
  mutate(Cluster = 1:k)
colnames(kmeans_centers_projected) <- c('PC1', 'PC2', 'old_cluster')

kmeans_centers_projected <- left_join(x = kmeans_centers_projected, 
                                      y = cluster_relabel %>% ungroup() %>%
                                        select(old_cluster, new_cluster), 
                                      by = c('old_cluster' = 'old_cluster')) %>%
  mutate(Cluster = new_cluster) %>% select(-new_cluster)

p <- pca_coords_points %>% 
  ggplot(aes(x = PC1, y = PC2, colour = as.factor(Cluster))) +
  geom_point(alpha = 0.3) +
  geom_point(data = kmeans_centers_projected, size = 5, shape = 24, 
             colour = 'black', aes(fill = as.factor(Cluster))) +
  theme(legend.position = 'top') +
  xlab('PC1 (68.14%)') + ylab('PC2 (15.31%)') +
  labs(colour = 'Cluster', fill = 'Cluster')
p  

#### Try a different version of Fig 1F indicating the cluster numbers ####

data_fig1f_new <- data_kmeans %>% ungroup() %>% 
  mutate(cluster = kmeans_fitness$cluster) %>%
  pivot_longer(cols = c('Weak', 'Suboptimal', 'Near-optimal', 
                        'Optimal', 'Overexpressed'), names_to = 'exp_level', 
               values_to = 'mean_sel_coeff')

data_fig1f_new <- left_join(x = data_fig1f_new, 
                            y = cluster_relabel %>% ungroup() %>%
                              select(old_cluster, new_cluster), 
                            by = c('cluster' = 'old_cluster')
                            ) %>%
  mutate(cluster = new_cluster) %>% select(-new_cluster)

p_fig1f <- data_fig1f_new %>% rowwise() %>%
  mutate(ID = str_c(WT_Residue, Position, Residue, sep = ''), 
         exp_level = factor(exp_level,
                            levels = c('Weak', 'Suboptimal', 'Near-optimal',
                                       'Optimal', 'Overexpressed'))) %>%
  # filter(Arabinose %in% c(0.01, 0.025, 0.05, 0.2, 0.4)) %>%
  # mutate(Arabinose = factor(Arabinose, levels = c(0.01, 0.025, 0.05, 0.2, 0.4))) %>%
  ggplot(aes(x = exp_level, y = mean_sel_coeff, fill = exp_level)) +
  geom_violin(alpha = 0.4) +
  stat_summary(fun="mean", geom="point", size = 5, colour = 'black') +
  geom_point() +
  geom_line(aes(group = ID, colour = as.factor(cluster)), alpha = 0.3) +
  # scale_alpha_manual(values = c(0.1, 0.4)) +
  # scale_colour_manual(values = c('grey', 'black')) +
  # scale_colour_manual(values = c('#1b9e77', '#ffbf80', '#7570b3', '#e7298a', '#66a61e')) +
  scale_colour_manual(
    values = c('#1b9e77', '#ffbf80', '#7570b3', '#e7298a')
    # values = c('#a6611a', '#80cdc1', '#018571', '#b3b3b3', '#dfc27d')
    ) +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'), 
                    guide = 'none') +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', 
        axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.justification = 0.5) +
  xlab('Expression level') +
  # ylab('s') + 
  labs(colour = 'Cluster', y = expression(bolditalic('s'))) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
p_fig1f

#### Draw figure 1 ####
# 
# p_fig1A <- ggdraw() + draw_image('Figures/2022-05-06_Figures_paper/Fig1A_DfrB1.png')
# p_fig1B <- ggdraw() + draw_image('Figures/2022-05-06_Figures_paper/Fig1B_DfrB1.png')
# 
# p_fig1AB <- plot_grid(p_fig1A, p_fig1B, ncol = 2, 
#                       labels = c('A', 'B'), label_fontface = 'bold', 
#                       label_size = 20)
# 
# legend_fig1CD <- get_legend(p_fig1C + 
#                               theme(legend.title = element_text(size = 18), 
#                                     legend.text = element_text(size = 18))
#                             )
# 
# p_fig1CD <- plot_grid(p_fig1C + theme(legend.position = 'none'), p_fig1D, ncol = 2,
#                       labels = c('C', 'D'), label_fontface = 'bold', 
#                       label_size = 20)
# 
# p_fig1EF <- plot_grid(p_fig1e + theme(legend.position = 'none'),
#                       p_fig1f, labels = c('E', 'F'), label_fontface = 'bold', 
#                       label_size = 20)
# 
# p_fig1 <- plot_grid(p_fig1AB, legend_fig1CD, p_fig1CD, p_fig1EF, nrow = 4, 
#                     rel_heights = c(0.5, 0.1, 1, 1))
# p_fig1
# ggsave(p_fig1, device = cairo_pdf, width = 17, height = 14, dpi = 300, 
#        filename = 'Figures/2022-05-06_Figures_paper/1.Fig1_ara0.4_new.pdf')
# ggsave(p_fig1, device = 'png', width = 17, height = 14, dpi = 300, 
#        filename = 'Figures/2022-05-06_Figures_paper/1.Fig1_ara0.4_new.png')

#### New version of Fig 1 ####

p_fig1A <- ggdraw() + draw_image('Figures/2022-05-06_Figures_paper/Fig1A_DfrB1.png')
p_fig1B <- ggdraw() + draw_image('Figures/2022-05-06_Figures_paper/Fig1B_DfrB1.png')

p_fig1AB <- plot_grid(p_fig1A, p_fig1B, ncol = 2, 
                      labels = c('A', 'B'), label_fontface = 'bold', 
                      label_size = 20)


p_fig1CD <- plot_grid(p_fig1C, p_fig1D, ncol = 2,
                      labels = c('C', 'D'), label_fontface = 'bold', 
                      label_size = 20)

p_fig1EF <- plot_grid(p_fig1e + theme(legend.position = 'none'),
                      p_fig1f, labels = c('E', 'F'), label_fontface = 'bold', 
                      label_size = 20)

p_fig1 <- plot_grid(p_fig1AB, p_fig1CD, p_fig1EF, nrow = 3, 
                    rel_heights = c(0.5, 1, 1))
p_fig1
ggsave(p_fig1, device = cairo_pdf, width = 20, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-06_Figures_paper/1.Fig1_ara0.4_new.pdf')
ggsave(p_fig1, device = 'png', width = 20, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-06_Figures_paper/1.Fig1_ara0.4_new.png')


#### This section has all been moved to the script for supplementary figures ####
# #### Let's check if the k-means clusters are different in terms of sites and ddGs ####
#  
# clusters_ddG <- left_join(x = all_data_complete, 
#                           y = data_fig1f_new %>% ungroup() %>%
#                             select(Position, WT_Residue, Residue, cluster) %>%
#                             unique(), 
#                           by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue',
#                                  'Residue' = 'Residue')) %>%
#   group_by(Position, WT_Residue, Residue, cluster) %>%
#   filter(!(is.na(cluster)), !(is.na(Mean_ddG_stab_HET))) %>%
#   summarise(ddG_stab= mean(Mean_ddG_stab_HET), 
#             ddG_dim_int = mean(Mean_ddG_int_HM_A_C), 
#             ddG_tet_int = mean(Mean_ddG_int_HM_A_D))
# 
# p <- clusters_ddG %>% ungroup() %>%
#   pivot_longer(cols = c('ddG_stab', 'ddG_dim_int', 'ddG_tet_int'), 
#                names_to = 'ddG_type', values_to = 'ddG_value') %>%
#   mutate(ddG_type = ifelse(ddG_type == 'ddG_stab', 'Subunit stability', 
#                            ifelse(ddG_type == 'ddG_dim_int',
#                                   'Dimerization interface', 
#                                   ifelse(ddG_type == 'ddG_tet_int', 
#                                          'Tetramerization interface', NA)))) %>%
#   ggplot(aes(x = as.factor(cluster), y = ddG_value, fill = ddG_type)) +
#   geom_violin(aes(fill = ddG_type), alpha = 0.6, scale = 'width',
#               position = 'dodge', width = 0.8) +
#   geom_boxplot(aes(fill = ddG_type), alpha = 0.4, position = 'dodge', width = 0.8) +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#                 panel.grid.minor = element_blank(),
#                 legend.position = 'top',
#                 axis.title = element_text(face = 'bold', size = 20),
#                 axis.text = element_text(size = 18),
#                 legend.text = element_text(size = 20),
#                 legend.title = element_text(size = 20),
#                 legend.justification = 0.5
#         ) +
#   ylim(-5, 10) +
#   xlab('Cluster') + ylab('\u0394\u0394G [kcal / mol]') + 
#   labs(fill = '')
# p
# ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/Fig_extra_ddG_clusters.pdf')
# 
# #### Check if the clusters belong preferentially to specific sites ####
# 
# # Load the table with protein sites
# data_interfaces_final <- read_delim('Data/Complete_datasets_TMP0_TMP10/data_annotation_2.txt', delim = '\t')
# 
# data_fig1f_new_summ <- data_fig1f_new %>% ungroup() %>% 
#   group_by(Position, WT_Residue, Residue) %>%
#   summarise(Cluster = mean(cluster))
# 
# # Add the data for protein regions
# data_fig1f_new_summ <- left_join(x = data_fig1f_new_summ, 
#                                  y = data_interfaces_final %>% rowwise() %>%
#                                    mutate(Unannotated = ifelse(sum(`A,C`, `A,D`, DHF, NADPH,
#                                                                    Cat_residues, Disordered_region,
#                                                                    Buried) == 0, 
#                                                                1, 0)
#                                           ), 
#                                  by = c('Position' = 'Position'))
# 
# ## A global analysis
# data_fig1f_new_glob <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site_check == 1)
# 
# cont_table <- table(data_fig1f_new_glob$Site, data_fig1f_new_glob$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Prepare a figure for the residuals
# residuals <- as.data.frame(chisq$residuals) %>%
#   pivot_wider(names_from = Var2, values_from = Freq)
# 
# colnames(residuals) <- c('Site', 'Cluster 1', 'Cluster 2',
#                          'Cluster 3', 'Cluster 4', 'Cluster 5')
# 
# needed_rownames <- residuals$Site
# residuals %<>% select(-Site)
# 
# residuals <- as.matrix(residuals)
# rownames(residuals) <- needed_rownames
# 
# heatmap_test <- Heatmap(residuals, cluster_columns = F, cluster_rows = F, 
#                         ## Color ramp for the old normalized scores
#                         # col = colorRamp2(
#                         #   breaks = seq(0, 4, length.out = 7),
#                         #   colors = viridis::viridis(7)),
#                         show_column_names = T, row_names_side = 'left',
#                         # width=unit(31, 'cm'), height = unit(31, 'cm'),
#                         width = unit(10, 'cm'), height = unit(7, 'cm'),
#                         border = T,
#                         row_title = 'Site',
#                         row_title_gp = gpar(fontsize=20, fontface = 'bold'),
#                         row_names_rot = 0, 
#                         row_names_centered = T,
#                         row_names_gp = gpar(fontsize=16, fontface = 'bold'),
#                         column_title = 'Cluster', 
#                         column_title_side = 'bottom',
#                         column_names_gp = gpar(fontsize=16,fontface='bold'),
#                         column_title_gp = gpar(fontsize=20, fontface = 'bold'),
#                         # right_annotation = ha_empty,
#                         show_heatmap_legend = TRUE
#                         # heatmap_legend_param = list(
#                         #   at = c(0, 2, 4),
#                         #   title = "Mutant fraction in library (%)", 
#                         #   title_gp = gpar(fontsize = 20),
#                         #   legend_height = unit(3.5, "cm"),
#                         #   legend_width = unit(2, "cm"),
#                         #   border='black',
#                         #   lwd=1.7,
#                         #   labels_gp = gpar(fontsize = 18),
#                         #   title_position = "leftcenter-rot"
#                         # )
# )
# heatmap_test
# 
# ## Only for the dimerization interface
# data_fig1f_new_dim_int <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'A,C')
# 
# cont_table <- table(data_fig1f_new_dim_int$Site_check, data_fig1f_new_dim_int$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# dim_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# ## Only for the tetramerization interface
# data_fig1f_new_tet <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'A,D')
# 
# cont_table <- table(data_fig1f_new_tet$Site_check, data_fig1f_new_tet$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# tet_log2fold <- log2((chisq$observed + 1)/ chisq$expected)[2,]
# 
# ## Only for DHF binding
# data_fig1f_dhf <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'DHF')
# 
# cont_table <- table(data_fig1f_dhf$Site_check, data_fig1f_dhf$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# dhf_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# ## Only for the catalytic residues
# data_fig1f_new_cat <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'Cat_residues')
# 
# cont_table <- table(data_fig1f_new_cat$Site_check, data_fig1f_new_cat$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# cat_log2fold <- log2((chisq$observed + 1)/ chisq$expected)[2,]
# 
# ## Only for the NADPH binding residues
# data_fig1f_new_nadph <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'NADPH')
# 
# cont_table <- table(data_fig1f_new_nadph$Site_check, data_fig1f_new_nadph$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# nadph_log2fold <- log2((chisq$observed + 1)/ chisq$expected)[2,]
# 
# ## Only for the disordered region
# data_fig1f_new_dis <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'Disordered_region')
# 
# cont_table <- table(data_fig1f_new_dis$Site_check, data_fig1f_new_dis$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# dis_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# # Prepare a figure for the residuals
# residuals <- as.data.frame(chisq$residuals) %>%
#   pivot_wider(names_from = Var2, values_from = Freq)
# 
# colnames(residuals) <- c('Disordered region', 'Cluster 1', 'Cluster 2',
#                          'Cluster 3', 'Cluster 4', 'Cluster 5')
# 
# needed_rownames <- residuals$`Disordered region`
# residuals %<>% select(-`Disordered region`)
# 
# residuals <- as.matrix(residuals)
# rownames(residuals) <- needed_rownames
# 
# heatmap_test <- Heatmap(residuals, cluster_columns = F, cluster_rows = F, 
#                         ## Color ramp for the old normalized scores
#                         # col = colorRamp2(
#                         #   breaks = seq(0, 4, length.out = 7),
#                         #   colors = viridis::viridis(7)),
#                         show_column_names = T, row_names_side = 'left',
#                         # width=unit(31, 'cm'), height = unit(31, 'cm'),
#                         width = unit(10, 'cm'), height = unit(7, 'cm'),
#                         border = T,
#                         row_title = 'Disordered region',
#                         row_title_gp = gpar(fontsize=20, fontface = 'bold'),
#                         row_names_rot = 0, 
#                         row_names_centered = T,
#                         row_names_gp = gpar(fontsize=16, fontface = 'bold'),
#                         column_title = 'Cluster', 
#                         column_title_side = 'bottom',
#                         column_names_gp = gpar(fontsize=16,fontface='bold'),
#                         column_title_gp = gpar(fontsize=20, fontface = 'bold'),
#                         # right_annotation = ha_empty,
#                         show_heatmap_legend = TRUE
#                         # heatmap_legend_param = list(
#                         #   at = c(0, 2, 4),
#                         #   title = "Mutant fraction in library (%)", 
#                         #   title_gp = gpar(fontsize = 20),
#                         #   legend_height = unit(3.5, "cm"),
#                         #   legend_width = unit(2, "cm"),
#                         #   border='black',
#                         #   lwd=1.7,
#                         #   labels_gp = gpar(fontsize = 18),
#                         #   title_position = "leftcenter-rot"
#                         # )
# )
# heatmap_test
# 
# ## Only for the buried residues
# data_fig1f_new_bur <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'Buried')
# 
# cont_table <- table(data_fig1f_new_bur$Site_check, data_fig1f_new_bur$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# bur_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# ## Only for the unannotated residues
# data_fig1f_new_unannot <- data_fig1f_new_summ  %>% 
#   pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
#                         'Cat_residues', 'Disordered_region', 'Buried', 
#                         'Unannotated'), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'Unannotated')
# 
# cont_table <- table(data_fig1f_new_unannot$Site_check, data_fig1f_new_unannot$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# unannot_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# # Prepare a figure for the residuals
# residuals <- as.data.frame(chisq$residuals) %>%
#   pivot_wider(names_from = Var2, values_from = Freq)
# 
# colnames(residuals) <- c('Unannotated', 'Cluster 1', 'Cluster 2',
#                          'Cluster 3', 'Cluster 4', 'Cluster 5')
# 
# needed_rownames <- residuals$Unannotated
# residuals %<>% select(-Unannotated)
# 
# residuals <- as.matrix(residuals)
# rownames(residuals) <- needed_rownames
# 
# heatmap_test <- Heatmap(residuals, cluster_columns = F, cluster_rows = F, 
#                         ## Color ramp for the old normalized scores
#                         # col = colorRamp2(
#                         #   breaks = seq(0, 4, length.out = 7),
#                         #   colors = viridis::viridis(7)),
#                         show_column_names = T, row_names_side = 'left',
#                         # width=unit(31, 'cm'), height = unit(31, 'cm'),
#                         width = unit(10, 'cm'), height = unit(7, 'cm'),
#                         border = T,
#                         row_title = 'Unannotated',
#                         row_title_gp = gpar(fontsize=20, fontface = 'bold'),
#                         row_names_rot = 0, 
#                         row_names_centered = T,
#                         row_names_gp = gpar(fontsize=16, fontface = 'bold'),
#                         column_title = 'Cluster', 
#                         column_title_side = 'bottom',
#                         column_names_gp = gpar(fontsize=16,fontface='bold'),
#                         column_title_gp = gpar(fontsize=20, fontface = 'bold'),
#                         # right_annotation = ha_empty,
#                         show_heatmap_legend = TRUE
#                         # heatmap_legend_param = list(
#                         #   at = c(0, 2, 4),
#                         #   title = "Mutant fraction in library (%)", 
#                         #   title_gp = gpar(fontsize = 20),
#                         #   legend_height = unit(3.5, "cm"),
#                         #   legend_width = unit(2, "cm"),
#                         #   border='black',
#                         #   lwd=1.7,
#                         #   labels_gp = gpar(fontsize = 18),
#                         #   title_position = "leftcenter-rot"
#                         # )
# )
# heatmap_test
# 
# #### Add the ddG data ####
# data_fig1f_new_summ_ddg <- left_join(x = data_fig1f_new_summ, 
#                  y = all_data_complete %>% select(Position, WT_Residue, Residue,
#                                                   Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, 
#                                                   Mean_ddG_int_HM_A_D) %>%
#                    unique(),
#                  by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
#                         'Residue' = 'Residue'))
# 
# data_fig1f_new_summ_ddg %<>% filter(!(is.na(Mean_ddG_stab_HET))) %>%
#   mutate(high_ddG_stab = ifelse(Mean_ddG_stab_HET >= 2, 1, 0), 
#          low_ddG_stab = ifelse(Mean_ddG_stab_HET < 2, 1, 0), 
#          
#          high_ddG_dim = ifelse(Mean_ddG_int_HM_A_C >= 2, 1, 0), 
#          low_ddG_dim = ifelse(Mean_ddG_int_HM_A_C < 2, 1, 0), 
#          
#          high_ddG_tet = ifelse(Mean_ddG_int_HM_A_D >= 2, 1, 0), 
#          low_ddG_tet = ifelse(Mean_ddG_int_HM_A_D < 2, 1, 0))
# 
# ## Only for the mutants with high ddG stability
# data_fig1f_high_ddg_stab <- data_fig1f_new_summ_ddg  %>% 
#   pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
#                         high_ddG_dim, low_ddG_dim, 
#                         high_ddG_tet, low_ddG_tet), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'high_ddG_stab')
# 
# cont_table <- table(data_fig1f_high_ddg_stab$Site_check, data_fig1f_high_ddg_stab$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# high_ddg_stab_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# ## Only for the mutants with low ddG stability
# data_fig1f_low_ddg_stab <- data_fig1f_new_summ_ddg  %>% 
#   pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
#                         high_ddG_dim, low_ddG_dim, 
#                         high_ddG_tet, low_ddG_tet), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'low_ddG_stab')
# 
# cont_table <- table(data_fig1f_low_ddg_stab$Site_check, data_fig1f_low_ddg_stab$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# low_ddg_stab_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# ## Only for the mutants with high ddG dimerization interface
# data_fig1f_high_ddg_dim <- data_fig1f_new_summ_ddg  %>% 
#   pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
#                         high_ddG_dim, low_ddG_dim, 
#                         high_ddG_tet, low_ddG_tet), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'high_ddG_dim')
# 
# cont_table <- table(data_fig1f_high_ddg_dim$Site_check, data_fig1f_high_ddg_dim$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# high_ddg_dim_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# ## Only for the mutants with low ddG stability
# data_fig1f_low_ddg_dim <- data_fig1f_new_summ_ddg  %>% 
#   pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
#                         high_ddG_dim, low_ddG_dim, 
#                         high_ddG_tet, low_ddG_tet), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'low_ddG_dim')
# 
# cont_table <- table(data_fig1f_low_ddg_dim$Site_check, data_fig1f_low_ddg_dim$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# low_ddg_dim_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# ## Only for the mutants with high ddG tetramerization interface
# data_fig1f_high_ddg_tet <- data_fig1f_new_summ_ddg  %>% 
#   pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
#                         high_ddG_dim, low_ddG_dim, 
#                         high_ddG_tet, low_ddG_tet), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'high_ddG_tet')
# 
# cont_table <- table(data_fig1f_high_ddg_tet$Site_check, data_fig1f_high_ddg_tet$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# high_ddg_tet_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# ## Only for the mutants with low ddG stability
# data_fig1f_low_ddg_tet <- data_fig1f_new_summ_ddg  %>% 
#   pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
#                         high_ddG_dim, low_ddG_dim, 
#                         high_ddG_tet, low_ddG_tet), names_to = 'Site',
#                values_to = 'Site_check') %>%
#   filter(Site == 'low_ddG_tet')
# 
# cont_table <- table(data_fig1f_low_ddg_tet$Site_check, data_fig1f_low_ddg_tet$Cluster)
# chisq <- chisq.test(cont_table)
# chisq$observed
# chisq$expected
# chisq$p.value
# 
# # Calculate log2(obs / exp)
# low_ddg_tet_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]
# 
# 
# #### Concatenate the data for all sites and put it in a heatmap ####
# 
# # c('A,C', 'A,D', 'DHF', 'NADPH',
# #   'Cat_residues', 'Disordered_region', 'Buried', 
# #   'Unannotated')
# data_enrichment <- rbind(dim_log2fold, tet_log2fold, dhf_log2fold, cat_log2fold, 
#                          nadph_log2fold, dis_log2fold, bur_log2fold, unannot_log2fold, 
#                          high_ddg_dim_log2fold, high_ddg_tet_log2fold, high_ddg_stab_log2fold, 
#                          low_ddg_dim_log2fold, low_ddg_tet_log2fold, low_ddg_stab_log2fold)
# 
# rownames(data_enrichment) <- c('Dimerization interface', 'Tetramerization interface', 
#                                'DHF binding', 'Catalytic residues', 'NADPH binding', 
#                                'Disordered region',
#                                'Buried residues', 'Unannotated residues', 
#                                '\u0394\u0394G dim >= 2', '\u0394\u0394G stab >= 2',
#                                '\u0394\u0394G tet >= 2', 
#                                '\u0394\u0394G dim < 2', '\u0394\u0394G stab < 2',
#                                '\u0394\u0394G tet < 2')
# 
# seq1 <- seq(-8, 0, length.out = 4)
# seq2 <- seq(0, 2, length.out = 4)
# 
# heatmap_test <- Heatmap(data_enrichment, cluster_columns = F, cluster_rows = T,
#                         clustering_distance_rows = 'pearson',
#                         ## Color ramp for the old normalized scores
#                         col = colorRamp2(
#                           breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
#                           colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
#                         show_column_names = T, row_names_side = 'left',
#                         # width=unit(31, 'cm'), height = unit(31, 'cm'),
#                         width = unit(14, 'cm'), height = unit(7, 'cm'),
#                         border = T,
#                         row_title = '',
#                         row_title_gp = gpar(fontsize=18, fontface = 'bold'),
#                         row_names_rot = 0,
#                         row_names_centered = T,
#                         row_names_gp = gpar(fontsize=14, fontface = 'bold'),
#                         column_title = 'Cluster',
#                         column_title_side = 'bottom',
#                         column_names_rot = 0,
#                         column_names_gp = gpar(fontsize=14,fontface='bold'),
#                         column_title_gp = gpar(fontsize=18, fontface = 'bold'),
#                         # right_annotation = ha_empty,
#                         show_heatmap_legend = TRUE,
#                         heatmap_legend_param = list(
#                           at = c(-8, -4, 0, 1, 2),
#                           title = "log2 fold enrichment",
#                           title_gp = gpar(fontsize = 20),
#                           legend_height = unit(3.5, "cm"),
#                           legend_width = unit(2, "cm"),
#                           border='black',
#                           lwd=1.7,
#                           labels_gp = gpar(fontsize = 18),
#                           title_position = "leftcenter-rot"
#                         )
# )
# heatmap_test
# 
# #### Another panel for the kmeans centroids ####
# 
# cluster_relabel_new <- cluster_relabel %>% select(-old_cluster) %>%
#   pivot_longer(cols = c('Weak', 'Suboptimal', 'Near-optimal', 'Optimal', 'Overexpressed'), 
#                names_to = 'exp_level', values_to = 'mean_sel_coeff')
# 
# p <- cluster_relabel_new %>% 
#   mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
#                                                   'Optimal', 'Overexpressed'))) %>%
#   ggplot(aes(x = exp_level, y = mean_sel_coeff, colour = as.factor(new_cluster))) +
#   geom_point(size = 3) +
#   geom_line(aes(group = as.factor(new_cluster)), size = 2) +
#   scale_colour_manual(values = c('#1b9e77', '#ffbf80', '#7570b3', '#e7298a', '#66a61e')) +
#   theme(axis.title = element_text(face = 'bold', size = 20), 
#                 axis.text = element_text(size = 18),
#                 legend.text = element_text(size = 20),
#                 legend.title = element_text(size = 20),
#                 legend.position = 'top',
#                 legend.justification = 0.5,
#                 panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#                 panel.grid.minor = element_blank()) +
#   xlab('Expression level') + ylab('s') + labs(colour = 'Cluster')
# p
# p_fig_supp_new <- plot_grid(grid.grabExpr(draw(heatmap_test)), p, nrow = 2, 
#                             labels = c('A', 'B'), label_size = 20, label_fontface = 'bold')
# 
# #### Maybe a big heatmap with clustering ####
# 
# ## Make sure WT appears only once
# row_wt <- data_fig1f_exp %>% ungroup() %>% filter(WT_Residue == Residue, Position == 2)
# data_big_heatmap <- rbind(data_fig1f_exp %>% ungroup() %>% filter(WT_Residue != Residue), 
#                           row_wt) %>% mutate(ID = ifelse(WT_Residue == Residue, 'WT', ID)) %>%
#   select(ID, exp_level, mean_sel_coeff)
# 
# data_big_heatmap_wide <- data_big_heatmap %>% ungroup() %>%
#   pivot_wider(names_from = exp_level, values_from = mean_sel_coeff)
# 
# # Format the matrix for the heatmap
# needed_rownames <- data_big_heatmap_wide$ID
# data_big_heatmap_wide %<>% select(-ID) %>% as.matrix()
# rownames(data_big_heatmap_wide) <- needed_rownames
# 
# # Prepare an annotation with the cluster IDs
# data_kmeans %<>% ungroup() %>% mutate(Cluster = kmeans_fitness$cluster) %>%
#   rowwise() %>% mutate(ID = str_c(Position, Residue, sep = '')) %>%
#   mutate(ID = ifelse(ID == '2E', 'WT', ID))
# 
# data_kmeans_annot <- data_kmeans %>% select(ID, Cluster)
# 
# ha_all_data <- HeatmapAnnotation(cluster = data_kmeans_annot$Cluster,
#                                  show_annotation_name = T, 
#                                  annotation_name_gp = gpar(fontface = 'bold', fontsize = 24),
#                                  simple_anno_size = unit(0.75, 'cm'),
#                                  # annotation_name_side = 'left',
#                                  col = list(
#                                    cluster = c('1' = "#a6cee3",
#                                                '2' = "#1f78b4",
#                                                '3' = "#b2df8a",
#                                                '4' = "#33a02c",
#                                                '5' = 'black')
#                                    ),
#                                  annotation_legend_param = list(
#                                    cluster = list(
#                                      title_gp = gpar(fontsize = 24, fontface = 'bold'), 
#                                      labels_gp = gpar(fontsize = 22)
#                                    )
#                                  ),
#                                  show_legend = TRUE,
#                                  which = 'row'
#                                  # gp = gpar(col = "black")
#                                  )
# 
# ## Colorbar from -0.8 to +0.4
# seq1 <- seq(-8, 0, length.out = 4) / 10
# seq2 <- seq(0, 4, length.out = 4) / 10
# 
# ## Define variables for heatmap size formatting
# row_title_size = 28
# row_name_size = 2
# column_name_size = 20
# legend_title_size = 28
# legend_label_size = 26
# heatmap_width = 14
# heatmap_height = 95
# legend_height = 7
# 
# p_all_exp_levels <- Heatmap(
#   data_big_heatmap_wide, cluster_columns = F, cluster_rows = T,
#   clustering_distance_rows = 'euclidean', clustering_distance_columns = 'euclidean',
#   col = colorRamp2(
#     breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
#     colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
#   show_column_names = T, row_names_side = 'left',
#   width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
#   border = T,
#   row_title = "Residue",
#   row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
#   row_names_rot = 0, 
#   row_names_centered = T,
#   row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
#   column_names_gp = gpar(fontsize=column_name_size,fontface='bold'),
#   show_heatmap_legend = T,
#   # right_annotation = ha_all_data,
#   left_annotation = ha_all_data,
#   heatmap_legend_param = list(
#     ## Ticks for the scale from -8 to +4
#     at = c(-8, -4,  0, 2, 4) / 10,
#     title = "s", 
#     title_gp = gpar(fontsize = legend_title_size),
#     legend_height = unit(legend_height, "cm"),
#     legend_width = unit(2, "cm"),
#     border='black',
#     lwd=1.7,
#     labels_gp = gpar(fontsize = legend_label_size),
#     title_position = "leftcenter-rot"
#   )
# )
# p_all_exp_levels
# ggsave(grid.grabExpr(draw(p_all_exp_levels)), device = cairo_pdf, width = 10,
#        height = 40, dpi = 300, 
#        filename = 'Figures/2022-05-06_Figures_paper/Heatmap_all_exp_levels.pdf')
# 
#### Comment in the meantime, while we decide what to do with old figures 2B and 2C ####

## Figure 2B ##

# data_fig2b <- data_fig_2_final %>% ungroup() %>%
#   mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff,
#          Arabinose_2 = str_c(Arabinose_2, '% arabinose', sep = ''))
# 
# 
# # Relevel the data set to make sure the boxplots appear in the proper order
# data_fig2b %<>% mutate(Arabinose_2 = factor(Arabinose_2, 
#                                             # levels = c('0.025% arabinose', '0.05% arabinose', '0.2% arabinose'))
#                                             levels = c('0.01% arabinose', '0.025% arabinose', '0.05% arabinose'))
# )
# 
# data_fig2b_exp <- left_join(x = data_fig2b %>% 
#                               separate(col = Arabinose_2, into = c('Arabinose_num', 'tmp'), sep = '% arabinose') %>%
#                               mutate(Arabinose_num = as.numeric(Arabinose_num)),
#                             y = expression_level, 
#                             by = c('Arabinose_num' = 'Arabinose')) %>%
#   mutate(exp_level = as.factor(exp_level))
# 
# ## Make sure WT appears only once
# data_fig2b_exp_nowt <- data_fig2b_exp %>% filter(WT_Residue != Residue)
# data_fig2b_exp_wt <- data_fig2b_exp %>% filter(WT_Residue == Residue, Position == 2)
# 
# data_fig2b_exp_plot <- bind_rows(data_fig2b_exp_nowt, data_fig2b_exp_wt)
# 
# comps <- compare_means(diffNormScore~exp_level, data = data_fig2b_exp_plot, paired = TRUE) %>%
#   mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16', 
#                            sprintf("p = %2.1e", as.numeric(p))
#   ),
#   y_pos = c(0.6, 0.7, 0.8)
#   )
# 
# p_fig2b <- data_fig2b_exp_plot %>% ungroup() %>% 
#   ggplot(aes(x = exp_level, y = diffNormScore,fill = exp_level)) +
#   geom_point(aes(colour = exp_level),alpha = 0.3, 
#              position = position_jitterdodge(jitter.width = 0.25)) +
#   geom_violin(aes(fill = exp_level),alpha = 0.75) +
#   geom_signif(data = as.data.frame(comps), inherit.aes = FALSE, aes(xmin = group1, xmax = group2,
#                                                                     annotations=p.format, y_position = y_pos), 
#               manual = TRUE, textsize = 5) +
#   stat_summary(fun="median", geom="point", size = 5, colour = 'black') +
#   # scale_fill_manual(values = c('#fd8d3c', '#e6002e', '#80001a')) +
#   # scale_colour_manual(values = c('#fd8d3c', '#e6002e', '#80001a')) +
#   scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026')) +
#   scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026')) +
#   ylab('\u0394s') +  # Delta fitness
#   xlab('Expression level') +
#   theme(axis.title = element_text(face = 'bold', size = 20), 
#         axis.text = element_text(size = 18), 
#         legend.text = element_text(size = 20),
#         legend.title = element_text(size = 20),
#         legend.position = 'none', 
#         legend.justification = 0.5, 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank())
# p_fig2b
# 
# ## Figure 2C ##
# 
# data_fig2c <- data_fig_2 %>% ungroup() %>%
#   group_by(Position, Arabinose) %>%
#   summarise(meanSelCoeff = mean(mean_sel_coeff))
# 
# data_fig2c_exp <- left_join(x = data_fig2c, y = expression_level, 
#                             by = c('Arabinose' = 'Arabinose')) %>%
#   mutate(exp_level = as.factor(exp_level))
# 
# p_fig2c <- data_fig2c_exp %>% 
#   filter(Arabinose %in% c(0.01, 0.025, 0.05, 0.2)) %>%
#   mutate(Arabinose = factor(Arabinose, levels = c(0.01, 0.025, 0.05, 0.2))) %>%
#   ggplot(aes(x = exp_level, y = meanSelCoeff, fill = exp_level)) +
#   geom_violin(alpha = 1) +
#   stat_summary(fun="median", geom="point", size = 5) +
#   geom_point() +
#   geom_line(aes(group = Position)) +
#   scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a')) +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'none', 
#         axis.title = element_text(face = 'bold', size = 20), 
#         axis.text = element_text(size = 18), 
#         legend.text = element_text(size = 20),
#         legend.title = element_text(size = 20),
#         legend.justification = 0.5) +
#   xlab('Expression level') + 
#   ylab('Mean s per position')
# p_fig2c

#### Correlation between fitness effects with and without TMP ####

all_data_complete_tmp0_tmp10 <- all_data_complete %>% 
  select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP, mean_sel_coeff) %>%
  filter(Timepoint == 10) %>%
  pivot_wider(names_from = TMP, values_from = mean_sel_coeff, names_prefix = 'TMP_')

temp_nowt <- all_data_complete_tmp0_tmp10 %>% ungroup() %>% group_by(Arabinose) %>%
  filter(Residue != WT_Residue)

wt_row <- all_data_complete_tmp0_tmp10 %>% filter(Residue == WT_Residue, Position == 2)

bind_rows(temp_nowt, wt_row) %>%
  summarise(cor_sel_coeff = cor.test(TMP_0, TMP_10, method = 'spearman')$estimate, 
            pval = cor.test(TMP_0, TMP_10, method = 'spearman')$p.value)


#### Save a dataframe of mutants with significant expression-dependent effects for RF ####

# data_rf_significant <- left_join(
#   x = data_fig1f_exp %>% rowwise() %>%
#     mutate(mut_check_diff = (abs(diffNormScore) > 0.1)) %>%
#     mutate(mut_check_final = and(mut_check_diff, mut_check)) %>%
#     filter(mut_check_final == TRUE, Arabinose == 0.2), # Only significant values, avoid repetition
#     # filter(mut_check == TRUE, Arabinose == 0.2),
#   y = all_data_complete %>% filter(Timepoint == 10, Arabinose == 0.2, TMP == 10) %>% # Avoid repetition
#     select(-mean_sel_coeff, -num_samples, -sd_sel_coeff, -sem_sel_coeff, -Timepoint, -Arabinose, -TMP),
#   by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
# )
# 
# write.table(x = data_rf_significant, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
#             file = 'Data/dataset_diffNorm_ara0.2_ara0.01_index_differences_significant.txt')
# 
# #### Let's try a classification dataset ####
# ## Three classes:
# # - Mutants with no fitness effects (WT-like) (abs(s_opt.) < 0.1, abs(s_weak) < 0.1)
# # - Mutants with fitness effects but insensitive to expression (abs(deltaS) < 0.1)
# # - Mutants that are compensated by higher expression (abs(deltaS) > 0.1)
# data_rf_significant_new <- data_fig1f_exp %>% filter(Arabinose %in% c(0.01, 0.2)) %>%
#   select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, diffNormScore) %>% rowwise() %>%
#   mutate(Arabinose = str_c('ara', Arabinose, sep = '')) %>%
#   pivot_wider(names_from = Arabinose, values_from = mean_sel_coeff)
# 
# 
# data_rf_significant_new %<>% rowwise() %>% 
#   mutate(chk1 = (abs(ara0.2) < 0.1), chk2 = (abs(diffNormScore) < 0.1)) %>%
#   mutate(Class = ifelse(and(chk1, chk2), 'WT-like',
#                         ifelse(chk2, 'Not sensitive to expression', 
#                                ifelse(and(diffNormScore < -0.1, ara0.2 < -0.1), 'Fitter at optimal expression', 
#                                       ifelse(and(diffNormScore < -0.1, ara0.2 > -0.1), 'Fully compensated at optimal expression',
#                                       'Fitter at lower expression')))), 
#     Class_num = ifelse(and(chk1, chk2), 0,
#                        ifelse(chk2, 1, 
#                               ifelse(and(diffNormScore < -0.1, ara0.2 < -0.1), 2, 
#                                      ifelse(and(diffNormScore < -0.1, ara0.2 > -0.1), 3,
#                                             4))))
#     ) %>% select(-chk1, -chk2)
# 
# # Add the rest of the columns
# data_rf_significant_class_final <- left_join(
#   x = data_rf_significant_new,
#   y = all_data_complete %>% filter(Timepoint == 10, Arabinose == 0.2, TMP == 10) %>% # Avoid repetition
#     select(-mean_sel_coeff, -num_samples, -sd_sel_coeff, -sem_sel_coeff, -Timepoint, -Arabinose, -TMP),
#   by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
# )
# 
# write.table(x = data_rf_significant_class_final, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
#             file = 'Data/dataset_diffNorm_ara0.2_ara0.01_index_differences_significant_class.txt')
# 
# ## Figure showing where the classes are:
# p <- data_rf_significant_class_final %>% 
#   ggplot(aes(x = ara0.2, y = diffNormScore, colour = Class)) +
#   geom_point(size = 2) +
#   geom_hline(yintercept = 0.1, linetype = 'dashed') +
#   geom_hline(yintercept = -0.1, linetype = 'dashed') +
#   geom_vline(xintercept = 0.1, linetype = 'dashed') +
#   geom_vline(xintercept = -0.1, linetype = 'dashed') +
#   xlab('Sel. coeff. (optimal expression)') + 
#   labs(y = expression(s[weak] - s[opt])) +
#   theme(legend.position = 'top', 
#         legend.justification = 'center')
# p
# ggsave(p, device = cairo_pdf, width = 14, height = 10, dpi = 300, 
#        filename = 'Figures/2022-04-03_extra_figures/Classes_RF.pdf')

#### Try ANOVA including protein sites ####

# all_data_all_reps_TMP10 <- all_data_all_reps %>% rowwise() %>%
#   filter(TMP == 10, Timepoint == 10, Arabinose != 0.4 
#          # WT_Residue != Residue
#   ) %>%
#   mutate(Genotype = str_c(WT_Residue, Residue, Position),
#          Sample = new_Sample,
#          Secondary_structure = ifelse(Secondary_structure == 'Missing', 'Disordered', Secondary_structure)) %>%
#   select(Sample, Genotype, Timepoint, TMP, Arabinose, sel_coeff, Secondary_structure,Sequencer)
# 
# data_anova <- all_data_all_reps_TMP10 %>% 
#   mutate(Sample = as.factor(Sample),
#          Arabinose = as.factor(Arabinose), 
#          Genotype = as.factor(Genotype), 
#          Secondary_structure = as.factor(Secondary_structure))
# 
# ## Takes a long time to run
# m <- art(sel_coeff ~ Arabinose*Genotype + Arabinose*Genotype*Secondary_structure + Error(Sample), data=data_anova)
# summary(m)
# anova_general <- anova(m)
# 
# anova_general_table <- data.frame(anova_general) %>%
#   mutate(Pr..F. = ifelse(Pr..F. == 0, '<2.2e-16', Pr..F.))
# write.table(x = anova_general_table, append = F, quote = F, sep = '\t', row.names = F, col.names = T,
#             file = 'Figures/2022-04-07_supp_figures/TableSX.ANOVA_general.csv')


#### Figure 2: The fitness landscape at different transcription levels ####

#### Prepare annotation for the bottom of the heatmaps ####

# Load the interface data for the annotation
data_interfaces <- read_delim('Data/Complete_datasets_TMP0_TMP10/interface_data_A.txt', delim = '\t')
data_interfaces %<>% filter(col_interfaces != 'A,B', Core_interface == 'Core')

data_interfaces_wide <- data_interfaces %>% select(Position, col_interfaces, Core_interface) %>%
  pivot_wider(names_from = col_interfaces, values_from = Core_interface)

df_left <- data.frame(Position = c(2:78))

# Load the data on residues that bind to ligand and cofactor
data_ligands <- data.frame(Position = c(32, 67, 68, 69,
                                        32, 35, 36, 50, 64, 65,
                                        66, 67, 68, 69, 70, 72, 73),
                           Ligand = c(rep('DHF', 4),
                                      rep('NADPH', 13)
                                      )
)

data_ligands_wide <- data_ligands %>% pivot_wider(names_from = Ligand, values_from = Ligand)

# Join all the data for the bottom annotation
data_interfaces_temp <- left_join(x = df_left, y = data_interfaces_wide, 
                                  by = c('Position' = 'Position'))

data_interfaces_final <- left_join(x = data_interfaces_temp, y = data_ligands_wide, 
                                   by = c('Position' = 'Position')) %>%
  mutate(`A,C` = ifelse(is.na(`A,C`), 0, 1),
         `A,D` = ifelse(is.na(`A,D`), 0, 1),
         DHF = ifelse(is.na(DHF), 0, 1),
         NADPH = ifelse(is.na(NADPH), 0, 1 ),
         Cat_residues = ifelse(Position %in% c(32, 67, 68, 69), 1, 0), 
         Disordered_region = ifelse(Position %in% c(2:20), 1, 0)
  )

# Load the data on solvent accessibility to add a class of buried residues
sec_struc <- read_delim('/media/axelle/afe8c733-963d-4db8-a2ee-551a0b73c9d7/Angel/PhD_projects/R67_DMS_December2020/Data/Structures/DHFR_2rk1_DSSP_table_bio_rSASA.txt', delim = '\t')

# Add the data on buried residues
data_interfaces_final <- left_join(x = data_interfaces_final, 
                                   y = sec_struc %>% 
                                     mutate(Buried = ifelse(Solvent_accessibility == 'Missing', 0, 
                                                            ifelse(rSASA <= 0.25, 1, 0))) %>%
                                     select(Position, Buried), 
                                   by = c('Position' = 'Position'))

# Save the annotation for the heatmaps
write.table(data_interfaces_final, file = 'Data/Complete_datasets_TMP0_TMP10/data_annotation_2.txt', quote = F, 
            sep = '\t', row.names = F, col.names = T)

#### End prepare annotation ####

#### Check the distribution of RSA values for interface residues ####

# residues_RSA <- left_join(x = all_data_complete %>% 
#                             select(Position, rSASA),
#                           y = data_interfaces_final, 
#                           by = c('Position' = 'Position'))
#   
# residues_RSA_long <- residues_RSA %>% pivot_longer(cols = c("A,C", "A,D", "DHF",
#                                                             "NADPH", "Cat_residues", "Disordered_region"),
#                                                    names_to = "Region", values_to = "val_check") %>%
#   filter(val_check == 1) %>% unique()
# 
# p <- residues_RSA_long %>% ggplot(aes(x = Region, y = rSASA)) +
#   geom_violin() +
#   geom_boxplot(alpha = 0.5, outlier.shape = NA) +
#   geom_jitter(alpha = 1) +
#   geom_hline(yintercept = 0.25, linetype = 'dashed') +
#   xlab('Region') + ylab('Relative solvent accessibility')
# p
# ggsave(file = 'Figures/2022-04-03_extra_figures/check_RSA_interfaces.pdf', plot = p,
#        device = cairo_pdf, width = 10, height = 7, dpi = 300)

#### End check on distributions of RSA values ####

# col_fun = colorRamp2(c(0, 1), c("white", "#595959"))

data_fig_2 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, Entropy)

# Relevel amino acids so that they appear in the correct order
aa_order <- c('*', 'G','A','V', 'L', 'I', 'M', 'C',
              'P', 'W', 'F', 'Y', 'S', 'T', 'N', 'Q', 
              'H', 'K', 'R', 'D', 'E')

data_fig_2$Residue <- factor(data_fig_2$Residue, levels = rev(aa_order))

entropy <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Arabinose == 0.01) %>%
  group_by(Position) %>%
  summarise(Entropy = mean(Entropy))

# Entropy annotation, remove the first value for position 1
ha1 <- HeatmapAnnotation(Entropy = anno_barplot(entropy$Entropy[2:78],
                                                bar_width = 1,
                                                gp = gpar(col = "white", fill = "black"), 
                                                border = FALSE,
                                                gap = unit(1, "points"),
                                                axis=FALSE,
                                                height = unit(3, "cm")
),
show_annotation_name = T,
annotation_name_gp = gpar(fontface = 'bold', fontsize = 24),
annotation_name_side = 'left',
annotation_name_rot = 90,
annotation_name_offset = c(Entropy = '0.15cm')
)

# Need to convert to wide formatted data
ara_0.01 <- data_fig_2 %>%
  filter(Arabinose == 0.01) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.01_final <- as.matrix(ara_0.01 %>% select(-Position))

rownames(ara_0.01_final) <- ara_0.01$Position

# Get a matrix of true/false values for the synonymous codons
ara_0.01_bool <- data_fig_2 %>%
  filter(Arabinose == 0.01) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-Entropy, -WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

ara_0.01_bool_final <- as.matrix(ara_0.01_bool %>% select(-Position))

rownames(ara_0.01_bool_final) <- ara_0.01_bool$Position

# Need to reorder the columns in the matrices
ara_0.01_final <- ara_0.01_final[1:nrow(ara_0.01_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
ara_0.01_bool_final <- ara_0.01_bool_final[1:nrow(ara_0.01_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

## Colorbar from -0.8 to +0.4
seq1 <- seq(-8, 0, length.out = 4) / 10
seq2 <- seq(0, 4, length.out = 4) / 10

## Define variables for heatmap size formatting
row_title_size = 28
# row_name_size = 23
# column_name_size = 24
row_name_size = 20
column_name_size = 20
legend_title_size = 28
legend_label_size = 26
heatmap_width = 39
heatmap_height = 13
legend_height = 7
column_labels = c('', '', '', '5', 
                  '', '', '', '', '10', 
                  '', '', '', '', '15',
                  '', '', '', '', '20', 
                  '', '', '', '', '25', 
                  '', '', '', '', '30', 
                  '', '', '', '', '35', 
                  '', '', '', '', '40', 
                  '', '', '', '', '45', 
                  '', '', '', '', '50', 
                  '', '', '', '', '55', 
                  '', '', '', '', '60', 
                  '', '', '', '', '65', 
                  '', '', '', '', '70', 
                  '', '', '', '', '75', 
                  '', '', ''
)

p2_ara0.01_new <- Heatmap(
  t(ara_0.01_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  # width=unit(30, 'cm'), height = unit(8.5, 'cm'),
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size,fontface='bold'),
  top_annotation = ha1,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[j,i]){
      # grid.text('.', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    ## Ticks for the scale from -8 to +4
    at = c(-8, -4,  0, 2, 4) / 10,
    title = "s", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)

p2_ara0.01_new

### A vertical version of the same heatmap ###
ha1_vert <- HeatmapAnnotation(Entropy = anno_barplot(entropy$Entropy[2:78],
                                                bar_width = 1,
                                                gp = gpar(col = "white", fill = "black"), 
                                                border = FALSE,
                                                gap = unit(1, "points"),
                                                axis=FALSE,
                                                width = unit(3, "cm")
),
show_annotation_name = T,
which = 'row',
annotation_name_gp = gpar(fontface = 'bold', fontsize = 24),
annotation_name_side = 'top',
annotation_name_rot = 0,
annotation_name_offset = c(Entropy = '0.15cm')
)

ha2_vert <- HeatmapAnnotation(# `Interface 1` = data_interfaces_final$`A,C`,
  # `Interface 2` = data_interfaces_final$`A,D`,
  `Dimerization interface` = data_interfaces_final$`A,C`,
  `Tetramerization interface` = data_interfaces_final$`A,D`,
  `DHF binding` = data_interfaces_final$DHF,
  `NADPH binding` = data_interfaces_final$NADPH,
  `Catalytic residues` = data_interfaces_final$Cat_residues,
  `Disordered region` = data_interfaces_final$Disordered_region,
  `Buried residues` = data_interfaces_final$Buried,
  show_annotation_name = T,
  which = 'row',
  annotation_name_gp = gpar(fontface = 'bold', fontsize = 24),
  annotation_name_rot = 90,
  simple_anno_size = unit(0.75, 'cm'),
  annotation_name_side = 'top',
  show_legend = FALSE,
  col = list(# `Interface 1` = colorRamp2(c(0, 1), c("white", "#009E73")), 
    # `Interface 2` = colorRamp2(c(0, 1), c("white", "#0072B2")),
    `Dimerization interface` = colorRamp2(c(0, 1), c("white", "#009E73")), 
    `Tetramerization interface` = colorRamp2(c(0, 1), c("white", "#0072B2")),
    `DHF binding` = colorRamp2(c(0, 1), c("white", "#E69F00")),
    `NADPH binding` = colorRamp2(c(0, 1), c("white", "#D55E00")),
    `Catalytic residues` = colorRamp2(c(0, 1), c("white", "#56B4E9")),
    `Disordered region` = colorRamp2(c(0, 1), c("white", "#CC79A7")), 
    `Buried residues` = colorRamp2(c(0, 1), c("white", "#9933ff"))
  ),
  gp = gpar(col = "black")
)



p2_ara0.01_new_vertical <- Heatmap(
  ara_0.01_final, cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  column_title_side = 'bottom', 
  # width=unit(30, 'cm'), height = unit(8.5, 'cm'),
  width=unit(heatmap_height, 'cm'), height = unit(heatmap_width, 'cm'),
  border = T,
  column_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  row_names_rot = 0, 
  row_names_centered = T,
  column_names_rot = 0,
  column_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size,fontface='bold'),
  # top_annotation = ha1,
  left_annotation = ha2_vert,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[i,j]){
      # grid.text('.', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    ## Ticks for the scale from -8 to +4
    at = c(-8, -4,  0, 2, 4) / 10,
    title = "s", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = column_labels
)

p2_ara0.01_new_vertical


## 0.025% arabinose ##

# Need to convert to wide formatted data
ara_0.025 <- data_fig_2 %>%
  filter(Arabinose == 0.025) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.025_final <- as.matrix(ara_0.025 %>% select(-Position))

rownames(ara_0.025_final) <- ara_0.025$Position

# Get a matrix of true/false values for the synonymous codons
ara_0.025_bool <- data_fig_2 %>%
  filter(Arabinose == 0.025) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-Entropy, -WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

ara_0.025_bool_final <- as.matrix(ara_0.025_bool %>% select(-Position))

rownames(ara_0.025_bool_final) <- ara_0.025_bool$Position

# Need to reorder the columns in the matrices
ara_0.025_final <- ara_0.025_final[1:nrow(ara_0.025_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
ara_0.025_bool_final <- ara_0.025_bool_final[1:nrow(ara_0.025_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p2_ara0.025_new <- Heatmap(
  t(ara_0.025_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size,fontface='bold'),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 2, 4) / 10,
    title = "s", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)
p2_ara0.025_new

#### New vertical version ####

p2_ara0.025_new_vertical <- Heatmap(
  ara_0.025_final, cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_height, 'cm'), height = unit(heatmap_width, 'cm'),
  border = T,
  column_title = "Residue",
  column_title_side = 'bottom',
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  row_names_rot = 0, 
  row_names_centered = T,
  column_names_rot = 0,
  column_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size,fontface='bold'),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[i,j]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 2, 4) / 10,
    title = "s", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = column_labels
)
p2_ara0.025_new_vertical

## Arabinose 0.05%
# Need to convert to wide formatted data
ara_0.05 <- data_fig_2 %>%
  filter(Arabinose == 0.05) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.05_final <- as.matrix(ara_0.05 %>% select(-Position))

rownames(ara_0.05_final) <- ara_0.05$Position

# Get a matrix of true/false values for the synonymous codons
ara_0.05_bool <- data_fig_2 %>%
  filter(Arabinose == 0.05) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-Entropy, -WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

ara_0.05_bool_final <- as.matrix(ara_0.05_bool %>% select(-Position))

rownames(ara_0.05_bool_final) <- ara_0.05_bool$Position

# Need to reorder the columns in the matrices
ara_0.05_final <- ara_0.05_final[1:nrow(ara_0.05_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
ara_0.05_bool_final <- ara_0.05_bool_final[1:nrow(ara_0.05_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p2_ara0.05_new <- Heatmap(
  t(ara_0.05_final), cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  row_names_rot = 90,
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size, fontface = 'bold'),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }
  },
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 2, 4) / 10,
    title = expression(bolditalic("s")),
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)

p2_ara0.05_new

#### New vertical version ####

p2_ara0.05_new_vertical <- Heatmap(
  ara_0.05_final, cluster_columns = F, cluster_rows = F,
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_height, 'cm'), height = unit(heatmap_width, 'cm'),
  border = T,
  column_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 0,
  column_names_rot = 0,
  column_names_centered = T,
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size, fontface = 'bold'),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[i,j]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }
  },
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 2, 4) / 10,
    title = expression(bolditalic("s")),
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = column_labels
)
p2_ara0.05_new_vertical



## Arabinose 0.2%
# Need to convert to wide formatted data
ara_0.2 <- data_fig_2 %>%
  filter(Arabinose == 0.2) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.2_final <- as.matrix(ara_0.2 %>% select(-Position))

rownames(ara_0.2_final) <- ara_0.2$Position

# Get a matrix of true/false values for the synonymous codons
ara_0.2_bool <- data_fig_2 %>%
  filter(Arabinose == 0.2) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-Entropy, -WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

ara_0.2_bool_final <- as.matrix(ara_0.2_bool %>% select(-Position))

rownames(ara_0.2_bool_final) <- ara_0.2_bool$Position

# Need to reorder the columns in the matrices
ara_0.2_final <- ara_0.2_final[1:nrow(ara_0.2_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
ara_0.2_bool_final <- ara_0.2_bool_final[1:nrow(ara_0.2_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p2_ara0.2_new <- Heatmap(
  t(ara_0.2_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=row_name_size,fontface='bold'),
  # bottom_annotation = ha2,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 2, 4) / 10,
    title = "s", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)

p2_ara0.2_new

#### New vertical version ####

p2_ara0.2_new_vertical <- Heatmap(
  ara_0.2_final, cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_height, 'cm'), height = unit(heatmap_width, 'cm'),
  border = T,
  column_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 0, 
  column_names_rot = 0,
  column_names_centered = T,
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=row_name_size,fontface='bold'),
  # bottom_annotation = ha2,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[i,j]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 2, 4) / 10,
    title = "s", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = column_labels
)
p2_ara0.2_new_vertical


## Arabinose 0.4%
# Need to convert to wide formatted data
ara_0.4 <- data_fig_2 %>%
  filter(Arabinose == 0.4) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.4_final <- as.matrix(ara_0.4 %>% select(-Position))

rownames(ara_0.4_final) <- ara_0.4$Position

# Get a matrix of true/false values for the synonymous codons
ara_0.4_bool <- data_fig_2 %>%
  filter(Arabinose == 0.4) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-Entropy, -WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

ara_0.4_bool_final <- as.matrix(ara_0.4_bool %>% select(-Position))

rownames(ara_0.4_bool_final) <- ara_0.4_bool$Position

# Need to reorder the columns in the matrices
ara_0.4_final <- ara_0.4_final[1:nrow(ara_0.4_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
ara_0.4_bool_final <- ara_0.4_bool_final[1:nrow(ara_0.4_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

# Prepare the annotation on interfaces
# Interface A-C will be interface 1
# Interface A-D will be interface 2
ha2 <- HeatmapAnnotation(# `Interface 1` = data_interfaces_final$`A,C`,
  # `Interface 2` = data_interfaces_final$`A,D`,
  `Dimerization interface` = data_interfaces_final$`A,C`,
  `Tetramerization interface` = data_interfaces_final$`A,D`,
  `DHF binding` = data_interfaces_final$DHF,
  `NADPH binding` = data_interfaces_final$NADPH,
  `Catalytic residues` = data_interfaces_final$Cat_residues,
  `Disordered region` = data_interfaces_final$Disordered_region,
  `Buried residues` = data_interfaces_final$Buried,
  show_annotation_name = T,
  annotation_name_gp = gpar(fontface = 'bold', fontsize = 24),
  simple_anno_size = unit(0.75, 'cm'),
  annotation_name_side = 'left',
  show_legend = FALSE,
  col = list(# `Interface 1` = colorRamp2(c(0, 1), c("white", "#009E73")), 
    # `Interface 2` = colorRamp2(c(0, 1), c("white", "#0072B2")),
    `Dimerization interface` = colorRamp2(c(0, 1), c("white", "#009E73")), 
    `Tetramerization interface` = colorRamp2(c(0, 1), c("white", "#0072B2")),
    `DHF binding` = colorRamp2(c(0, 1), c("white", "#E69F00")),
    `NADPH binding` = colorRamp2(c(0, 1), c("white", "#D55E00")),
    `Catalytic residues` = colorRamp2(c(0, 1), c("white", "#56B4E9")),
    `Disordered region` = colorRamp2(c(0, 1), c("white", "#CC79A7")), 
    `Buried residues` = colorRamp2(c(0, 1), c("white", "#9933ff"))
  ),
  gp = gpar(col = "black")
)


p2_ara0.4_new <- Heatmap(
  t(ara_0.4_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size,fontface='bold'),
  bottom_annotation = ha2,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.4_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 2, 4) / 10,
    title = "s", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)
p2_ara0.4_new

#### New vertical version ####

p2_ara0.4_new_vertical <- Heatmap(
  ara_0.4_final, cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(
    breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
    colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
  show_column_names = T, row_names_side = 'left',
  width=unit(heatmap_height, 'cm'), height = unit(heatmap_width, 'cm'),
  border = T,
  column_title = "Residue",
  row_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_gp = gpar(fontsize=row_title_size, fontface = 'bold'),
  column_title_side = 'bottom',
  row_names_rot = 0, 
  row_names_centered = T,
  column_names_rot = 0,
  column_names_centered = T,
  row_names_gp = gpar(fontsize=row_name_size, fontface = 'bold'),
  column_names_gp = gpar(fontsize=column_name_size,fontface='bold'),
  # left_annotation = ha2_vert,
  right_annotation = ha1_vert,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.4_bool_final[i,j]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 2, 4) / 10,
    title = "s", 
    title_gp = gpar(fontsize = legend_title_size),
    legend_height = unit(legend_height, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = legend_label_size),
    title_position = "leftcenter-rot"
  ), 
  row_labels = column_labels
)
p2_ara0.4_new_vertical


# Put the heatmaps together
ht_list = p2_ara0.01_new %v% p2_ara0.025_new%v% p2_ara0.05_new %v% p2_ara0.2_new %v% p2_ara0.4_new
p2_final <- grid.grabExpr(
  draw(ht_list, row_title_gp = gpar(fontsize=20, fontface = 'bold'),
       ht_gap = unit(1, "cm"))
)

#### Save the vertical versions of the heatmaps ####
ht_list_vert = p2_ara0.01_new_vertical + p2_ara0.025_new_vertical + 
  p2_ara0.05_new_vertical + p2_ara0.2_new_vertical + p2_ara0.4_new_vertical
p2_final_vert <- grid.grabExpr(
  draw(ht_list_vert, column_title_gp = gpar(fontsize=20, fontface = 'bold'),
       ht_gap = unit(1, "cm"))
)
ggsave(p2_final_vert, width = 33, height = 21, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-06_Figures_paper/Heatmaps_without_secStruc/2.Fig2_new_all_heatmaps_buried_vertical.pdf')


## Save files for visualization of mean fitness effects per position on the structure
## 0.01% arabinose
ara0.01_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.01, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(mean_s = mean(mean_sel_coeff))
write.table(ara0.01_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/meanFitEff_ara0.01_new.txt')

## 0.025% arabinose
ara0.025_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.025, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(mean_s = mean(mean_sel_coeff))
write.table(ara0.025_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/meanFitEff_ara0.025_new.txt')

## 0.05% arabinose
ara0.05_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.05, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(mean_s = mean(mean_sel_coeff))
write.table(ara0.05_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/meanFitEff_ara0.05_new.txt')

## 0.2% arabinose
ara0.2_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.2, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(mean_s = mean(mean_sel_coeff))
write.table(ara0.2_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/meanFitEff_ara0.2_new.txt')

## 0.4% arabinose
ara0.4_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.4, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(mean_s = mean(mean_sel_coeff))
write.table(ara0.4_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/meanFitEff_ara0.4_new.txt')

## Save tables for visualization of maximum and minimum effects on the DfrB1 structure ####

## 0.01% arabinose
ara0.01_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.01, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(max_s = max(mean_sel_coeff))
write.table(ara0.01_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/max/maxFitEff_ara0.01_new.txt')

ara0.01_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.01, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(min_s = min(mean_sel_coeff))
write.table(ara0.01_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/min/minFitEff_ara0.01_new.txt')

## 0.025% arabinose
ara0.025_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.025, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(max_s = max(mean_sel_coeff))
write.table(ara0.025_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/max/maxFitEff_ara0.025_new.txt')

ara0.025_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.025, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(min_s = min(mean_sel_coeff))
write.table(ara0.025_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/min/minFitEff_ara0.025_new.txt')

## 0.05% arabinose
ara0.05_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.05, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(max_s = max(mean_sel_coeff))
write.table(ara0.05_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/max/maxFitEff_ara0.05_new.txt')

ara0.05_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.05, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(min_s = min(mean_sel_coeff))
write.table(ara0.05_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/min/minFitEff_ara0.05_new.txt')

## 0.2% arabinose
ara0.2_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.2, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(max_s = max(mean_sel_coeff))
write.table(ara0.2_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/max/maxFitEff_ara0.2_new.txt')

ara0.2_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.2, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(min_s = min(mean_sel_coeff))
write.table(ara0.2_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/min/minFitEff_ara0.2_new.txt')

## 0.4% arabinose
ara0.4_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.4, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(max_s = max(mean_sel_coeff))
write.table(ara0.4_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/max/maxFitEff_ara0.4_new.txt')

ara0.4_data <- all_data_complete %>% ungroup() %>%
  filter(Arabinose == 0.4, TMP == 10, Timepoint == 10) %>%
  group_by(Position) %>% summarise(min_s = min(mean_sel_coeff))
write.table(ara0.4_data, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/min/minFitEff_ara0.4_new.txt')

### Load figures of projected s values onto structures
### Mean fitness effects
# ## Load figure for ara0.01
# ara0.01_struc <- ggdraw() + draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/meanFitEff_ara0.01_af2.png')
# ara0.01_struc
# 
# ## Load figure for ara0.025
# ara0.025_struc <- ggdraw() + draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/meanFitEff_ara0.025_af2.png')
# ara0.025_struc
# 
# ## Load figure for ara0.05
# ara0.05_struc <- ggdraw() + draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/meanFitEff_ara0.05_af2.png')
# ara0.05_struc
# 
# # Load figure for ara0.2
# ara0.2_struc <- ggdraw() + draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/meanFitEff_ara0.2_af2.png')
# ara0.2_struc

# ### Maximum fitness effects per position
# ## Load figure for ara0.01
# ara0.01_struc_max <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/max/maxFitEff_ara0.01_af2.png')
# ara0.01_struc_max
# 
# ## Load figure for ara0.025
# ara0.025_struc_max <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/max/maxFitEff_ara0.025_af2.png')
# ara0.025_struc_max
# 
# ## Load figure for ara0.05
# ara0.05_struc_max <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/max/maxFitEff_ara0.05_af2.png')
# ara0.05_struc_max
# 
# ## Load figure for ara0.2
# ara0.2_struc_max <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/max/maxFitEff_ara0.2_af2.png')
# ara0.2_struc_max
# 
# ## Load figure for ara0.4
# ara0.4_struc_max <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/max/maxFitEff_ara0.4_af2.png')
# ara0.4_struc_max
# 
# # Use cowplot to put the figures together
# p_fig2_structures_max <- plot_grid(ara0.01_struc_max + theme(plot.margin = margin(t = 3, r = -10, b = 0, l = 0, 'cm')),
#                                ara0.025_struc_max + theme(plot.margin = margin(t = 1, r = -10, b = 2, l = 0, 'cm')),
#                                ara0.05_struc_max + theme(plot.margin = margin(t = 0, r = -10, b = 3, l = 0, 'cm')),
#                                ara0.2_struc_max + theme(plot.margin = margin(t = -0, r = -10, b = 3, l = 0, 'cm')),
#                                ara0.4_struc_max + theme(plot.margin = margin(t = -2, r = -10, b = 5, l = 0, 'cm')),
#                                nrow = 5)
# 
# ### Repeat for the structures with the minimum score at each position
# ## Load figure for ara0.01
# ara0.01_struc_min <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/min/minFitEff_ara0.01_af2.png')
# ara0.01_struc_min
# 
# ## Load figure for ara0.025
# ara0.025_struc_min <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/min/minFitEff_ara0.025_af2.png')
# ara0.025_struc_min
# 
# ## Load figure for ara0.05
# ara0.05_struc_min <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/min/minFitEff_ara0.05_af2.png')
# ara0.05_struc_min
# 
# ## Load figure for ara0.2
# ara0.2_struc_min <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/min/minFitEff_ara0.2_af2.png')
# ara0.2_struc_min
# 
# ## Load figure for ara0.4
# ara0.4_struc_min <- ggdraw() + draw_image(
#   'Figures/2022-03-24_Figures_paper/chimerax/Figures/min/minFitEff_ara0.4_af2.png')
# ara0.4_struc_min
# 
# # Use cowplot to put the figures together
# p_fig2_structures_min <- plot_grid(ara0.01_struc_min + theme(plot.margin = margin(t = 3, r = 0, b = 0, l = 0.1, 'cm')),
#                                    ara0.025_struc_min + theme(plot.margin = margin(t = 1, r = 0, b = 2, l = 0.1, 'cm')),
#                                    ara0.05_struc_min + theme(plot.margin = margin(t = 0, r = 0, b = 3, l = 0.1, 'cm')),
#                                    ara0.2_struc_min + theme(plot.margin = margin(t = 0, r = 0, b = 3, l = 0.1, 'cm')),
#                                    ara0.4_struc_min + theme(plot.margin = margin(t = -2, r = 0, b = 5, l = 0.1, 'cm')),
#                                    nrow = 5)
# Create labels for each row
# Arabinose 0.01
par(mar = c(0,0,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.3, y = 0.45, paste('Weak\nexpression'), 
     cex = 2.5, col = '#fed976', srt = 90, font = 2)
text_fig_ara0.01 <- ggdraw(recordPlot())
text_fig_ara0.01

# Arabinose 0.025
par(mar = c(0,0,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.3, y = 0.55, paste('Suboptimal\nexpression'), 
     cex = 2.5, col = '#fd8d3c', srt = 90, font = 2)
text_fig_ara0.025 <- ggdraw(recordPlot())
text_fig_ara0.025

# Arabinose 0.05
par(mar = c(0,0,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.3, y = 0.50, paste('Near-optimal\nexpression'), 
     cex = 2.5, col = '#bd0026', srt = 90, font = 2)
text_fig_ara0.05 <- ggdraw(recordPlot())
text_fig_ara0.05

# Arabinose 0.2
par(mar = c(0,0,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.3, y = 0.6, paste('Optimal\nexpression'), 
     cex = 2.5, col = '#80001a', srt = 90, font = 2)
text_fig_ara0.2 <- ggdraw(recordPlot())
text_fig_ara0.2

# Arabinose 0.4
par(mar = c(0,0,0,0))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.3, y = 0.70, paste('Overexpression'), 
     cex = 2.5, col = 'black', srt = 90, font = 2)
text_fig_ara0.4 <- ggdraw(recordPlot())
text_fig_ara0.4


p_fig2_text <- plot_grid(text_fig_ara0.01, text_fig_ara0.025, 
                         text_fig_ara0.05, text_fig_ara0.2,
                         text_fig_ara0.4,
                         nrow = 5, rel_heights = c(1.1, 0.85, 1, 1, 1.1))

## Old figure with the PDB structures on the side added here
# p_fig2_new <- plot_grid(p_fig2_text, 
#                         p2_final, 
#                         p_fig2_structures_max,
#                         p_fig2_structures_min,
#                         ncol = 4,
#                         rel_widths = c(0.8, 1.3, 1, 1))
# ggsave(p_fig2_new, width = 32, height = 27, dpi = 300, device = cairo_pdf, 
# filename = 'Figures/2022-05-06_Figures_paper/Heatmaps_without_secStruc/2.Fig2_new_all.pdf')

## New version (structures will be added with Inkscape)
# p_fig2_new <- plot_grid(p_fig2_text,
#                         p2_final,
#                         ncol = 2,
#                         rel_widths = c(0.5, 1.3))
# ggsave(p_fig2_new, width = 24, height = 27, dpi = 300, device = cairo_pdf,
# filename = 'Figures/2022-05-06_Figures_paper/Heatmaps_without_secStruc/2.Fig2_new_all_new.pdf')

## Try saving the heatmaps and their titles separately
ggsave(p2_final, width = 21, height = 31.5, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-06_Figures_paper/Heatmaps_without_secStruc/2.Fig2_new_all_heatmaps_buried.pdf')
ggsave(p_fig2_text, width = 3, height = 31.5, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-06_Figures_paper/Heatmaps_without_secStruc/2.Fig2_new_all_text.pdf')

#### Figure 3: Heatmap of differences in fitness effects ####

data_fig_3 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose)

# Separate the data for ara 0.2
data_part_1 <- data_fig_3 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_3 %>%
  filter(Arabinose == 0.01)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_3_final <- left_join(x = data_part_1, y = data_part_2, 
                               by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

summary(data_fig_3_final$diffNormScore)

#### Block of code to save the data for the ML analysis and the chimerax figure ####
aa_changes_block <- all_data_complete %>% filter(TMP == 10, Arabinose == 0.01, Timepoint == 10)

## For the old dataframe
# aa_changes <- aa_changes_block[8:ncol(aa_changes_block)]
## For the new dataframe
aa_changes <- aa_changes_block[11:ncol(aa_changes_block)]

data_aa_diffNorm_0.2_0.01 <- bind_cols(data_fig_3_final, 
                                       aa_changes) #%>% select(-Position, -num_samples))

write.table(data_aa_diffNorm_0.2_0.01, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/dataset_diffNorm_ara0.2_ara0.01_index_differences_new.txt')


## Calculate maximum, minimum and mean deltaS at each position
deltaS_max <- data_aa_diffNorm_0.2_0.01 %>% ungroup() %>%
  group_by(Position) %>%
  summarise(max_deltaS = max(diffNormScore))

deltaS_min <- data_aa_diffNorm_0.2_0.01 %>% ungroup() %>%
  group_by(Position) %>%
  summarise(max_deltaS = min(diffNormScore))

deltaS_mean <- data_aa_diffNorm_0.2_0.01 %>% ungroup() %>%
  group_by(Position) %>%
  summarise(max_deltaS = mean(diffNormScore))

# Write the tables
write.table(deltaS_max, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/deltas_weak_opt/max_deltaS_ara0.2_ara0.01.txt')
write.table(deltaS_min, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/deltas_weak_opt/min_deltaS_ara0.2_ara0.01.txt')
write.table(deltaS_mean, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/deltas_weak_opt/mean_deltaS_ara0.2_ara0.01.txt')


#### Need to save the differences in fitness effects for the data without TMP

data_fig_3_notmp <- all_data_complete %>% ungroup() %>%
  filter(TMP == 0, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose)

# Separate the data for ara 0.2
data_part_1_notmp <- data_fig_3_notmp %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2_notmp <- data_fig_3_notmp %>%
  filter(Arabinose == 0.01)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2_notmp) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_3_final_notmp <- left_join(x = data_part_1_notmp, y = data_part_2_notmp, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

summary(data_fig_3_final_notmp$diffNormScore)

aa_changes_block_notmp <- all_data_complete %>% 
  filter(TMP == 0, Arabinose == 0.01, Timepoint == 10) # To avoid repetition

## For the old dataframe
# aa_changes <- aa_changes_block[8:ncol(aa_changes_block)]
## For the new dataframe
aa_changes_notmp <- aa_changes_block_notmp[11:ncol(aa_changes_block_notmp)]

data_aa_diffNorm_0.2_0.01_notmp <- bind_cols(data_fig_3_final_notmp, 
                                       aa_changes_notmp) #%>% select(-Position, -num_samples))

write.table(data_aa_diffNorm_0.2_0.01_notmp, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/dataset_diffNorm_ara0.2_ara0.01_index_differences_new_notmp.txt')


#### End save data without TMP ####

# Get the mean difference per position
mean_diffNorm_0.2_0.01 <- data_aa_diffNorm_0.2_0.01 %>% group_by(Position) %>%
  summarise(Mean_of_differences_fit_eff = mean(diffNormScore))

write.table(mean_diffNorm_0.2_0.01, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/meanDiffFitEff_ara0.2_ara0.01_new.txt')



# Save the standard deviation of selection coefficients for each position
sd_diffNorm_0.2_0.01 <- data_aa_diffNorm_0.2_0.01 %>% group_by(Position) %>%
  summarise(Sd_of_differences_fit_eff = sd(diffNormScore))
write.table(sd_diffNorm_0.2_0.01, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Figures/2022-03-24_Figures_paper/chimerax/Tables/sdDiffFitEff_ara0.2_ara0.01_new.txt')

# Check deltaS for Q67C
data_fig_3_final %>% filter(WT_Residue == 'Q', Residue == 'C', Position == 67)

# Need to convert to wide formatted data
data_fig_3_final_df <- data_fig_3_final %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_3_final <- as.matrix(data_fig_3_final_df %>% select(-Position))

rownames(data_fig_3_final) <- data_fig_3_final_df$Position

# Get a matrix of true/false values for the synonymous codons
fig_3_bool <- data_fig_3 %>%
  filter(Arabinose == 0.2) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_3_bool_final <- as.matrix(fig_3_bool %>% select(-Position))

rownames(fig_3_bool_final) <- fig_3_bool$Position

# Need to reorder the columns in the matrices
data_fig_3_final <- data_fig_3_final[1:nrow(data_fig_3_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_3_bool_final <- fig_3_bool_final[1:nrow(fig_3_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

## Adjust the size of ha1
# Entropy annotation, remove the first value for position 1
ha1 <- HeatmapAnnotation(Entropy = anno_barplot(entropy$Entropy[2:78],
                                                bar_width = 1,
                                                gp = gpar(col = "white", fill = "black"), 
                                                border = FALSE,
                                                gap = unit(1, "points"),
                                                axis=FALSE,
                                                height = unit(2, "cm")
),
show_annotation_name = T,
annotation_name_gp = gpar(fontface = 'bold', fontsize = 16),
annotation_name_side = 'left',
annotation_name_rot = 90,
annotation_name_offset = c(Entropy = '0.15cm')
)


## Adjust the size the annotation on interfaces
# Interface A-C will be interface 1
# Interface A-D will be interface 2
ha2 <- HeatmapAnnotation(# `Interface 1` = data_interfaces_final$`A,C`,
  # `Interface 2` = data_interfaces_final$`A,D`,
  `Dimerization interface` = data_interfaces_final$`A,C`,
  `Tetramerization interface` = data_interfaces_final$`A,D`,
  `DHF binding` = data_interfaces_final$DHF,
  `NADPH binding` = data_interfaces_final$NADPH,
  `Catalytic residues` = data_interfaces_final$Cat_residues,
  `Disordered region` = data_interfaces_final$Disordered_region,
  `Buried residues` = data_interfaces_final$Buried,
  show_annotation_name = T,
  annotation_name_gp = gpar(fontface = 'bold', fontsize = 16),
  # simple_anno_size = unit(0.75, 'cm'),
  annotation_name_side = 'left',
  show_legend = FALSE,
  col = list(# `Interface 1` = colorRamp2(c(0, 1), c("white", "#009E73")), 
    # `Interface 2` = colorRamp2(c(0, 1), c("white", "#0072B2")),
    `Dimerization interface` = colorRamp2(c(0, 1), c("white", "#009E73")), 
    `Tetramerization interface` = colorRamp2(c(0, 1), c("white", "#0072B2")),
    `DHF binding` = colorRamp2(c(0, 1), c("white", "#E69F00")),
    `NADPH binding` = colorRamp2(c(0, 1), c("white", "#D55E00")),
    `Catalytic residues` = colorRamp2(c(0, 1), c("white", "#56B4E9")),
    `Disordered region` = colorRamp2(c(0, 1), c("white", "#CC79A7")),
    `Buried residues` = colorRamp2(c(0, 1), c("white", "#9933ff"))
  ),
  gp = gpar(col = "black")
)


p_fig3a <- Heatmap(
  t(data_fig_3_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-4, 4, length.out = 7) / 10, 
                   colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  row_title = 'Residue',
  row_title_gp = gpar(fontsize=18, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  top_annotation = ha1,
  bottom_annotation = ha2,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_3_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-4, -2, 0, 2, 4) / 10,
    # title = "\u0394s", 
    # title = expression(s[weak] - s[opt]),
    title = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                             bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                             sep = '')),
    title_gp = gpar(fontsize = 24),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 20),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)

p_fig3a

#### Check the distribution of delta(DMS scores) for critical sites ####

data_fig_3 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, Secondary_structure, rSASA)

# Separate the data for ara 0.2
data_part_1 <- data_fig_3 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_3 %>%
  filter(Arabinose == 0.01)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2",
                           "Arabinose_2", "Secondary_structure", "rSASA")

data_fig_3_final_new <- left_join(x = data_part_1, y = data_part_2, 
                                  by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue',
                                         'Residue' = 'Residue', 'rSASA' = 'rSASA',
                                         'Secondary_structure' = 'Secondary_structure')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)


data_fig_3_final_new <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue', 
                                     'rSASA' = 'rSASA', 'Secondary_structure' = 'Secondary_structure')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Merge the dataframes
critical_sites_deltaDMS <- left_join(x = data_fig_3_final_new, y = data_interfaces_final, 
                                     by = c('Position' = 'Position'))

# Pivot to a longer format
critical_sites_deltaDMS %<>% rowwise() %>% mutate(
    `Only NADPH` = ifelse(and(NADPH == 1, 
                              sum(`A,C`, `A,D`, DHF, NADPH, Cat_residues, Buried, Disordered_region) == 1), 1, 0), 
    Unannotated = ifelse(sum(`A,C`, `A,D`, DHF, NADPH, Cat_residues, Buried, Disordered_region) == 0, 1, 0)
  ) %>%
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH', 'Cat_residues', 
                        'Only NADPH', 'Disordered_region', 'Buried', 'Unannotated'), 
                                          names_to = 'Site', values_to = 'Site_check') %>%
  filter(Site_check == 1) %>% 
  mutate(# Site = ifelse(Site == 'A,C', 'Interface 1' , 
         #              ifelse(Site == 'A,D', 'Interface 2', 
         Site = ifelse(Site == 'A,C', 'Dimerization interface' , 
                      ifelse(Site == 'A,D', 'Tetramerization interface', 
                              ifelse(Site == 'Cat_residues', 'Catalytic residues', 
                                     ifelse(Site == 'Disordered_region', 'Disordered region', 
                                            ifelse(Site == 'Buried', 'Buried residues',
                                                   ifelse(Site == 'DHF', 'DHF binding', 
                                                          ifelse(Site == 'NADPH', 'NADPH binding', 
                                                                 Site))))))))
                                                   

residues_site <- critical_sites_deltaDMS %>% select(Position, Site) %>%
  group_by(Site) %>%
  summarise(Position = toString(unique(Position)))

# Color for a gradient background
g <- rasterGrob(brewer.pal(n = 7, name = 'BrBG'), width=unit(1,"npc"), height = unit(1,"npc"), 
                interpolate = TRUE)

# Statistical comparisons
comps <- compare_means(diffNormScore~Site, data = critical_sites_deltaDMS,
                       paired = F) %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
                           ifelse(p > 0.01, str_c('p = ', round(as.numeric(p), 2), sep = ''),
                                  sprintf("p = %2.1e", as.numeric(p)))
  )
  # y_pos = rep(c(0.6, 0.7, 0.8), 8)
  )

p_fig3b <- critical_sites_deltaDMS %>% rowwise() %>%
  mutate(mut_id = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  filter(Site != 'Only NADPH') %>%
  mutate(Site = factor(Site, levels = c('Disordered region', 'Catalytic residues',
                                        'DHF binding', 'NADPH binding',
                                        # 'Interface 1', 'Interface 2',
                                        'Buried residues',
                                        'Tetramerization interface', 'Dimerization interface', 
                                        'Unannotated')
                       )
         ) %>%
  ggplot(aes(x = Site, y = diffNormScore, colour = Site, fill = Site)) +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  scale_colour_manual(values = c('#CC79A7', '#56B4E9', '#E69F00', '#D55E00',
                                 '#9933ff',
                                 '#0072B2', '#009E73',
                                 '#000000')) +
  scale_fill_manual(values = c('#CC79A7', '#56B4E9', '#E69F00', '#D55E00',
                               '#9933ff',
                               '#0072B2', '#009E73',
                               '#000000')) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.8) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(legend.position = 'none', 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(face = 'bold', size = 20),
        axis.title.y = element_text(face = 'bold', size = 28),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 18), 
        plot.margin = margin(t = 0, l = 5, b = 0, r = 4, 'cm')) +
  xlab('') + # ylab('\u0394s') +
  # labs(y = expression(s[weak] - s[opt])) +
  labs(y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                            bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                            sep = ''))) +
  ylim(-0.4, 0.4) +  
  annotate('text', x = 2, y = 0.3,
           # label = 's[weak] > s[opt]',
           label = expression(italic(s[weak] > s[opt])),
           parse = T, size = 10) +
  annotate('text', x = 2, y = -0.3,
           # label = 's[weak] < s[opt]',
           label = expression(italic(s[weak] > s[opt])),
           parse = T, size = 10)
p_fig3b

#### Figure 3C: Load the PDB structures with deltaS values ####

p_fig3c_min <- ggdraw() + 
    draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/deltas_weak_opt/2rk1_min_deltaS_ara0.01_ara0.2.png')
p_fig3c_mean <- ggdraw() + 
  draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/deltas_weak_opt/2rk1_mean_deltaS_ara0.01_ara0.2.png')
p_fig3c_max <- ggdraw() + 
  draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/deltas_weak_opt/2rk1_max_deltaS_ara0.01_ara0.2.png')


#### Figure 3D: deltaS vs s (optimal), protein sites ####

## Load data from validations to label them in figure 3C
# Read data for selected mutants
file.od <- 'Validations_July2021/AFFC_growthrate_script/14_07_21_bact.xlsx'
plate.ind <- 'Validations_July2021/AFFC_growthrate_script/Mutant_id_validation_DMS_growth_curves_IGA_14_07_21_no_empty.csv'

# function to process plates ----------------------------------------------
read.my.gc <- function(file, plate.index){
  pl <- read.xlsx(file,sheetIndex = 1, rowIndex = 4:99, stringsAsfators = FALSE,
                  header = F)
  # ind <- read.csv(plate.index, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  ind <- read.csv2(plate.index, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  # time <- seq(0,0.25*(ncol(pl)-2), 0.25)
  time <- seq(0,0.30*(ncol(pl)-2), 0.30)
  
  colnames(pl)[1] <- "Well"
  colnames(pl)[2:ncol(pl)] <- time
  pl %<>% select(1:(ncol(pl)-2)) 
  
  data.pl <- gather(pl, key = "time", value = "OD",2:ncol(pl), convert = F)
  data.pl$time <- as.numeric(data.pl$time)
  data.pl$OD <- as.numeric(data.pl$OD)
  data.pl %<>% left_join(ind, by = "Well")
  
  d <-  select(data.pl, 1:3)
  d %<>% spread(key = "Well", value = "OD", convert = F) 
  colnames(d)[1] <- "time"
  
  ## Use Growthcurver as follows
  ## I stop at t = 14.4 because that is the time it takes the last WT culture to get to its max OD
  gc_out <- SummarizeGrowthByPlate(d, t_trim =14.4)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
}

data.od1 <- read.my.gc(file.od, plate.ind)

## Figure 3D ##
all_data_complete_new <- all_data_complete %>% filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose) %>%
  rowwise() %>%
  mutate(ID = str_c(Position, WT_Residue, Residue, '_', Arabinose, sep = '')) %>%
  mutate(ID = ifelse(ID == '2EE_0.01', 'WT_0.01', ID)) %>%
  mutate(ID = ifelse(ID == '2EE_0.2', 'WT_0.2', ID))

data.od1.new <- data.od1 %>% 
  rowwise() %>%
  mutate(ID = str_c(Mutant, Arabinose, sep = '_')) %>%
  group_by(ID, TMP) %>%
  summarise(k = mean(k),
            n0 = mean(n0),
            r = mean(r),
            t_mid = mean(t_mid),
            t_gen = mean(t_gen),
            auc_l = mean(auc_l),
            auc_e = mean(auc_e),
            sigma = mean(sigma))

joined_sets <- inner_join(x = all_data_complete_new, 
                          y = data.od1.new %>% filter(TMP == 10),
                          by = c('ID' = 'ID')) %>%
  mutate(Arabinose = as.factor(Arabinose))

joined_sets %<>% mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                                           ifelse(Arabinose == 0.025, 'Suboptimal', 
                                                  ifelse(Arabinose == 0.05, 'Near-optimal', 
                                                         ifelse(Arabinose == 0.2, 'Optimal', NA)))))

# Mutants to annotate
list_mut <- joined_sets %>% select(Position, WT_Residue, Residue) %>% ungroup() %>%
  unique() %>% mutate(Genotype = str_c(WT_Residue, Position, Residue, sep = ''))

# Separate the data for ara 0.2
data_part_1 <- data_fig_3 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_3 %>%
  filter(Arabinose == 0.01)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2", 
                           "Secondary_structure", "rSASA")


# Join
data_fig_3_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue', 
                                     'Secondary_structure' = 'Secondary_structure', 'rSASA' = 'rSASA')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

mutants_highlight <- data_fig_3_final %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop',
                                   # 'Missense')), 
                                   'Amino acid substitution')),
         Genotype = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  filter(Genotype %in% list_mut$Genotype)

# Rename WT
list_mut %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))
mutants_highlight %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))

# Get data to represent the stop and WT codons
summary_wt_stop <- data_fig_3_final %>% rowwise() %>% ungroup() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop',
                                   # 'Missense'))) %>% 
                                   'Amino acid substitution'))) %>%
  filter(mut_check != 'Amino acid substitution', Position >= 30, Position <= 70) %>%
  group_by(mut_check) %>%
  summarise(mean_s = mean(mean_sel_coeff), sem_s = sd(mean_sel_coeff) / sqrt(n()), 
            mean_ds = mean(diffNormScore), sem_ds = sd(diffNormScore) / sqrt(n()))

data_fig3d <- data_fig_3_final %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop',
                                   'Amino acid substitution'))) %>%
  filter(mut_check == 'Amino acid substitution') %>% 
  mutate(sem_s = 0, sem_ds = 0) %>% # These are single points
  select(mut_check, mean_sel_coeff, sem_s, diffNormScore, sem_ds)

colnames(data_fig3d) <- c('mut_check', 'mean_s','sem_s', 'mean_ds', 'sem_ds')

data_fig_3d <- rbind(data_fig3d, summary_wt_stop)

## Rename the column for mutants to highlight
mutants_highlight %<>% mutate(mean_s = mean_sel_coeff, mean_ds = diffNormScore)

# Draw the figure
p_fig3d <- data_fig_3d %>% rowwise() %>%
  mutate(mut_check = factor(mut_check, 
                            levels = c('Amino acid substitution', 'Stop', 'WT'))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(alpha = mut_check, size = mut_check, 
                 colour = mut_check, shape = mut_check)) +
  geom_errorbar(aes(ymax = mean_ds + sem_ds, ymin = mean_ds - sem_ds, 
                    colour = mut_check), show.legend = F) +
  geom_errorbarh(aes(xmax = mean_s + sem_s, xmin = mean_s - sem_s, 
                     colour = mut_check), show.legend = F) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 7, fontface = 'bold') +
  # xlab('s (optimal expression)') + 
  scale_size_manual(values = c(3, 5, 5)) +
  guides(alpha = 'none') +
  # scale_colour_manual(values = c('grey', 'red', 'blue')) +
  scale_colour_manual(values = c('grey', 'black', 'blue')) +
  scale_shape_manual(values = c(16, 8, 16)) +
  scale_alpha_manual(values = c(0.5, 1, 1)) +
  labs(colour = 'Mutation type', shape = 'Mutation type', size = 'Mutation type',
       # y = expression(s[weak] - s[opt])) +
       y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                                 bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                                 sep = '')), 
       x =  expression(paste(bolditalic('s '), bold('(optimal expression)'), 
                             sep = ''))
      ) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 20), 
    axis.title.y = element_text(face = 'bold', size = 24),
    axis.text = element_text(size = 18), 
    legend.position = 'top',
    legend.justification = 'center', 
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           # label = 'Higher fitness at low expression s[weak] > s[opt]', parse = T) +
           # label = 's[weak] > s[opt]',
           label = expression(italic(s[weak] > s[opt])),
           parse = T, size = 10) +
  annotate('text', x = -0.7, y = -0.5, 
           # label = 'Higher fitness at optimal expression sweak < sopt', parse = T)
           # label = 's[weak] < s[opt]',
           label = expression(italic(s[weak] > s[opt])),
           parse = T, size = 10)

p_fig3d

## Figure 3E
data_regions <- read_delim('Data/Complete_datasets_TMP0_TMP10/data_annotation_2.txt', delim = '\t') %>%
  rowwise() %>%
  mutate(Unannotated = ifelse(sum(DHF, NADPH, `A,C`, `A,D`, Disordered_region, Cat_residues) == 0, 1, 0)) %>%
  pivot_longer(cols = c('DHF', 'NADPH', 'A,C', 'A,D', 'Disordered_region', 'Cat_residues', 'Buried', 'Unannotated'), 
               names_to = 'Region', values_to = 'Value') %>%
  filter(Value == 1)

data_regions_new <- left_join(x = data_fig_3_final, 
                              y = data_regions %>% select(Position, Region),
                              by = c('Position' = 'Position')) %>%
  unique()

# Add regions to highlighted mutants
mutants_highlight_regions <- left_join(x = mutants_highlight,
                                       y = data_regions_new %>% select(Position, WT_Residue, Residue, 
                                                                       Region), 
                                       by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                                              'Residue' = 'Residue'))  %>%
  mutate(# Region = ifelse(Region == 'A,C', 'Interface 1', 
                         # ifelse(Region == 'A,D', 'Interface 2', 
         Region = ifelse(Region == 'A,C', 'Dimerization interface', 
                        ifelse(Region == 'A,D', 'Tetramerization interface', 
                                ifelse(Region == 'Disordered_region', 'Disordered region',
                                       ifelse(Region == 'Cat_residues', 'Catalytic residues',
                                       # ifelse(Region == 'Cat_residues', 'DHF binding / catalytic',
                                              ifelse(Region == 'Buried', 'Buried residues', 
                                                     Region))))))

# p_fig3d <- data_regions_new %>% rowwise() %>%
#   mutate(# Region = ifelse(Region == 'A,C', 'Interface 1', 
#     # ifelse(Region == 'A,D', 'Interface 2', 
#     Region = ifelse(Region == 'A,C', 'Dimerization interface', 
#                     ifelse(Region == 'A,D', 'Tetramerization interface',
#                                 ifelse(Region == 'Disordered_region', 'Disordered region',
#                                        ifelse(Region == 'Cat_residues', 'Catalytic residues',
#                                               Region))))) %>% 
#   # If I don't want to include unannotated
#   filter(Region != 'Unannotated', Region != 'Buried') %>%
#   mutate(Region_stop = ifelse(Residue == '*', 'Stop', Region)) %>%
#   mutate(Region = factor(Region, levels = c('Disordered region', 'Catalytic residues', 'DHF', 'NADPH',
#                                             # 'Interface 1', 'Interface 2',
#                                             'Tetramerization interface', 'Dimerization interface'
#   ))) %>%
#   mutate(Region_stop = factor(Region_stop, levels = c('Disordered region', 'Catalytic residues',
#                                                     'DHF', 'NADPH',
#                                                     'Tetramerization interface', 'Dimerization interface', 
#                                                     'Stop'))) %>%
#   arrange(Region_stop) %>%
#   ggplot(aes(x = mean_sel_coeff, y = diffNormScore)) +
#   geom_point(# size = 3, 
#              aes(colour = Region_stop, size = Region_stop, 
#                  shape = Region_stop), alpha = 0.7) +
#   # geom_hex() +
#   geom_label_repel(data = mutants_highlight_regions %>% 
#                      # When I don't want to include the unannotated
#                      filter(Region != 'Unannotated', Genotype != 'WT') %>%
#                      mutate(Region = factor(Region, levels = c('Disordered region', 'Catalytic residues',
#                                                                'DHF', 'NADPH',
#                                                                # 'Interface 1', 'Interface 2',
#                                                                'Tetramerization interface', 'Dimerization interface' 
#                                                                 
#                      ))),
#                    aes(label = Genotype, group = Region), 
#                    box.padding = 0.4, size = 5, fontface = 'bold') +
#   facet_wrap(~Region, nrow = 3, scales = 'free') +
#   xlab('s (optimal expression)') +
#   labs(y = expression(s[weak] - s[opt])) +
#   scale_colour_manual(values = c('#CC79A7', '#56B4E9', '#E69F00', '#D55E00',
#                                  '#0072B2', '#009E73',
#                                  'black'
#                                  # 'red'# ,
#                                  # '#000000'
#   )) +
#   # scale_alpha_manual(values = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 1)) +
#   scale_size_manual(values = c(3, 3, 3, 3, 3, 3, 4)) +
#   scale_shape_manual(values = c(16, 16, 16, 16, 16, 16, 8)) +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         axis.title.x = element_text(face = 'bold', size = 20), 
#         axis.title.y = element_text(face = 'bold', size = 28),
#         axis.text = element_text(size = 18), 
#         legend.position = 'none',
#         legend.justification = 'center', 
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 18), 
#         strip.text = element_text(size = 18, face = 'bold'), 
#         strip.background = element_rect(fill = 'white'), 
#         axis.line=element_line()
#   ) + 
#   labs(colour = 'Site') +
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4)
# p_fig3d

#### New version that merges catalytic residues with DHF binding and adds buried residues ####
data_regions <- read_delim('Data/Complete_datasets_TMP0_TMP10/data_annotation_2.txt', delim = '\t') %>%
  rowwise() %>%
  mutate(Unannotated = ifelse(sum(DHF, NADPH, `A,C`, `A,D`, Disordered_region, Cat_residues) == 0, 1, 0)) %>%
  pivot_longer(cols = c('DHF', 'NADPH', 'A,C', 'A,D', 'Disordered_region', 'Cat_residues', 'Buried', 'Unannotated'), 
               names_to = 'Region', values_to = 'Value') %>%
  filter(Value == 1)

data_regions_new <- left_join(x = data_fig_3_final, 
                              y = data_regions %>% select(Position, Region),
                              by = c('Position' = 'Position')) %>%
  unique()

# Add regions to highlighted mutants
mutants_highlight_regions <- left_join(x = mutants_highlight,
                                       y = data_regions_new %>% select(Position, WT_Residue, Residue, 
                                                                       Region), 
                                       by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                                              'Residue' = 'Residue'))  %>%
  filter(Region != 'DHF') %>%
  mutate(# Region = ifelse(Region == 'A,C', 'Interface 1', 
    # ifelse(Region == 'A,D', 'Interface 2', 
    Region = ifelse(Region == 'A,C', 'Dimerization interface', 
                    ifelse(Region == 'A,D', 'Tetramerization interface', 
                           ifelse(Region == 'Disordered_region', 'Disordered region',
                                  # ifelse(Region == 'Cat_residues', 'Catalytic residues',
                                  ifelse(Region == 'Cat_residues', 'DHF binding / catalytic',
                                         ifelse(Region == 'Buried', 'Buried residues', 
                                                ifelse(Region == 'NADPH', 'NADPH binding',
                                                       Region)))))))

data_regions_new %<>% rowwise() %>%
  mutate(# Region = ifelse(Region == 'A,C', 'Interface 1', 
    # ifelse(Region == 'A,D', 'Interface 2', 
    Region = ifelse(Region == 'A,C', 'Dimerization interface', 
                    ifelse(Region == 'A,D', 'Tetramerization interface',
                           ifelse(Region == 'Disordered_region', 'Disordered region',
                                  ifelse(Region == 'Cat_residues', 'DHF binding / catalytic',
                                         ifelse(Region == 'NADPH', 'NADPH binding',
                                                ifelse(Region == 'Buried', 'Buried residues',
                                                       Region))))))) %>% 
  # If I don't want to include unannotated
  filter(Region != 'Unannotated', Region != 'DHF') %>%
  mutate(Region_stop = ifelse(Residue == '*', 'Stop', Region)) %>%
  mutate(Region = factor(Region, levels = c('Disordered region',#  'Catalytic residues',
                                            'DHF binding / catalytic', 'NADPH binding',
                                            'Buried residues',
                                            # 'Interface 1', 'Interface 2',
                                            'Tetramerization interface', 'Dimerization interface'
  ))) %>%
  mutate(Region_stop = factor(Region_stop, levels = c('Disordered region', # 'Catalytic residues',
                                                      'DHF binding / catalytic', 'NADPH binding',
                                                      'Buried residues',
                                                      'Tetramerization interface', 'Dimerization interface', 
                                                      'Stop'))) %>%
  arrange(Region_stop)

mutants_highlight_regions %<>% 
  # When I don't want to include the unannotated
  filter(Region != 'Unannotated', Genotype != 'WT') %>%
  mutate(Region = factor(Region, levels = c('Disordered region',#  'Catalytic residues',
                                            'DHF binding / catalytic', 'NADPH binding',
                                            'Buried residues',
                                            # 'Interface 1', 'Interface 2',
                                            'Tetramerization interface', 'Dimerization interface' 
                                            
  )))


p_fig3e <- data_regions_new %>%
  ggplot(aes(x = mean_sel_coeff, y = diffNormScore)) +
  geom_point(# size = 3, 
    aes(colour = Region_stop, size = Region_stop, 
        shape = Region_stop), alpha = 0.7) +
  # geom_hex() +
  geom_label_repel(data = mutants_highlight_regions %>% mutate(Region_stop = Region),
                   aes(label = Genotype, group = Region_stop), 
                   box.padding = 0.4, size = 5, fontface = 'bold') +
  facet_wrap(~Region, nrow = 3, scales = 'free') +
  # xlab('s (optimal expression)') +
  # labs(y = expression(s[weak] - s[opt])) +
  labs(
  y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                       bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                       sep = '')), 
  x =  expression(paste(bolditalic('s '), bold('(optimal expression)'), 
                      sep = ''))
  ) +
  scale_colour_manual(values = c('#CC79A7', #'#56B4E9',
                                 '#E69F00', '#D55E00',
                                 '#9933ff',
                                 '#0072B2', '#009E73',
                                 'black'
                                 # 'red'# ,
                                 # '#000000'
  )) +
  # scale_alpha_manual(values = c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 1)) +
  scale_size_manual(values = c(3, 3, 3, 3, 3, 3, 4)) +
  scale_shape_manual(values = c(16, 16, 16, 16, 16, 16, 8)) +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(face = 'bold', size = 20), 
        axis.title.y = element_text(face = 'bold', size = 28),
        axis.text = element_text(size = 18), 
        legend.position = 'none',
        legend.justification = 'center', 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        strip.text = element_text(size = 18, face = 'bold'), 
        strip.background = element_rect(fill = 'white'), 
        axis.line=element_line()
  ) + 
  labs(colour = 'Site') +
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4)
p_fig3e


# Put the figures together
p_fig3c <- plot_grid(p_fig3c_min, p_fig3c_mean, p_fig3c_max, ncol = 3, 
                     labels = c('C', '', ''), label_size = 20, label_fontface = 'bold')

p_fig3de <- plot_grid(p_fig3d, p_fig3e, ncol = 2, labels = c('D', 'E'), label_fontface = 'bold', 
                    label_size = 20)



#### Join all the panels of figure 3 ####
p_fig3_new <- plot_grid(grid.grabExpr(draw(p_fig3a)), 
                        p_fig3b,
                        p_fig3c,
                        p_fig3de,
                        nrow = 4, labels = c('A', 'B', ''), 
                        rel_heights = c(1, 1, 1, 1), label_size = 20, 
                        label_fontface = 'bold')

ggsave(p_fig3_new, device = cairo_pdf, width = 17, height = 30, dpi = 300, 
       filename = 'Figures/2022-05-06_Figures_paper/Heatmaps_without_secStruc/Fig3_buried_2022-06-03.pdf')
       
#### Figure 4: Slightly destabilizing mutants are masked by higher expression ####
#### Figure 4A ####

all_data_complete_new <- all_data_complete %>% filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D) %>%
  rowwise() %>%
  mutate(ID = str_c(Position, WT_Residue, Residue, '_', Arabinose, sep = '')) %>%
  mutate(ID = ifelse(ID == '2EE_0.01', 'WT_0.01', ID)) %>%
  mutate(ID = ifelse(ID == '2EE_0.2', 'WT_0.2', ID))

data.od1.new <- data.od1 %>% 
  rowwise() %>%
  mutate(ID = str_c(Mutant, Arabinose, sep = '_')) %>%
  group_by(ID, TMP) %>%
  summarise(k = mean(k),
            n0 = mean(n0),
            r = mean(r),
            t_mid = mean(t_mid),
            t_gen = mean(t_gen),
            auc_l = mean(auc_l),
            auc_e = mean(auc_e),
            sigma = mean(sigma))

joined_sets <- inner_join(x = all_data_complete_new, 
                          y = data.od1.new %>% filter(TMP == 10),
                          by = c('ID' = 'ID')) %>%
  mutate(Arabinose = as.factor(Arabinose))

joined_sets %<>% mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                                           ifelse(Arabinose == 0.025, 'Suboptimal', 
                                                  ifelse(Arabinose == 0.05, 'Near-optimal', 
                                                         ifelse(Arabinose == 0.2, 'Optimal', NA)))))

# Mutants to annotate
list_mut <- joined_sets %>% select(Position, WT_Residue, Residue) %>% ungroup() %>%
  unique() %>% mutate(Genotype = str_c(WT_Residue, Position, Residue, sep = ''))

data_fig_4a <- all_data_complete_new %>% ungroup() %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D)

# Separate the data for ara 0.2
data_part_1 <- data_fig_4a %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_4a %>%
  filter(Arabinose == 0.01)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2", 
                           "Mean_ddG_stab_HET", "Mean_ddG_int_HM_A_C", "Mean_ddG_int_HM_A_D")

# Join
data_fig_4a_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue',
                                     'Residue' = 'Residue', 'Mean_ddG_stab_HET' = 'Mean_ddG_stab_HET',
                                     'Mean_ddG_int_HM_A_C' = 'Mean_ddG_int_HM_A_C', 
                                     'Mean_ddG_int_HM_A_D' = 'Mean_ddG_int_HM_A_D')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

mutants_highlight <- data_fig_4a_final %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense')), 
         Genotype = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  filter(Genotype %in% list_mut$Genotype)

# Rename WT
list_mut %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))
mutants_highlight %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))

# Get data to represent the stop and WT codons
summary_wt_stop <- data_fig_4a_final %>% rowwise() %>% ungroup() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense'))) %>% 
  filter(mut_check != 'Missense') %>%
  group_by(mut_check) %>%
  summarise(mean_s = mean(mean_sel_coeff), sem_s = sd(mean_sel_coeff) / sqrt(n()), 
            mean_ds = mean(diffNormScore), sem_ds = sd(diffNormScore) / sqrt(n()),
            ddG_stab = mean(Mean_ddG_stab_HET), 
            ddG_dim_int = mean(Mean_ddG_int_HM_A_C), 
            ddG_tet_int = mean(Mean_ddG_int_HM_A_D))

data_fig4a <- data_fig_4a_final %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
  filter(mut_check %in% c('Missense', 'WT')) %>% 
  mutate(sem_s = 0, sem_ds = 0) %>% # These are single points
  select(mut_check, mean_sel_coeff, sem_s, diffNormScore, sem_ds, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D) %>%
  ungroup()

colnames(data_fig4a) <- c('mut_check', 'mean_s','sem_s', 'mean_ds', 'sem_ds', 
                          'ddG_stab', 'ddG_dim_int', 'ddG_tet_int')

data_fig_4a <- rbind(data_fig4a, summary_wt_stop)

## Rename the column for mutants to highlight
mut_ht_wt <- mutants_highlight %>% filter(Genotype == 'WT') %>%
  rowwise() %>%
  # Set mean_s and mean_ds for the WT to zero to annotate the corresponding point
  mutate(mean_s = ifelse(Genotype == 'WT', 0, mean_s), 
         mean_ds = ifelse(Genotype == 'WT', 0, mean_ds))
mutants_highlight %<>% mutate(mean_s = mean_sel_coeff, mean_ds = diffNormScore) %>% 
  filter(!(is.na(Mean_ddG_stab_HET))) 

mutants_highlight <- bind_rows(mut_ht_wt, mutants_highlight)

# Add bins of ddG values
data_fig_4a %<>%
  mutate(bins_stab = cut(ddG_stab,
                         breaks = quantile(ddG_stab, probs = c(0, 0.25, 0.5, 0.75, 1), 
                                           include.lowest = TRUE, na.rm = TRUE)), 
         bins_dim_int = cut(ddG_dim_int,
                            breaks = quantile(ddG_dim_int, probs = c(0, 0.25, 0.5, 0.75, 1), 
                                              include.lowest = TRUE, na.rm = TRUE)),
         bins_tet_int = cut(ddG_tet_int,
                            breaks = quantile(ddG_tet_int, probs = c(0, 0.25, 0.5, 0.75, 1), 
                                              include.lowest = TRUE, na.rm = TRUE)),
  )

table(data_fig_4a$bins_stab)
table(data_fig_4a$bins_dim_int)
table(data_fig_4a$bins_tet_int)

# data_fig_4a %<>%
#   mutate(bins_stab = cut(ddG_stab,
#                          breaks = c(-4, -0.5, 0, 2, 42)), 
#          bins_dim_int = cut(ddG_dim_int,
#                             breaks = c(-4, -0.5, 0, 2, 42)),
#          bins_tet_int = cut(ddG_tet_int,
#                             breaks = c(-4, -0.5, 0, 2, 42)),
#   )


# # Draw the figure for ddG stability
# p_fig4a_stab <- data_fig_4a %>% rowwise() %>%
#   filter(!(is.na(ddG_stab))) %>%
#   mutate(ddG_stab = ifelse(ddG_stab > 10, 10, ddG_stab)) %>%
#   ggplot(aes(x = mean_s, y = mean_ds)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
#   geom_point(aes(# alpha = mut_check, 
#                  size = mut_check, 
#                  colour = ddG_stab)) +
#   # geom_errorbar(aes(ymax = mean_ds + sem_ds, ymin = mean_ds - sem_ds, 
#   #                   colour = mut_check), show.legend = F) +
#   # geom_errorbarh(aes(xmax = mean_s + sem_s, xmin = mean_s - sem_s, 
#   #                    colour = mut_check), show.legend = F) +
#   geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
#                    box.padding = 0.4, size = 4, fontface = 'bold') +
#   xlab('Sel. coeff. (optimal expression)') + 
#   scale_size_manual(values = c(3, 5, 5)) +
#   guides(size = 'none', alpha = 'none') +
#   # scale_colour_manual(values = c('grey', 'red', 'blue')) +
#   # scale_colour_gradient2(mid = 'orange', high = 'red', midpoint = 0) +
#   # scale_colour_gradient2() +
#   scale_colour_viridis(option = 'A') +
#   # scale_alpha_manual() +
#   labs(colour = '\u0394\u0394G subunit stability [kcal/mol]', y = expression(s[weak] - s[opt])) +
#   theme(
#     axis.title.x = element_text(face = 'bold', size = 20), 
#     axis.title.y = element_text(face = 'bold', size = 28),
#     axis.text = element_text(size = 18), 
#     legend.position = 'top',
#     # panel.background = element_rect(fill = 'grey'),
#     legend.justification = 'center', 
#     legend.key.width = unit(1.15, 'cm'),
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 16)
#   ) + 
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
#   annotate('text', x = -0.7, y = 0.25, 
#            # label = 'Higher fitness at low expression s[weak] > s[opt]', parse = T) +
#            label = 's[weak] > s[opt]', parse = T, size = 8) +
#   annotate('text', x = -0.7, y = -0.5, 
#            # label = 'Higher fitness at optimal expression sweak < sopt', parse = T)
#            label = 's[weak] < s[opt]', parse = T, size = 8)
# 
# p_fig4a_stab
# # ggsave(p_fig3c, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
# #        filename = 'Figures/2022-05-06_Figures_paper/Fig4_test_ddG_stab.pdf')
# 
# ## Figure for ddG dimerization interface
# p_fig4a_dim_int <- data_fig_4a %>% rowwise() %>%
#   filter(!(is.na(ddG_dim_int))) %>%
#   mutate(ddG_dim_int = ifelse(ddG_dim_int > 10, 10, ddG_dim_int)) %>%
#   ggplot(aes(x = mean_s, y = mean_ds)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
#   geom_point(aes(# alpha = mut_check, 
#     size = mut_check, 
#     colour = ddG_dim_int)) +
#   # geom_errorbar(aes(ymax = mean_ds + sem_ds, ymin = mean_ds - sem_ds, 
#   #                   colour = mut_check), show.legend = F) +
#   # geom_errorbarh(aes(xmax = mean_s + sem_s, xmin = mean_s - sem_s, 
#   #                    colour = mut_check), show.legend = F) +
#   geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
#                    box.padding = 0.4, size = 4, fontface = 'bold') +
#   xlab('Sel. coeff. (optimal expression)') + 
#   scale_size_manual(values = c(3, 5, 5)) +
#   guides(size = 'none', alpha = 'none') +
#   # scale_colour_manual(values = c('grey', 'red', 'blue')) +
#   # scale_colour_gradient2(mid = 'orange', high = 'red', midpoint = 0) +
#   # scale_colour_gradient2() +
#   # scale_alpha_manual() +
#   scale_colour_viridis(option = 'A') +
#   labs(colour = '\u0394\u0394G dim. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
#   theme(
#     axis.title.x = element_text(face = 'bold', size = 20), 
#     axis.title.y = element_text(face = 'bold', size = 28),
#     axis.text = element_text(size = 18), 
#     legend.position = 'top',
#     # panel.background = element_rect(fill = 'grey'),
#     legend.justification = 'center', 
#     legend.key.width = unit(1.15, 'cm'),
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 16)
#   ) + 
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
#   annotate('text', x = -0.7, y = 0.25, 
#            # label = 'Higher fitness at low expression s[weak] > s[opt]', parse = T) +
#            label = 's[weak] > s[opt]', parse = T, size = 8) +
#   annotate('text', x = -0.7, y = -0.5, 
#            # label = 'Higher fitness at optimal expression sweak < sopt', parse = T)
#            label = 's[weak] < s[opt]', parse = T, size = 8)
# 
# p_fig4a_dim_int
# # ggsave(p_fig3c, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
# #        filename = 'Figures/2022-05-06_Figures_paper/Fig4_test_ddG_dim_int.pdf')
# 
# ## ddG tetramerization interface
# p_fig4a_tet_int <- data_fig_4a %>% rowwise() %>%
#   filter(!(is.na(ddG_tet_int))) %>%
#   mutate(ddG_tet_int = ifelse(ddG_tet_int > 10, 10, ddG_tet_int)) %>%
#   ggplot(aes(x = mean_s, y = mean_ds)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
#   geom_point(aes(# alpha = mut_check, 
#     size = mut_check, 
#     colour = ddG_tet_int)) +
#   # geom_errorbar(aes(ymax = mean_ds + sem_ds, ymin = mean_ds - sem_ds, 
#   #                   colour = mut_check), show.legend = F) +
#   # geom_errorbarh(aes(xmax = mean_s + sem_s, xmin = mean_s - sem_s, 
#   #                    colour = mut_check), show.legend = F) +
#   geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
#                    box.padding = 0.4, size = 4, fontface = 'bold') +
#   xlab('Sel. coeff. (optimal expression)') + 
#   scale_size_manual(values = c(3, 5, 5)) +
#   guides(size = 'none', alpha = 'none') +
#   # scale_colour_manual(values = c('grey', 'red', 'blue')) +
#   # scale_colour_gradient2(mid = 'orange', high = 'red', midpoint = 0) +
#   # scale_colour_gradient2() +
#   # scale_alpha_manual() +
#   scale_colour_viridis(option = 'A') +
#   labs(colour = '\u0394\u0394G tet. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
#   theme(
#     axis.title.x = element_text(face = 'bold', size = 20), 
#     axis.title.y = element_text(face = 'bold', size = 28),
#     axis.text = element_text(size = 18), 
#     legend.position = 'top',
#     # panel.background = element_rect(fill = 'grey'),
#     legend.justification = 'center', 
#     legend.key.width = unit(1.15, 'cm'),
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 16)
#   ) + 
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
#   annotate('text', x = -0.7, y = 0.25, 
#            # label = 'Higher fitness at low expression s[weak] > s[opt]', parse = T) +
#            label = 's[weak] > s[opt]', parse = T, size = 8) +
#   annotate('text', x = -0.7, y = -0.5, 
#            # label = 'Higher fitness at optimal expression sweak < sopt', parse = T)
#            label = 's[weak] < s[opt]', parse = T, size = 8)
# 
# p_fig4a_tet_int
# # ggsave(p_fig3c, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
# #        filename = 'Figures/2022-05-06_Figures_paper/Fig4_test_ddG_tet_int.pdf')

#### New panels with bins of ddG values ####
# Draw the figure for ddG stability
p_fig4a_stab <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_stab))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_stab), size = 3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 4, fontface = 'bold') +
  xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  labs(colour = '\u0394\u0394G subunit stability [kcal/mol]', y = expression(s[weak] - s[opt])) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 20), 
    axis.title.y = element_text(face = 'bold', size = 28),
    axis.text = element_text(size = 18), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = 's[weak] > s[opt]', parse = T, size = 8) +
  annotate('text', x = -0.7, y = -0.5, 
           label = 's[weak] < s[opt]', parse = T, size = 8) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))
p_fig4a_stab

# Draw the figure for ddG dim int
p_fig4a_dim_int <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_dim_int))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_dim_int), size = 3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 4, fontface = 'bold') +
  xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  labs(colour = '\u0394\u0394G dim. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 20), 
    axis.title.y = element_text(face = 'bold', size = 28),
    axis.text = element_text(size = 18), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = 's[weak] > s[opt]', parse = T, size = 8) +
  annotate('text', x = -0.7, y = -0.5, 
           label = 's[weak] < s[opt]', parse = T, size = 8) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))

p_fig4a_dim_int

# Draw the figure for ddG tet int
p_fig4a_tet_int <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_tet_int))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_tet_int), size =3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 4, fontface = 'bold') +
  xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  labs(colour = '\u0394\u0394G tet. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 20), 
    axis.title.y = element_text(face = 'bold', size = 28),
    axis.text = element_text(size = 18), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = 's[weak] > s[opt]', parse = T, size = 8) +
  annotate('text', x = -0.7, y = -0.5, 
           label = 's[weak] < s[opt]', parse = T, size = 8) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))
p_fig4a_tet_int

#### Repeat but with meaningful bins ####

# Add bins of ddG values
data_fig_4a %<>%
  mutate(bins_stab = cut(ddG_stab, breaks = c(-40, 0, 1, 2, 5, 100),
                         include.lowest = TRUE, na.rm = TRUE), 
         bins_dim_int = cut(ddG_dim_int, breaks = c(-40, 0, 1, 2, 5, 100),
                            include.lowest = TRUE, na.rm = TRUE),
         bins_tet_int = cut(ddG_tet_int, breaks = c(-40, 0, 1, 2, 5, 100),
                            include.lowest = TRUE, na.rm = TRUE)
  ) %>% mutate(bins_stab = as.character(bins_stab), 
               bins_dim_int = as.character(bins_dim_int), 
               bins_tet_int = as.character(bins_tet_int)) %>%
  # Rename bins
  mutate(bins_stab = ifelse(bins_stab == '[-40,0]', '< 0', 
                            ifelse(bins_stab == '(5,100]', '> 5', bins_stab)), 
         bins_dim_int = ifelse(bins_dim_int == '[-40,0]', '< 0',
                               ifelse(bins_dim_int == '(5,100]', '> 5', bins_dim_int)), 
         bins_tet_int = ifelse(bins_tet_int == '[-40,0]', '< 0',
                               ifelse(bins_tet_int == '(5,100]', '> 5', bins_tet_int))) %>%
  mutate(bins_stab = factor(bins_stab, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5')),
         bins_dim_int = factor(bins_dim_int, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5')), 
         bins_tet_int = factor(bins_tet_int, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5'))
         )

table(data_fig_4a$bins_stab)
table(data_fig_4a$bins_dim_int)
table(data_fig_4a$bins_tet_int)

# Check correlation between deltaS and destabilizing effects
corr_ddg_stab_deltaS <- data_fig_4a %>% filter(!(is.na(ddG_stab)), mut_check != 'WT') %>%
  select(ddG_stab, ddG_dim_int, ddG_tet_int, mean_ds)

cor.test(corr_ddg_stab_deltaS$ddG_stab, corr_ddg_stab_deltaS$mean_ds, method = 'spearman')
cor.test(corr_ddg_stab_deltaS$ddG_dim_int, corr_ddg_stab_deltaS$mean_ds, method = 'spearman')
cor.test(corr_ddg_stab_deltaS$ddG_tet_int, corr_ddg_stab_deltaS$mean_ds, method = 'spearman')

cor(corr_ddg_stab_deltaS, method = 'spearman')

# Draw the figure for ddG stability
p_fig4a_stab <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_stab))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_stab), size = 3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 6, fontface = 'bold') +
  xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  # scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  scale_colour_manual(values = c('#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  labs(colour = '\u0394\u0394G subunit stability [kcal/mol]', y = expression(s[weak] - s[opt])) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 24), 
    axis.title.y = element_text(face = 'bold', size = 30),
    axis.text = element_text(size = 22), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = 's[weak] > s[opt]', parse = T, size = 10) +
  annotate('text', x = -0.7, y = -0.5, 
           label = 's[weak] < s[opt]', parse = T, size = 10) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))
p_fig4a_stab

# Draw the figure for ddG dim int
p_fig4a_dim_int <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_dim_int))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_dim_int), size = 3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 6, fontface = 'bold') +
  xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  # scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  scale_colour_manual(values = c('#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  labs(colour = '\u0394\u0394G dim. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 24), 
    axis.title.y = element_text(face = 'bold', size = 30),
    axis.text = element_text(size = 22), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = 's[weak] > s[opt]', parse = T, size = 10) +
  annotate('text', x = -0.7, y = -0.5, 
           label = 's[weak] < s[opt]', parse = T, size = 10) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))

p_fig4a_dim_int

# Draw the figure for ddG tet int
p_fig4a_tet_int <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_tet_int))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_tet_int), size =3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 6, fontface = 'bold') +
  xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  # scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  scale_colour_manual(values = c('#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  labs(colour = '\u0394\u0394G tet. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 24), 
    axis.title.y = element_text(face = 'bold', size = 30),
    axis.text = element_text(size = 22), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = 's[weak] > s[opt]', parse = T, size = 10) +
  annotate('text', x = -0.7, y = -0.5, 
           label = 's[weak] < s[opt]', parse = T, size = 10) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))
p_fig4a_tet_int


#### Figure 4B: Destabilizing mutants ####

# data_compensated <- score_significant %>% rowwise() %>%
#   mutate(signif_check = ifelse(and(diffNormScore > 0.1, p.adj <= 0.05), 'Fitter at low expression', 
#                                ifelse(and(diffNormScore < -0.1, p.adj <= 0.05), 'Fitter at high expression', 
#                                       'Expression-independent effects'))) %>%
#   mutate(signif_check = factor(signif_check, levels = c('Fitter at low expression', 
#                                                         'Fitter at high expression', 
#                                                         'Expression-independent effects')))
# 
# table(data_compensated$signif_check)
# 
# data_compensated <- left_join(x = data_compensated, 
#                               y = all_data_complete %>% 
#                                 filter(Arabinose == 0.2, TMP == 10, Timepoint == 10) %>% # Avoid repetition
#                                 select(Position, WT_Residue, Residue, Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, 
#                                        Mean_ddG_int_HM_A_D), 
#                               by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue') 
# )
# 
# # Panel for stability
# p_stab <- data_compensated %>% filter(!(is.na(Mean_ddG_stab_HET))) %>%
#   ggplot(aes(x = Mean_ddG_stab_HET, colour = signif_check)) + 
#   stat_ecdf(size = 1.5) +
#   theme(axis.text.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18, angle = 90),
#         axis.title = element_text(size = 20, face = 'bold'), 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'top', 
#         legend.justification = 'center') +
#   scale_colour_manual(values = c('#8C510A', '#01665E', '#bfbfbf')) +
#   xlab('ddG on subunit stability [kcal /mol]') + ylab('Cumulative proportion') +
#   labs(colour = '')
# p_stab
# 
# # Panel for binding affinity (interface 1)
# p_int1 <- data_compensated %>% filter(!(is.na(Mean_ddG_int_HM_A_C))) %>%
#   ggplot(aes(x = Mean_ddG_int_HM_A_C, colour = signif_check)) + 
#   stat_ecdf(size = 1.5) +
#   theme(axis.text.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18, angle = 90),
#         axis.title = element_text(size = 20, face = 'bold'), 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'none') +
#   scale_colour_manual(values = c('#8C510A', '#01665E', '#bfbfbf')) +
#   xlab('ddG on binding (dim. interface) [kcal/mol]') +
#   ylab('Cumulative proportion') +
#   labs(colour = '')
# p_int1
# 
# # Panel for binding affinity (interface 2)
# p_int2 <- data_compensated %>% filter(!(is.na(Mean_ddG_int_HM_A_D))) %>%
#   ggplot(aes(x = Mean_ddG_int_HM_A_D, colour = signif_check)) + 
#   stat_ecdf(size = 1.5) +
#   theme(axis.text.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18, angle = 90),
#         axis.title = element_text(size = 20, face = 'bold'), 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'none') +
#   scale_colour_manual(values = c('#8C510A', '#01665E', '#bfbfbf')) +
#   xlab('ddG on binding (tet. interface) [kcal /mol]') +
#   ylab('Cumulative proportion') +
#   labs(colour = '')
# p_int2

ddg_sel_coeff_data <- all_data_complete %>% 
  select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP, mean_sel_coeff, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D) %>%
  filter(!(is.na(Mean_ddG_stab_HET)), Timepoint == 10, TMP == 10)

summary(ddg_sel_coeff_data$Mean_ddG_int_HM_A_C)
summary(ddg_sel_coeff_data$Mean_ddG_int_HM_A_D)
summary(ddg_sel_coeff_data$Mean_ddG_stab_HET)

# Try with bins for stability
ddg_sel_coeff_data %<>% 
  mutate(bins_stab = cut(Mean_ddG_stab_HET, 
                         quantile(Mean_ddG_stab_HET, probs = c(0, 0.25, 0.5, 0.75, 1)), 
                         include.lowest = TRUE), 
         bins_dim_int = cut(Mean_ddG_int_HM_A_C, 
                         quantile(Mean_ddG_int_HM_A_C, probs = c(0, 0.25, 0.5, 0.75, 1)), 
                         include.lowest = TRUE), 
         bins_tet_int = cut(Mean_ddG_int_HM_A_D, 
                         quantile(Mean_ddG_int_HM_A_D, probs = c(0, 0.25, 0.5, 0.75, 1)), 
                         include.lowest = TRUE)
         ) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                                  'Optimal', 'Overexpressed')))


# Use the same meaningful bins as above
ddg_sel_coeff_data %<>% 
  mutate(bins_stab = cut(Mean_ddG_stab_HET, breaks = c(-40, 0, 1, 2, 5, 100),
                         include.lowest = TRUE, na.rm = TRUE), 
         bins_dim_int = cut(Mean_ddG_int_HM_A_C, breaks = c(-40, 0, 1, 2, 5, 100),
                            include.lowest = TRUE, na.rm = TRUE),
         bins_tet_int = cut(Mean_ddG_int_HM_A_D, breaks = c(-40, 0, 1, 2, 5, 100),
                            include.lowest = TRUE, na.rm = TRUE)
  ) %>% mutate(bins_stab = as.character(bins_stab), 
               bins_dim_int = as.character(bins_dim_int), 
               bins_tet_int = as.character(bins_tet_int)) %>%
  # Rename bins
  mutate(bins_stab = ifelse(bins_stab == '[-40,0]', '< 0', 
                            ifelse(bins_stab == '(5,100]', '> 5', bins_stab)), 
         bins_dim_int = ifelse(bins_dim_int == '[-40,0]', '< 0',
                               ifelse(bins_dim_int == '(5,100]', '> 5', bins_dim_int)), 
         bins_tet_int = ifelse(bins_tet_int == '[-40,0]', '< 0',
                               ifelse(bins_tet_int == '(5,100]', '> 5', bins_tet_int))) %>%
  # mutate(bins_stab = factor(bins_stab, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5')),
  #        bins_dim_int = factor(bins_dim_int, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5')), 
  #        bins_tet_int = factor(bins_tet_int, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5'))
  # ) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                                  'Optimal', 'Overexpressed')))

# Add nonsense mutants for reference
nonsense_mutants <- all_data_complete %>% 
  filter(TMP == 10, Residue == '*', Position >= 30, Position <= 70) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA))))), 
         bins_stab = 'Stop', bins_dim_int = 'Stop', bins_tet_int = 'Stop') %>% 
  select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP, mean_sel_coeff, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D, 
         bins_stab, bins_dim_int, bins_tet_int, exp_level)

ddg_sel_coeff_data <- bind_rows(ddg_sel_coeff_data %>% ungroup() %>%
                                  select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C,
                                         -Mean_ddG_int_HM_A_D),
                                nonsense_mutants %>% ungroup() %>%
                                  select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C,
                                         -Mean_ddG_int_HM_A_D)) %>% 
  mutate(
    bins_stab = factor(bins_stab, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')),
    bins_dim_int = factor(bins_dim_int, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    bins_tet_int = factor(bins_tet_int, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                             'Optimal', 'Overexpressed'))
  )
  

# Draw the figure for subunit stability
p_fig4B_stab <- ddg_sel_coeff_data %>% 
  ggplot(aes(x = bins_stab, y = mean_sel_coeff, fill = exp_level, colour = exp_level)) + 
  # geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.2)) +
  geom_boxplot(alpha = 0.4) +
  xlab('\u0394\u0394G subunit stability [kcal /mol]') + ylab('s') +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'),
                    name = 'Expression level') +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'), 
                      name = 'Expression level') +
  theme(axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.title = element_text(size = 24, face = 'bold'), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top',
        legend.justification = 'center',
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)
        )
p_fig4B_stab

# Draw the figure for the dimerization interface
p_fig4B_dim_int <- ddg_sel_coeff_data %>% 
  ggplot(aes(x = bins_dim_int, y = mean_sel_coeff, fill = exp_level, colour = exp_level)) + 
  # geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.2)) +
  geom_boxplot(alpha = 0.4) +
  xlab('\u0394\u0394G dim. interface [kcal /mol]') + ylab('s') +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'),
                    name = 'Expression level') +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'),
                    name = 'Expression level') +
  theme(axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.title = element_text(size = 24, face = 'bold'), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'none',
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)
  )
p_fig4B_dim_int


# Draw the figure for the dimerization interface
p_fig4B_tet_int <- ddg_sel_coeff_data %>% 
  ggplot(aes(x = bins_tet_int, y = mean_sel_coeff, fill = exp_level, colour = exp_level)) + 
  # geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.2)) +
  geom_boxplot(alpha = 0.4) +
  xlab('\u0394\u0394G tet. interface [kcal /mol]') + ylab('s') +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'),
                    name = 'Expression level') +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'),
                    name = 'Expression level') +
  theme(axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22),
        axis.title = element_text(size = 24, face = 'bold'), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'none',
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)
  )
p_fig4B_tet_int


## Put the panels together
p_fig4a <- plot_grid(p_fig4a_stab, p_fig4a_dim_int, p_fig4a_tet_int, ncol = 3)

legend_fig4b <- get_legend(p_fig4B_stab + theme(legend.position = 'top', 
                                                legend.title = element_text(size = 24), 
                                                legend.text = element_text(size = 22)))
p_fig4b_panels <- plot_grid(p_fig4B_stab + theme(legend.position = 'none'),
                            p_fig4B_dim_int, p_fig4B_tet_int, ncol = 3)

p_fig4b <- plot_grid(legend_fig4b, p_fig4b_panels, nrow = 2, rel_heights = c(0.2, 1))

p_fig4 <- plot_grid(p_fig4a, p_fig4b, nrow = 2, labels = c('A', 'B'), 
                    label_size = 20, label_fontface = 'bold')
p_fig4

ggsave(p_fig4, device = cairo_pdf, width = 24, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-06_Figures_paper/Fig4.pdf')
ggsave(p_fig4, device = 'png', width = 24, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-06_Figures_paper/Fig4.png')


#### Figures on mechanisms ####

## Fitter at lower expression
p_x <- seq(from = 1, to = 100, by = 0.1)
p_y <- log2(p_x)

p_y1 <- log2(p_x + 10)

df_mechanism1 <- data.frame(induction = p_x, WT = p_y, 
                            mutant = p_y1) %>%
  pivot_longer(cols = c(WT, mutant), names_to = 'Type', 
               values_to = 'Activity')

p <- df_mechanism1 %>% 
  ggplot(aes(x = induction, y = Activity, colour = Type)) +
  geom_line(size = 1.5) +
  theme(axis.title = element_text(size = 20, face = 'bold'), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(size = 2)) +
  xlab('Induction') + ylab('Activity payoff') +
  ggtitle('Fitter at lower expression')
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-06_Figures_paper/Fig_mechanism_fitLower.pdf')

## Fitter at higher expression
p_x <- seq(from = 1, to = 100, by = 0.1)
p_y <- log2(p_x)

p_y1 <- logistic(p_x, z = log2(60), d = 50, a = 0.1)

df_mechanism2 <- data.frame(induction = p_x, WT = p_y, 
                            mutant = p_y1) %>%
  pivot_longer(cols = c(WT, mutant), names_to = 'Type', 
               values_to = 'Activity')

p <- df_mechanism2 %>% 
  ggplot(aes(x = induction, y = Activity, colour = Type)) +
  geom_line(size = 1.5) +
  theme(axis.title = element_text(size = 20, face = 'bold'), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(size = 2)) +
  xlab('Induction') + ylab('Activity payoff') +
  ggtitle('Fitter at higher expression')
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-06_Figures_paper/Fig_mechanism_fitHigher.pdf')

#### Probably not using the rest ####
       
#### Figure S13: 2ER increases abundance, expression level ####
#### compensates weakly destabilizing mutants ####

# Load the main table with the sample IDs
annotation_df <- read_delim('Data/2022-04-25_exp_increase/Cytometry/April_2022/Stats_run_21_04_22_all_samples_clean.csv',
                            delim = ';', locale = locale(decimal_mark = ','))

# Use a loop to load the rest of the data
path_files <- 'Data/2022-04-25_exp_increase/Cytometry/April_2022/All_raw_files/'
list_files <- list.files(path_files)

all_data <- c()

for(infile in list_files){
  # Load the data
  new_data <- read_delim(file.path(path_files, infile), delim = ',', col_names = T) %>%
    select(TIME, 'GRN-B-HLin', 'FSC-HLin', 'SSC-HLin', 'FSC-HLog', 'SSC-HLog')
  
  
  # Add the name of the file as an ID (remove the csv)
  new_data %<>% mutate(ID = substr(x = infile, start = 1, stop = (nchar(infile) - 4)))
  
  all_data <- bind_rows(all_data, new_data)
}

# Add the annotation data to identify each sample
# Join by well
all_data <- left_join(x = all_data %>% separate(col = ID, into = c('ID', 'Well'), sep = 'am.'), 
                      y = annotation_df,
                      by = c('Well' = 'Well'))


# Add the logarithms of the green fluorescence
all_data %<>% mutate(GRNBHLog = ifelse(log10(`GRN-B-HLin`) < 0, 0, log10(`GRN-B-HLin`)),
                     `FSC-HLog` = ifelse(log10(`FSC-HLin`) < 0, 0, log10(`FSC-HLin`)),
                     `SSC-HLog` = ifelse(log10(`SSC-HLin`) < 0, 0, log10(`SSC-HLin`))
)

all_data_processed <- all_data %>% rowwise() %>%
  mutate(circle_test = (`FSC-HLog` - 1.25)^2 + (`SSC-HLog` - 1.25)^2) %>%
  filter(circle_test < 1) %>%
  mutate(
    Arabinose = Arabinoseconcentration
  ) %>% select(-Arabinoseconcentration)

all_data_processed %<>% 
  mutate(Expression_level = ifelse(Arabinose == 0, 'Absent', 
                                   ifelse(Arabinose == 0.01, 'Weak expression', 
                                          ifelse(Arabinose == 0.2, 'Optimal expression', NA)))
  ) %>%
  mutate(Expression_level = factor(Expression_level,
                                   levels = c('Absent', 'Weak expression', 'Optimal expression')))

all_data_processed_summary <- all_data_processed %>% 
  ungroup() %>% group_by(Expression_level, Timepoint, Replicate, Name) %>%
  summarise(med_fluo = median(GRNBHLog))

# Draw the boxplots
p_fig4a_box <- all_data_processed %>%
  filter(Timepoint == 'Exponential phase') %>%
  mutate(Name = factor(Name, levels = c('DfrB1', 'DfrB1-GFP', 'DfrB1-M26-GFP',
                                        'DfrB1-2ER-GFP', 'DfrB1-2ER-M26-GFP','GFP'))) %>%
  ggplot(aes(x = Name, y = as.numeric(GRNBHLog))) +
  # geom_jitter(alpha = 0.3, width = 0.2) +
  # geom_violin(alpha = 0.6) + 
  geom_boxplot(alpha = 1, outlier.shape = NA, colour = 'black', fill = '#bfbfbf') + 
  facet_wrap(~Expression_level, scales = 'free') +
  # scale_colour_manual(values = c('#000000', # DfrB1
  #                                '#D55E00', # DfrB1-GFP
  #                                '#CC79A7', # DfrB1-M26-GFP
  #                                '#56B4E9', # DfrB1-2ER-GFP
  #                                '#0072B2', # DfrB1-2ER-M26-GFP
  #                                '#009E73' # GFP
  # )) +
  # scale_fill_manual(values = c('#000000', # DfrB1
  #                                '#D55E00', # DfrB1-GFP
  #                                '#CC79A7', # DfrB1-M26-GFP
  #                                '#56B4E9', # DfrB1-2ER-GFP
#                                '#0072B2', # DfrB1-2ER-M26-GFP
#                                '#009E73' # GFP
# )) +
theme(legend.position = 'none',
      strip.background = element_rect(fill = 'white'), strip.text = element_text(face = 'bold', size = 16),
      axis.text.x =  element_text(size = 16, angle = 45, hjust = 1), 
      axis.text.y = element_text(size = 16),
      axis.title = element_text(face = 'bold', size = 20),
      panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
      panel.grid.minor = element_blank(), 
      plot.title = element_text(hjust = 0.5), 
      axis.line=element_line(), 
      legend.justification = 'center') +
  xlab('Construct') + ylab('GFP fluorescence (log10)') +
  labs(fill = 'Construct', colour = 'Construct')
p_fig4a_box

# # Prepare the panel for the constructs
# p_fig4a_constructs <- ggdraw() + 
#   # draw_image('Figures/2022-04-19_extra_figures/Figure_S12_expressionIncrease_2ER_constructs.png')
#   draw_image('Figures/2022-04-19_extra_figures/Figure_S12_expressionIncrease_2ER_nocolors.png')
# p_fig4a_constructs
# 
# p_fig4a_final <- plot_grid(p_fig4a_constructs, p_fig4a, nrow = 2, rel_heights = c(0.5, 1))
# p_fig4a_final


#### Test the correlation between solvent accessibility and deltaS for residues at the interface ####

p <- critical_sites_deltaDMS %>% # filter(Site %in% c('Interface 1', 'Interface 2')) %>%
  ggplot(aes(x = rSASA, y = diffNormScore)) + 
  facet_wrap(~Site) +
  geom_point() +
  xlab('Relative solvent accessibility') + ylab('\u0394s') +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = 0.25, linetype = 'dashed')
p
ggsave(p, device = cairo_pdf, width = 17, height = 10, dpi = 300, 
       filename = 'Figures/2022-04-03_extra_figures/RSA_deltaS_sites.pdf')

p <- critical_sites_deltaDMS %>% ungroup() %>% 
  group_by(Position, Site) %>% 
  summarise(rSASA = mean(rSASA), sd_dS = sd(diffNormScore)) %>%
  ggplot(aes(x = rSASA, y = sd_dS)) + 
  facet_wrap(~Site) +
  geom_point() +
  xlab('Relative solvent accessibility') + ylab('Std. dev. (\u0394s)') +
  geom_vline(xintercept = 0.25, linetype = 'dashed')
p
ggsave(p, device = cairo_pdf, width = 17, height = 10, dpi = 300, 
       filename = 'Figures/2022-04-03_extra_figures/RSA_sd_deltaS_sites.pdf')

p <- critical_sites_deltaDMS %>% ungroup() %>% filter(Site != 'Disordered region') %>%
  group_by(Position, Site) %>% 
  summarise(rSASA = mean(rSASA), sd_dS = sd(diffNormScore)) %>%
  ggplot(aes(x = rSASA, y = sd_dS)) + 
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.8, size = 6,
           label.y.npc = 0.15, method = 'spearman', show.legend = F, cor.coef.name = 'rho', 
  ) +
  geom_smooth(method = 'lm', show.legend = F) +
  xlab('Relative solvent accessibility') + ylab('Std. dev. (\u0394s)') +
  geom_vline(xintercept = 0.25, linetype = 'dashed')
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-04-03_extra_figures/RSA_sd_deltaS.pdf')

#### Let's see what the maximum and minimum deltaS at each position look like when comparing to RSA ####

p <- critical_sites_deltaDMS %>% ungroup() %>% filter(Site != 'Disordered region') %>%
  group_by(Position, Site) %>% 
  summarise(rSASA = mean(rSASA), max_dS = max(diffNormScore)) %>%
  ggplot(aes(x = rSASA, y = max_dS)) + 
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.8, size = 6,
           label.y.npc = 0.15, method = 'spearman', show.legend = F, cor.coef.name = 'rho', 
  ) +
  geom_smooth(method = 'lm', show.legend = F) +
  xlab('Relative solvent accessibility') + ylab('Max (\u0394s)') +
  geom_vline(xintercept = 0.25, linetype = 'dashed')
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-04-03_extra_figures/RSA_max_deltaS.pdf')

p <- critical_sites_deltaDMS %>% ungroup() %>% filter(Site != 'Disordered region') %>%
  group_by(Position, Site) %>% 
  summarise(rSASA = mean(rSASA), min_dS = min(diffNormScore)) %>%
  ggplot(aes(x = rSASA, y = min_dS)) + 
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.8, size = 6,
           label.y.npc = 0.15, method = 'spearman', show.legend = F, cor.coef.name = 'rho', 
  ) +
  geom_smooth(method = 'lm', show.legend = F) +
  xlab('Relative solvent accessibility') + ylab('Min (\u0394s)') +
  geom_vline(xintercept = 0.25, linetype = 'dashed')
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-04-03_extra_figures/RSA_min_deltaS.pdf')


#### Test a similar figure but with secondary structures instead of protein regions ####

p <- data_fig_4_final_new %>% rowwise() %>%
  mutate(mut_id = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  ggplot(aes(x = Secondary_structure, y = diffNormScore)) +
  geom_boxplot(outlier.shape = NA, fill = 'grey') +
  geom_jitter(width = 0.2, alpha = 0.3) + 
  theme(legend.position = 'none', 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(face = 'bold', size = 20), 
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 18), 
        plot.margin = margin(0, 3.5, 0, 3.5, 'cm')) +
  xlab('') + ylab('\u0394s')
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-04-03_extra_figures/Secondary_structure_deltaS.pdf')

#### Figure 5: Relationship between fitness effects and relative solvent accessibility ####

## Load the ChimeraX figures
# ara_diff_struc <- ggdraw() + draw_image('Figures/chimerax/ara0.2_ara0.01_meanDiffEffects_new.png')
ara_diff_struc <- ggdraw() + draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/meanDiffFitEff_ara0.2_ara0.01_af2.png')
ara_diff_struc

sd_ara_diff_struc <- ggdraw() + draw_image('Figures/2022-03-24_Figures_paper/chimerax/Figures/sdDiffFitEff_ara0.2_ara0.01_af2.png')
sd_ara_diff_struc

### Load the data from the random forest results ###

rel_importance_final <- read_delim('Results_ML/diffNorm_ara0.2_ara0.01/model_diffFit_relImportancesDropCol_bestVariables.txt',
                                   delim =  '\t') %>% arrange(Importance)

## Rename the tags (new normalization method)
rel_importance_final %<>% mutate(Feature = ifelse(Feature == 'random_var', 'Random variable', 
                                                 ifelse(Feature == 'Mean_ddG_int_HM_A_D', 'ddG interface 2', 
                                                        ifelse(Feature == 'Mean_ddG_int_HM_A_C', 'ddG interface 1',
                                                               ifelse(Feature == 'polarity_zimmerman', 'Polarity Zimmerman', 
                                                                      ifelse(Feature == 'hydrophilicity_hopp', 'Hydrophilicity Hopp', 
                                                                             ifelse(Feature == 'Mean_ddG_stab_HET', 'ddG subunit stability', 
                                                                                    ifelse(Feature == 'rSASA', 'Rel. solvent accessibility', NA))))))))


## Reorder the factors (new normalization method)
rel_importance_final %<>% mutate(Feature = factor(Feature, levels = c('Random variable',
                                                                      'Polarity Zimmerman',
                                                                      'Hydrophilicity Hopp',
                                                                     'ddG subunit stability',
                                                                     'Rel. solvent accessibility'
)
)
)

p_rel_importance <- rel_importance_final %>% 
  ggplot(aes(x = Importance, y = Feature)) +
  geom_bar(stat = 'identity') +
  theme(axis.title = element_text(size = 20, face = 'bold'), 
        axis.text = element_text(size = 18)) +
  geom_vline(xintercept = 0)
p_rel_importance

## Load the data on observations vs predictions
rf_pred <- read_delim('Results_ML/diffNorm_ara0.2_ara0.01/pred_rf_bestVariables.txt',
                                   delim =  '\t')

p_rf_pred <- rf_pred %>% ggplot(aes(x = pred_data, y = test_data)) +
  geom_point(size = 3) +
  theme(axis.title = element_text(size = 19, face = 'bold'), 
        axis.text = element_text(size = 18), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  geom_smooth(method = 'lm') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.7, size = 6,
           label.y.npc = 0.15, method = 'pearson', cor.coef.name = 'r'
  ) +
  xlab('Predicted \u0394s (RF)') + ylab('Observed \u0394s (test set)')
p_rf_pred

## Draw the figures together
p_fig5_new <- plot_grid(ara_diff_struc, sd_ara_diff_struc, p_rel_importance, p_rf_pred, 
                        labels = c('A', 'B', 'C', 'D'), label_size = 20, ncol = 2)
p_fig5_new
ggsave(p_fig5_new, device = cairo_pdf, width = 17, height = 14, dpi = 300, 
       filename = 'Figures/2022-04-06_Figures_paper/Fig5_new.pdf')
ggsave(p_fig5_new, device = 'png', width = 17, height = 14, dpi = 300, 
       filename = 'Figures/2022-04-06_Figures_paper/5.Fig5_new.png')
