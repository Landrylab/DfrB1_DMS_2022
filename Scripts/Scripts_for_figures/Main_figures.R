###########################################################
####        DfrB1_DMS_main_figures                     ####
#### This script will organize the data and prepare    ####
#### the figures for the paper.                        ####
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
library(scales)

library(car)
library(ggrepel)
library(nlme)
library(lme4)
library(rlme)
library(FDRestimation)
library(gridGraphics)
library(shadowtext)
library(viridis)

library(agricolae)
library(ggtext)
library(glue)

theme_set(theme_cowplot() + 
            theme(panel.background = element_rect(fill = 'white'),
                  plot.background = element_rect(fill = 'white'))
          ) 

## Set working directory as the DfrB1_DMS_2022 home directory
setwd('/path/to/DfrB1_DMS_2022')

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

# Load the complete dataset
all_data_complete <- read_delim(
  'Data/Complete_datasets/complete_dataset_avgBothSequencers.txt', delim = '\t')

#### Figure 1: Study system ####

#### Figure 1C: Effect of arabinose on fitness ####

annotation_df <- read_delim(
  'Data/Cytometry/DfrB1_sfGFP_expression_level/Stats_run2_05_05_22_matched.csv',
  locale = locale(decimal_mark = ','), delim = ',')

# Use a loop to load the rest of the data
path_files <- 'Data/Cytometry/DfrB1_sfGFP_expression_level/All_raw_files/'
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
all_data_expression %<>% 
  mutate(GRNBHLog = ifelse(log2(`GRN-B-HLin`) < 0, 0, log2(`GRN-B-HLin`)),
                     `FSC-HLog` = ifelse(log10(`FSC-HLin`) < 0, 0, log10(`FSC-HLin`)),
                     `SSC-HLog` = ifelse(log10(`SSC-HLin`) < 0, 0, log10(`SSC-HLin`))
)

all_data_processed_expression <- all_data_expression %>% rowwise() %>%
  mutate(circle_test = (`FSC-HLog` - 1.25)^2 + (`SSC-HLog` - 1.25)^2) %>%
  filter(circle_test < 1) %>%
  mutate(
    Arabinose = str_c(toString(Arabinose), '% arabinose', sep = '')
  )

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
  group_by(Arabinose, Well) %>%
  summarise(
    med_fluo = mean(GRNBHLog),
    sem_fluo = sd(GRNBHLog) / sqrt(n()),
    num_cells = n()
    )

all_data_processed_summary_expression %<>%
  mutate(exp_level = ifelse(Arabinose == '0.2% arabinose', 'Optimal',
                            ifelse(Arabinose == '0.05% arabinose', 'Near-optimal',
                                   ifelse(Arabinose == '0.025% arabinose', 'Suboptimal',
                                          ifelse(Arabinose == '0.01% arabinose', 'Weak',
                                                 ifelse(Arabinose == '0.4% arabinose', 'Above-optimal',
                                                        'Absent')))))
  )



#### Figures for growth curve validations ####
plate.ind <- 'Data/Growth_curves/Fitness_TMP_ID_IGA_23_02_22.xlsx'
file.od <- 'Data/Growth_curves/Data_23_02_2022_bact.xlsx'

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
fig1C_data <- left_join(x = data.od1.summary.wide, 
                        y = all_data_processed_summary_expression %>% 
                          separate(Arabinose, into = c('Arabinose', 'tmp'), sep = '%') %>%
                          mutate(Arabinose = as.numeric(Arabinose)) %>% ungroup() %>%
                          group_by(Arabinose) %>% 
                          summarise(mean_fluo = round(mean(med_fluo), 2)),
                        by = c('Arabinose' = 'Arabinose'))

# Add a column for the maximum difference (0% arabinose)
max_diff <- fig1C_data %>% filter(Arabinose == 0) %>% group_by(Arabinose) %>%
  summarise(max_diff = median(auc_diff))

fig1C_data %<>% mutate(pct_recovery = 100*(`10` / `0`))

fig1C_data %<>%
  mutate(exp_level = ifelse(Arabinose == 0.2, 'Optimal', 
                            ifelse(Arabinose == 0.05, 'Near-optimal',
                                   ifelse(Arabinose == 0.025, 'Suboptimal', 
                                          ifelse(Arabinose == 0.01, 'Weak', 
                                                 ifelse(Arabinose == 0.4, 'Above-optimal', 
                                                        'Absent')))))
  )

fig1C_data_sum <- fig1C_data %>% ungroup() %>% 
  group_by(Arabinose, exp_level) %>%
  summarise(mean_recovery = mean(pct_recovery), 
            sem_recovery = sd(pct_recovery) / sqrt(n()), 
            num_samples = n())

# Add an annotation box for the legend of figure 1C
annotation_box <- roundrectGrob(x = unit(0.7, 'npc'), y = unit(0.3, 'npc'), 
                                width = unit(0.3, 'npc'), height = unit(0.4, 'npc'), 
                                gp = gpar(lwd = 2, col = 'black',
                                          # fill = '#737373',
                                          fill = '#a6a6a6',
                                          lty = 2))

p_fig1C <- fig1C_data_sum %>% ungroup() %>%
  mutate(Arabinose = as.factor(Arabinose), 
         exp_level = factor(exp_level, 
                            levels = c('Absent', 'Weak', 'Suboptimal',
                                       'Near-optimal', 'Optimal', 'Above-optimal'))
         ) %>%
  ggplot(aes(x = Arabinose, y = mean_recovery,
             ymax = mean_recovery + 2.5*sem_recovery, 
             ymin = mean_recovery - 2.5*sem_recovery, 
             colour = exp_level)
         ) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.2, size = 0.7) +
  xlab('Arabinose concentration (% m/v)') +
  ylab('Recovered growth (%)') +
  labs(colour = 'Promoter activity') +
  scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20, face = 'bold'),
        legend.text = element_text(size = 18),
        legend.position = c(0.7, 0.3),
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  annotation_custom(annotation_box)
p_fig1C

write.table(fig1C_data_sum, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/growth_recovery_pct.tsv')

## Try a new version of figure 1C ##
comps <- compare_means(pct_recovery~exp_level, data = fig1C_data,
                       paired = FALSE, method = 't.test') %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16', 
                           ifelse(p >= 0.01, str_c("p = ", round(as.numeric(p), 2), sep = ''), 
                            sprintf("p = %2.1e", as.numeric(p)))
  )
  ) %>%
  filter(or(group1 == 'Optimal', group2 == 'Optimal')) %>%
  mutate(y_pos = c(105, 115, 125, 135, 145), 
         col1 = ifelse(group1 == 'Absent', 0, 
                       ifelse(group1 == 'Weak', 0.01, 
                              ifelse(group1 == 'Suboptimal', 0.025,
                                     ifelse(group1 == 'Near-optimal', 0.05, 
                                            ifelse(group1 == 'Optimal', 0.2, 0.4))))), 
         col2 = ifelse(group2 == 'Absent', 0, 
                       ifelse(group2 == 'Weak', 0.01, 
                              ifelse(group2 == 'Suboptimal', 0.025,
                                     ifelse(group2 == 'Near-optimal', 0.05, 
                                            ifelse(group2 == 'Optimal', 0.2, 0.4)))))
         )


m <- aov(pct_recovery~exp_level, data=fig1C_data)
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'exp_level')
tukey_test2 <- TukeyHSD(m)
write.table(tukey_test2$exp_level, file = 'Data/Tukey_test_fig1C.tsv', 
            row.names = T, col.names = T, sep = '\t', append = F)

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(105, 105, 105, 105, 105, 105))
tukey_annotation$Site_tukey <- rownames(tukey_annotation)

tukey_annotation %<>% 
  mutate(Arabinose = ifelse(Site_tukey == 'Absent', 0, 
                       ifelse(Site_tukey == 'Weak', 0.01, 
                              ifelse(Site_tukey == 'Suboptimal', 0.025,
                                     ifelse(Site_tukey == 'Near-optimal', 0.05, 
                                            ifelse(Site_tukey == 'Optimal', 0.2, 0.4)))))
         )

p_fig1C <- fig1C_data %>% ungroup() %>%
  mutate(Arabinose = as.factor(Arabinose), 
         exp_level = factor(exp_level, 
                            levels = c('Absent', 'Weak', 'Suboptimal',
                                       'Near-optimal', 'Optimal', 'Above-optimal'))
  ) %>%
  ggplot(aes(x = Arabinose, y = pct_recovery, colour = exp_level)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  geom_point(data = fig1C_data_sum, size = 3, inherit.aes = F, 
             aes(colour = exp_level, x = as.factor(Arabinose), y = mean_recovery)) +
  geom_errorbar(width = 0.2, size = 0.7, data = fig1C_data_sum, inherit.aes = F, 
                aes(colour = exp_level, x = as.factor(Arabinose), 
                    ymax = mean_recovery + 2.5*sem_recovery, 
                    ymin = mean_recovery - 2.5*sem_recovery)) +
  xlab('Arabinose concentration (% m/v)') +
  ylab('Recovered growth (%)') +
  labs(colour = 'Promoter activity') +
  scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18),
        legend.title = element_text(size = 20, face = 'bold'),
        legend.text = element_text(size = 18),
        legend.position = c(0.7, 0.3),
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  annotation_custom(annotation_box) +
  geom_text(data = tukey_annotation, size = 5, inherit.aes = F,
            aes(label = groups, x = as.factor(Arabinose), y = y))
p_fig1C

#### Figure 1D: Showing the medians of the distribution of GFP fluorescence ####

all_data_processed_expression %<>% ungroup() %>% 
  mutate(
    exp_level = ifelse(Arabinose == '0.2% arabinose', 'Optimal', 
                       ifelse(Arabinose == '0.05% arabinose', 'Near-optimal',
                              ifelse(Arabinose == '0.025% arabinose', 'Suboptimal', 
                                     ifelse(Arabinose == '0.01% arabinose', 'Weak', 
                                            ifelse(Arabinose == '0.4% arabinose', 'Above-optimal',
                                                   'Absent')))))
  ) %>%
  mutate(exp_level = factor(exp_level, levels = c('Absent', 'Weak', 'Suboptimal',
                                                  'Near-optimal', 'Optimal', 'Above-optimal')))

m <- aov(GRNBHLog~exp_level, data=all_data_processed_expression)
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'exp_level')
tukey_test2 <- TukeyHSD(m)
write.table(tukey_test2$exp_level, file = 'Data/Tukey_test_fig1D.tsv', 
            row.names = T, col.names = T, sep = '\t', append = F)

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(10, 10, 10, 10, 10, 10))
tukey_annotation$Site_tukey <- rownames(tukey_annotation)

tukey_annotation %<>% 
  mutate(Arabinose = ifelse(Site_tukey == 'Absent', 0, 
                            ifelse(Site_tukey == 'Weak', 0.01, 
                                   ifelse(Site_tukey == 'Suboptimal', 0.025,
                                          ifelse(Site_tukey == 'Near-optimal', 0.05, 
                                                 ifelse(Site_tukey == 'Optimal', 0.2, 0.4)))))
  )

p_fig1D <- all_data_processed_expression %>%
  ggplot(aes(x = exp_level, y = GRNBHLog, 
             colour = exp_level, group = Well)) + 
  geom_violin(aes(fill = exp_level), alpha = 0.6, scale = 'width',
              position = 'dodge', width = 0.8) +
  geom_boxplot(aes(fill = exp_level), alpha = 0.4, position = 'dodge', width = 0.8) +
  scale_colour_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  scale_fill_manual(values = c('grey', '#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  xlab('Promoter activity') + 
  ylab('GFP fluorescence (log2)') +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20, face = 'bold'), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'none') +
  geom_text(data = tukey_annotation, size = 5, inherit.aes = F,
            aes(label = groups, x = Site_tukey, y = y))
p_fig1D

all_data_processed_expression_summary <- all_data_processed_expression %>%
  group_by(Arabinose) %>%
  summarise(med_fluo = median(`GRN-B-HLin`))

#### Panels E and F ####

expression_level <- fig1C_data %>% ungroup() %>%
  group_by(Arabinose, exp_level) %>%
  summarise(mean_fluo = mean(mean_fluo)) %>%
  mutate(exp_level = factor(exp_level, levels = c('Absent', 'Weak', 'Suboptimal', 'Near-optimal', 'Optimal', 'Above-optimal')))

write.table(expression_level, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Data/expression_levels.tsv')

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

annotation_box <- roundrectGrob(x = unit(0.205, 'npc'), y = unit(0.755, 'npc'), 
                                width = unit(0.26, 'npc'), height = unit(0.3, 'npc'), 
                                gp = gpar(lwd = 2, col = 'black',
                                          # fill = '#737373',
                                          fill = '#a6a6a6',
                                          lty = 2))

p_fig1e <- data_fig_1e_final_exp_plot %>%
  ggplot(aes(x = mean_sel_coeff, y = mean_sel_coeff_2, 
             colour = as.factor(exp_level))
  ) +
  geom_point(alpha = 1, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'black') +
  labs(colour = 'Expression level in y-axis', 
       x = expression(paste(bolditalic('s '), bold('(optimal expression)'), sep = ' ')), 
       y = expression(paste(bolditalic('s '), bold('(non-optimal expression)'), sep = ' '))
       ) +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', 'black')) +
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
all_data_all_reps <- read_delim(
  'Data/Complete_datasets/all_data_all_reps_bothSequencers.txt', 
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
  
  if(and(all(!(is.na(data_genotype$sel_coeff))), # Make sure there are no NAs
         length(unique(data_genotype$Arabinose)) == 5) # Make sure they are represented at all expression levels
  ){
  
    m <- aov(sel_coeff ~ Arabinose, data=data_genotype)
    anova_test <- anova(m)
    p_value <- anova_test$`Pr(>F)`[1]
    
    new_row <- data.frame(Genotype = c(genotype), p_val = c(p_value))
    all_anovas <- bind_rows(all_anovas, new_row)
    }
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
  mutate(mut_check_diff = ifelse(abs(diffNormScore) > 0.1, TRUE, FALSE))

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

#### Use k-means to cluster the data from Figure 1F ####

data_fig1f_wide <- data_fig1f_exp %>% ungroup() %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, exp_level) %>%
  group_by(Position, WT_Residue, Residue) %>%
  pivot_wider(names_from = exp_level, values_from = mean_sel_coeff)

## Make sure WT appears only once
wt_row <- data_fig1f_wide %>% filter(WT_Residue == Residue, Position == 2)

data_kmeans <- bind_rows(data_fig1f_wide %>% filter(WT_Residue != Residue), 
                         wt_row)

lines_remove <- c()
## Remove lines that have NAs
for(i in 1:nrow(data_kmeans)){
  if(any(is.na(data_kmeans[i,]))){
    lines_remove <- c(lines_remove, i)
  }
}

if(!(is.null(lines_remove))){
  data_kmeans <- data_kmeans[-lines_remove, ]
}

set.seed(100)
kmeans_ss <- c()
for(k in 1:10){
  print(k)
  kmeans_fitness <- kmeans(data_kmeans %>% ungroup() %>% rowwise() %>% 
                             filter(WT_Residue != Residue) %>%
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

## Use k = 4 as a good compromise between parsimony and interpretability
k = 4
kmeans_fitness <- kmeans(data_kmeans %>% ungroup() %>% rowwise() %>% 
                           select(-WT_Residue, -Position, -Residue), 
                         centers = k, iter.max = 10, nstart = 25
)

kmeans_fitness$centers

# Rearrange the clusters automatically to ensure the cluster IDs are always
# the same.
temp_clusters <- kmeans_fitness$centers
temp_clusters %<>% as.data.frame() %>% mutate(old_cluster = row_number())

cluster_relabel <- temp_clusters %>% arrange(Weak) %>% mutate(new_cluster = row_number())


data_fig1f_new <- data_kmeans %>% ungroup() %>% 
  mutate(cluster = kmeans_fitness$cluster) %>%
  pivot_longer(cols = c('Weak', 'Suboptimal', 'Near-optimal', 
                        'Optimal', 'Above-optimal'), names_to = 'exp_level', 
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
                                       'Optimal', 'Above-optimal'))) %>%
  ggplot(aes(x = exp_level, y = mean_sel_coeff, fill = exp_level)) +
  geom_violin(alpha = 0.4) +
  stat_summary(fun="mean", geom="point", size = 5, colour = 'black') +
  geom_point() +
  geom_line(aes(group = ID, colour = as.factor(cluster)), alpha = 0.3) +
  scale_colour_manual(
    values = c('#1b9e77', '#ffbf80', '#7570b3', '#e7298a')
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
  labs(colour = 'Cluster', y = expression(bolditalic('s'))) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
p_fig1f

## Check the number of mutations in each cluster
cluster_count_check <- data_fig1f_new %>% filter(exp_level == 'Weak')
table(cluster_count_check$cluster)

#### Draw figure 1 ####

p_fig1A <- ggdraw() + draw_image('Figures/Main_figures/Fig1A_DfrB1.png')
p_fig1B <- ggdraw() + draw_image('Figures/Main_figures/Fig1B_DfrB1.png')

p_fig1AB <- plot_grid(p_fig1A, p_fig1B, ncol = 2, 
                      labels = c('A', 'B'), label_fontface = 'bold', 
                      label_size = 20)


p_fig1CD <- plot_grid(p_fig1C, p_fig1D, ncol = 2,
                      labels = c('C', 'D'), label_fontface = 'bold', 
                      label_size = 20)

## Keep figure 1-A-D, change scale in Fig 1D from log10 to log2
p_fig1 <- plot_grid(p_fig1AB, p_fig1CD, nrow = 2, 
                    rel_heights = c(0.5, 1))
p_fig1

ggsave(p_fig1, device = cairo_pdf, width = 20, height = 10, dpi = 300, 
       filename = 'Figures/Main_figures/1.Fig1.pdf')
       

ggsave(p_fig1, device = 'png', width = 20, height = 10, dpi = 300, 
       filename = 'Figures/Main_figures/1.Fig1.png')
       


#### Figure 2 ####

p_fig2a_new <- data_fig_1e_final_exp_plot %>%
  ggplot(aes(x = mean_sel_coeff, y = mean_sel_coeff_2, 
             colour = as.factor(exp_level))
  ) +
  geom_point(alpha = 1, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'black') +
  labs(
    x = expression(paste(bolditalic('s '), bold('(optimal promoter activity)'), sep = ' ')),
    y = expression(paste(bolditalic('s '), bold('(non-optimal promoter activity)'), sep = ' '))
  ) +
  geom_errorbarh(aes(xmax = mean_sel_coeff + sem_sel_coeff, xmin = mean_sel_coeff - sem_sel_coeff), 
                 alpha = 0.5) +
  geom_errorbar(aes(ymax = mean_sel_coeff_2 + sem_sel_coeff_2, ymin = mean_sel_coeff_2 - sem_sel_coeff_2), 
                alpha = 0.5) +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', 'black')) +
  facet_wrap(~exp_level, scales = 'free') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, colour = 'black',
           label.x.npc = 0.05,
           label.y.npc = 0.9,
           size = 8,
           vjust = 1,
           method = 'spearman', show.legend = F, cor.coef.name = 'rho') +
  geom_smooth(method = 'loess', show.legend = F) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'none', 
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        strip.text = element_text(size = 18, face = 'bold'), 
        strip.background = element_rect(fill = "#d9d9d9")
  ) + xlim(-1, 0.5) + ylim(-1, 0.5) +
  guides(fill = 'none')
p_fig2a_new  

## New figure 2b to show the violin plots and the clusters
p_fig2b_new <- data_fig1f_new %>% rowwise() %>%
  mutate(ID = str_c(WT_Residue, Position, Residue, sep = ''), 
         exp_level = factor(exp_level,
                            levels = c('Weak', 'Suboptimal', 'Near-optimal',
                                       'Optimal', 'Above-optimal'))) %>%
  ggplot(aes(x = exp_level, y = mean_sel_coeff, fill = exp_level, colour = exp_level)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_violin(alpha = 0.4, scale = 'width') +
  stat_summary(fun="mean", geom="point", size = 5, 
               colour = 'blue') +
  stat_summary(fun="mean", geom="point", size = 5, 
               colour = 'white', shape = 4) +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'), 
                    guide = 'none') +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'), 
                    guide = 'none') +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', 
        axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.justification = 0.5) +
  xlab('Promoter activity') + 
  labs(y = expression(bolditalic('s')))
p_fig2b_new


p_fig2c_new <- data_fig1f_new %>% rowwise() %>%
  mutate(ID = str_c(WT_Residue, Position, Residue, sep = ''), 
         exp_level = factor(exp_level,
                            levels = c('Weak', 'Suboptimal', 'Near-optimal',
                                       'Optimal', 'Above-optimal'))) %>%
  ggplot(aes(x = exp_level, y = mean_sel_coeff, fill = exp_level)) +
  # geom_violin(alpha = 0.4, scale = 'width') +
  stat_summary(fun="mean", geom="point", size = 5, colour = 'black') +
  geom_point() +
  geom_line(aes(group = ID, colour = as.factor(cluster)), alpha = 0.3) +
  scale_colour_manual(
    values = c('#1b9e77', '#ffbf80', '#7570b3', '#e7298a')
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
  xlab('Promoter activity') +
  labs(colour = 'Cluster', y = expression(bolditalic('s'))) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
p_fig2c_new

# Put the figures together
p_fig2b_2c <- plot_grid(p_fig2b_new, p_fig2c_new, ncol = 2,
                        labels = c('B', 'C'), label_size = 20, label_fontface = 'bold')

p_fig2_new <- plot_grid(p_fig2a_new, p_fig2b_2c, nrow = 2, 
                        labels = c('A', ''), label_size = 20, label_fontface = 'bold', 
                        rel_heights = c(2, 1))

ggsave(p_fig2_new, device = cairo_pdf, width = 20, height = 16, dpi = 300, 
       filename = 'Figures/Main_figures/2.Fig2.pdf')


ggsave(p_fig2_new, device = 'png', width = 20, height = 16, dpi = 300, 
       filename = 'Figures/Main_figures/2.Fig2.png')

#### Figure 3: The fitness landscape at different transcription levels ####

#### Prepare annotation for the bottom of the heatmaps ####

# Load the interface data for the annotation
data_interfaces <- read_delim('Data/interface_data_A.txt', delim = '\t')
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
sec_struc <- read_delim('Data/Structural_data/DHFR_2rk1_DSSP_table_bio_rSASA.txt', delim = '\t')

# Add the data on buried residues
data_interfaces_final <- left_join(x = data_interfaces_final, 
                                   y = sec_struc %>% 
                                     mutate(Buried = ifelse(Solvent_accessibility == 'Missing', 0, 
                                                            ifelse(rSASA <= 0.25, 1, 0))) %>%
                                     select(Position, Buried), 
                                   by = c('Position' = 'Position'))

# Save the annotation for the heatmaps
write.table(data_interfaces_final, file = 'Data/data_annotation_2.txt', quote = F, 
            sep = '\t', row.names = F, col.names = T)


data_fig_3 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, Entropy)

# Relevel amino acids so that they appear in the correct order
aa_order <- c('*', 'G','A','V', 'L', 'I', 'M', 'C',
              'P', 'W', 'F', 'Y', 'S', 'T', 'N', 'Q', 
              'H', 'K', 'R', 'D', 'E')

data_fig_3$Residue <- factor(data_fig_3$Residue, levels = rev(aa_order))

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
ara_0.01 <- data_fig_3 %>%
  filter(Arabinose == 0.01) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.01_final <- as.matrix(ara_0.01 %>% select(-Position))

rownames(ara_0.01_final) <- ara_0.01$Position

residue_list <- c('A', 'R', 'D', 'N', 'C', 
                 'E', 'Q', 'G', 'H', 'I', 
                 'L', 'K', 'M', 'F', 'P',
                 'S', 'T', 'W', 'Y', 'V', '*')
ara_0.01_bool_tmp <- data_fig_3 %>% ungroup() %>%
  select(Position, WT_Residue) %>%
  unique()

residue_options <- as.data.frame(cbind(rep(ara_0.01_bool_tmp$Position, each = length(residue_list)), 
                         rep(residue_list, nrow(ara_0.01_bool_tmp))))
colnames(residue_options) <- c('Position', 'Residue')

ara_0.01_bool <- left_join(x = ara_0.01_bool_tmp, 
                           y = residue_options %>% mutate(Position = as.numeric(Position)), 
                           by = c('Position' = 'Position')) %>%
  arrange(Position, Residue) %>% 
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue) %>%
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

p3_ara0.01_new <- Heatmap(
  t(ara_0.01_final), cluster_columns = F, cluster_rows = F, 
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
  top_annotation = ha1,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[j,i]){
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

p3_ara0.01_new

## 0.025% arabinose ##

# Need to convert to wide formatted data
ara_0.025 <- data_fig_3 %>%
  filter(Arabinose == 0.025) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.025_final <- as.matrix(ara_0.025 %>% select(-Position))

rownames(ara_0.025_final) <- ara_0.025$Position

# Need to reorder the columns in the matrices
ara_0.025_final <- ara_0.025_final[1:nrow(ara_0.025_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p3_ara0.025_new <- Heatmap(
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
p3_ara0.025_new

## Arabinose 0.05%
# Need to convert to wide formatted data
ara_0.05 <- data_fig_3 %>%
  filter(Arabinose == 0.05) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.05_final <- as.matrix(ara_0.05 %>% select(-Position))

rownames(ara_0.05_final) <- ara_0.05$Position

# Need to reorder the columns in the matrices
ara_0.05_final <- ara_0.05_final[1:nrow(ara_0.05_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p3_ara0.05_new <- Heatmap(
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

p3_ara0.05_new

## Arabinose 0.2%
# Need to convert to wide formatted data
ara_0.2 <- data_fig_3 %>%
  filter(Arabinose == 0.2) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.2_final <- as.matrix(ara_0.2 %>% select(-Position))

rownames(ara_0.2_final) <- ara_0.2$Position

# Need to reorder the columns in the matrices
ara_0.2_final <- ara_0.2_final[1:nrow(ara_0.2_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p3_ara0.2_new <- Heatmap(
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
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (ara_0.01_bool_final[j,i]){
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

p3_ara0.2_new

## Arabinose 0.4%
# Need to convert to wide formatted data
ara_0.4 <- data_fig_3 %>%
  filter(Arabinose == 0.4) %>%
  select(-Entropy, -WT_Residue, -Arabinose) %>%
  pivot_wider(names_from = Residue, values_from = mean_sel_coeff)

# Need to convert the dataframe to a matrix
ara_0.4_final <- as.matrix(ara_0.4 %>% select(-Position))

rownames(ara_0.4_final) <- ara_0.4$Position

# Need to reorder the columns in the matrices
ara_0.4_final <- ara_0.4_final[1:nrow(ara_0.4_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

# Prepare the annotation on interfaces
ha2 <- HeatmapAnnotation(
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
  col = list(
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


p3_ara0.4_new <- Heatmap(
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
    if (ara_0.01_bool_final[j,i]){
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
p3_ara0.4_new

# Put the heatmaps together
ht_list = p3_ara0.01_new %v% p3_ara0.025_new%v% p3_ara0.05_new %v% p3_ara0.2_new %v% p3_ara0.4_new
p3_final <- grid.grabExpr(
  draw(ht_list, row_title_gp = gpar(fontsize=20, fontface = 'bold'),
       ht_gap = unit(1, "cm"))
)

#### New labels ####

## Arabinose 0.01
text_fig_ara0.01 <- ggplot() + 
  draw_label('Weak\npromoter activity',
    x = 0.7, y = 0.3,
    fontface = 'bold', size = 35, angle = 90, colour = '#fed976') +
  theme(axis.line = element_blank())
text_fig_ara0.01

# Arabinose 0.025
text_fig_ara0.025 <- ggplot() + draw_label(
  'Suboptimal\npromoter activity',
  x = 0.7, y = 0.4,
  fontface = 'bold', size = 35, angle = 90, colour = '#fd8d3c') +
  theme(axis.line = element_blank())
text_fig_ara0.025

# Arabinose 0.05
text_fig_ara0.05 <- ggplot() + draw_label(
  'Near-optimal\npromoter activity',
  x = 0.7, y = 0.4,
  fontface = 'bold', size = 35, angle = 90, colour = '#bd0026') +
  theme(axis.line = element_blank())
text_fig_ara0.05

# Arabinose 0.2
text_fig_ara0.2 <- ggplot() + draw_label(
  'Optimal\npromoter activity',
  x = 0.7, y = 0.4,
  fontface = 'bold', size = 35, angle = 90, colour = '#80001a') +
  theme(axis.line = element_blank())
text_fig_ara0.2

# Arabinose 0.4
text_fig_ara0.4 <- ggplot() + draw_label(
  'Above-optimal\npromoter activity',
  x = 0.7, y = 0.35,
  fontface = 'bold', size = 35, angle = 90, colour = 'black') +
  theme(axis.line = element_blank())
text_fig_ara0.4

p_fig3_text <- plot_grid(text_fig_ara0.01, text_fig_ara0.025, 
                         text_fig_ara0.05, text_fig_ara0.2,
                         text_fig_ara0.4,
                         nrow = 5, rel_heights = c(1.1, 0.85, 1, 1, 1.1))

## Save the heatmaps and their titles separately
ggsave(p3_final, width = 21, height = 31.5, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/Main_figures/3.Fig3_new_all_heatmaps_buried_nolabels.pdf')

ggsave(p_fig3_text, width = 3, height = 31.5, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/Main_figures/3.Fig3_new_all_text_labels.pdf')


#### Figure 4 ####

## Figure 4A: Distribution of deltaS weak ##

data_fig_4 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, Secondary_structure, rSASA)

# Separate the data for ara 0.2
data_part_1 <- data_fig_4 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_4 %>%
  filter(Arabinose == 0.01)

## Subtract the scores
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2",
                           "Arabinose_2", "Secondary_structure", "rSASA")

data_fig_4_final_new <- left_join(x = data_part_1, y = data_part_2, 
                                  by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue',
                                         'Residue' = 'Residue', 'rSASA' = 'rSASA',
                                         'Secondary_structure' = 'Secondary_structure')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)


data_fig_4_final_new <- left_join(x = data_part_1, y = data_part_2, 
                                  by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue', 
                                         'rSASA' = 'rSASA', 'Secondary_structure' = 'Secondary_structure')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff, 
         Mutation = str_c(WT_Residue, Position, Residue)
         )

data_labels_e2v_e2r <- data_fig_4_final_new %>%
  filter(Mutation %in% c('E2R', 'E2V'))

data_fig_4_final_new %<>%
  mutate(mut_check = Mutation %in% c('E2R', 'E2V'))

p_fig4a <- ggdensity(data = data_fig_4_final_new,
                     x = 'diffNormScore',
                     fill = 'lightgray', rug = TRUE) +
  geom_label_repel(data = data_labels_e2v_e2r, inherit.aes = F, y = 0,
             aes(x = diffNormScore, label = Mutation), 
             box.padding = 2, size = 7,
             direction = 'y',
             position = position_nudge_repel(x = 0, y = 100)) +
  theme(axis.title = element_text(size = 20, face = 'bold'), 
        axis.text = element_text(size = 18), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()
  ) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  ylab('Density') +
  labs(x = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                            bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                            sep = ''))
  )
p_fig4a

## Figure 4B ##

# Load the main table with the sample IDs
annotation_df <- read_delim('Data/Cytometry/E2V_E2R_expression_increase/Stats_run_21_04_22_all_samples_clean.csv',
                            delim = ';', locale = locale(decimal_mark = ','))

# Use a loop to load the rest of the data
path_files <- 'Data/Cytometry/E2V_E2R_expression_increase/All_raw_files/'
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

## Load the data on enzymatic activity
infile = 'Data/DfrB1_activity/WT,E2R,E2V activity.xlsx'
activity_data <- read.xlsx(infile, sheetIndex = 1, rowIndex = 4:6, colIndex = 1:11,
                           stringsAsfactors = FALSE, header = F)
colnames(activity_data) <- c('Mutant', 'Slope_rep1', 'Slope_rep2', 'Slope_rep3', 'Avg_slopes', 
                             'Sd_slopes', 'Rel_activity_rep1', 'Rel_activity_rep2',
                             'Rel_activity_rep3', 'Avg_rel_activity', 'Sd_rel_activity')
activity_data <- activity_data[, 1:11]

# Calculate the standard error of the mean
activity_data %<>% mutate(sem_rel_activity = Sd_rel_activity / sqrt(3))

activity_data_long <- activity_data %>% ungroup() %>% 
  select(Mutant, Rel_activity_rep1, Rel_activity_rep2, Rel_activity_rep3) %>%
  pivot_longer(cols = c(Rel_activity_rep1, Rel_activity_rep2, Rel_activity_rep3),
               names_to = 'Replicate', values_to = 'Rel_activity')

# Calculate p values for comparisons
comps <- compare_means(Rel_activity~Mutant, 
                       data = activity_data_long,
                       paired = F, 
                       method = 't.test') %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
                           ifelse(p > 0.01, str_c('p = ', round(as.numeric(p), 2), sep = ''),
                                  sprintf("p = %2.1e", as.numeric(p)))
  ),
  y_pos = c(270, 295, 320)
  )

p_fig4b <- activity_data_long %>% 
  mutate(Mutant = factor(Mutant, levels = c('WT', 'E2R', 'E2V'))) %>%
  ggplot(aes(x = Mutant, y = Rel_activity)) +
  geom_jitter(colour = '#737373', stat = 'identity',
              width = 0.1, alpha = 0.7, size = 2) +
  geom_point(data = activity_data, inherit.aes = F,
             aes(x = Mutant, y = Avg_rel_activity), colour = 'black', 
             size = 7) +
  geom_errorbar(data = activity_data, inherit.aes = F, width = 0.2,
                aes(ymax = Avg_rel_activity + sem_rel_activity,
                    ymin = Avg_rel_activity - sem_rel_activity, 
                    x = Mutant)) +
  theme(axis.title = element_text(size = 20, face = 'bold'), 
        axis.text = element_text(size = 18), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  geom_signif(data = as.data.frame(comps), inherit.aes = FALSE,
              aes(xmin = group1, xmax = group2, annotations=p.format, y_position = y_pos), 
              manual = TRUE, textsize = 5) +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  xlab('') + ylab('Relative activity (%)') +
  ylim(0, 350)
p_fig4b

## Figure 4C ##

all_data_processed_final <- all_data_processed  %>%
  filter(Timepoint == 'Exponential phase') %>%
  ## Adjust names
  mutate(Name = ifelse(Name == 'DfrB1-2ER-GFP', 'pBAD-dfrB1(E2R)-sfGFP', 
                       ifelse(Name == 'DfrB1-2ER-M26-GFP', 'pBAD-dfrB1[1-25](E2R)-sfGFP',
                              ifelse(Name == 'DfrB1', 'pBAD-dfrB1',
                                     ifelse(Name == 'DfrB1-GFP', 'pBAD-dfrB1-sfGFP',
                                            ifelse(Name == 'DfrB1-M26-GFP', 'pBAD-dfrB1[1-25]-sfGFP',
                                                   ifelse(Name == 'GFP', 'pBAD-sfGFP', Name))))))) %>%
  mutate(Name = factor(Name, levels = c('pBAD-dfrB1', 'pBAD-dfrB1-sfGFP', 
                                        'pBAD-dfrB1[1-25]-sfGFP',
                                        'pBAD-dfrB1(E2R)-sfGFP',
                                        'pBAD-dfrB1[1-25](E2R)-sfGFP',
                                        'pBAD-sfGFP'
  )),
  Arabinose = factor(Arabinose, levels = c(0, 0.01, 0.2)))

## Will need to plot the three panels separately
## Arabinose 0
data_ara0 <- all_data_processed_final %>% filter(Arabinose == 0)

## Do the ANOVA
m <- aov(GRNBHLog ~ Name, data=data_ara0)
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'Name')
tukey_test2 <- TukeyHSD(m)
## Save table
write.table(tukey_test2$Name, append = F, quote = F, sep = '\t', row.names = T, col.names = T, 
            file = 'Data/Tukey_FigS13B_ara0.tsv')

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5))
tukey_annotation$Site_tukey <- rownames(tukey_annotation)

p_fig4c_new_0 <-data_ara0 %>%
  ggplot(aes(x = Name, y = as.numeric(GRNBHLog))) +
  geom_boxplot(colour = 'grey', fill = 'grey', alpha = 0.4, outlier.shape = NA) + 
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x =  element_text(size = 16, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 16),
        axis.title = element_text(face = 'bold', size = 20),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.line=element_line(), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18),
        legend.justification = 'center') +
  xlab('') +
  ylab('GFP fluorescence (log10)') +
  ggtitle('0% arabinose') +
  geom_text(data = tukey_annotation,
            inherit.aes = F,
            aes(x = Site_tukey, y = y, label = groups), 
            size = 7)
p_fig4c_new_0

## Arabinose 0.01
data_ara0.01 <- all_data_processed_final %>% filter(Arabinose == 0.01)

## Do the ANOVA
m <- aov(GRNBHLog ~ Name, data=data_ara0.01)
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'Name')
tukey_test2 <- TukeyHSD(m)
## Save table
write.table(tukey_test2$Name, append = F, quote = F, sep = '\t', row.names = T, col.names = T, 
            file = 'Data/Tukey_FigS13B_ara0.01.tsv')

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5))
tukey_annotation$Site_tukey <- rownames(tukey_annotation)

p_fig4c_new_0.01 <-data_ara0.01 %>%
  ggplot(aes(x = Name, y = as.numeric(GRNBHLog))) +
  geom_boxplot(colour = '#FED976', fill = '#FED976', alpha = 0.4, outlier.shape = NA) + 
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x =  element_text(size = 16, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 16),
        axis.title = element_text(face = 'bold', size = 20),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.line=element_line(), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18),
        legend.justification = 'center') +
  xlab('') + 
  ylab('GFP fluorescence (log10)') +
  ggtitle('0.01% arabinose') +
  geom_text(data = tukey_annotation,
            inherit.aes = F,
            aes(x = Site_tukey, y = y, label = groups), 
            size = 7)
p_fig4c_new_0.01

## Arabinose 0.2
data_ara0.2 <- all_data_processed_final %>% filter(Arabinose == 0.2)

## Do the ANOVA
m <- aov(GRNBHLog ~ Name, data=data_ara0.2)
anova_test <- anova(m)
tukey_test <- HSD.test(m, trt = 'Name')
tukey_test2 <- TukeyHSD(m)
## Save table
write.table(tukey_test2$Name, append = F, quote = F, sep = '\t', row.names = T, col.names = T, 
            file = 'Data/Tukey_FigS13B_ara0.02.tsv')

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5))
tukey_annotation$Site_tukey <- rownames(tukey_annotation)

p_fig4c_new_0.2 <-data_ara0.2 %>%
  ggplot(aes(x = Name, y = as.numeric(GRNBHLog))) +
  geom_boxplot(colour = '#80001A', fill = '#80001A', alpha = 0.4, outlier.shape = NA) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.x =  element_text(size = 16, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 16),
        axis.title = element_text(face = 'bold', size = 20),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        axis.line=element_line(), 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18),
        legend.justification = 'center') +
  xlab('') +
  ylab('GFP fluorescence (log10)') +
  ggtitle('0.2% arabinose') +
  geom_text(data = tukey_annotation,
            inherit.aes = F,
            aes(x = Site_tukey, y = y, label = groups), 
            size = 7)
p_fig4c_new_0.2

p_fig4c <- plot_grid(p_fig4c_new_0, p_fig4c_new_0.01, p_fig4c_new_0.2, ncol = 3)

## Figure 4D ##

# Load the growth curves with TMP
plate.ind <- 'Data/Growth_curves/Index_growth_curves_11_05_2022.xlsx'
file.od <- 'Data/Growth_curves/Growth_curves_11_05_2022.xlsx'

## Define function to read plate data
# function to process plates ----------------------------------------------
read.my.gc <- function(file, plate.index){
  pl <- read.xlsx(file,sheetIndex = 1, rowIndex = 3:62, stringsAsfactors = FALSE,
                  header = F)
  ind <- read.xlsx(plate.index, sheetIndex = 1, rowIndex = 1:64, header = T) 
  
  time <- seq(0,0.30*(ncol(pl)-2), 0.30)
  
  colnames(pl)[1] <- "Well"
  colnames(pl)[2:ncol(pl)] <- time
  pl %<>% select(1:(ncol(pl)-2)) 
  
  data.pl <- gather(pl, key = "time", value = "OD",2:ncol(pl), convert = F)
  data.pl$time <- as.numeric(data.pl$time)
  data.pl$OD <- as.numeric(data.pl$OD)
  data.pl %<>% left_join(ind, by = "Well")
  
  # Subtract the blank for this experiment and multiply by 5 to make it OD / mL
  # (experiment was carried out in 0.2 mL)
  data.pl %<>% mutate(OD = (OD - 0.085) * 5)
  
  d <-  select(data.pl, 1:3)
  d %<>% spread(key = "Well", value = "OD", convert = F) 
  colnames(d)[1] <- "time"
  
  ## Use Growthcurver
  gc_out <- SummarizeGrowthByPlate(d, t_trim =13.5)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
}

data.od1 <- read.my.gc(file.od, plate.ind)

## Load the data without TMP
plate.ind <- 'Data/Growth_curves/Index_growth_curves_20_05_2022.xlsx'
file.od <- 'Data/Growth_curves/Growth_curves_20_05_2022.xlsx'

read.my.gc <- function(file, plate.index){
  pl <- read.xlsx(file,sheetIndex = 1, rowIndex = 3:62, stringsAsfactors = FALSE,
                  header = F)
  ind <- read.xlsx(plate.index, sheetIndex = 1, rowIndex = 1:64, header = T) 
  
  time <- seq(0,0.30*(ncol(pl)-2), 0.30)
  
  colnames(pl)[1] <- "Well"
  colnames(pl)[2:ncol(pl)] <- time
  pl %<>% select(1:(ncol(pl)-2)) 
  
  data.pl <- gather(pl, key = "time", value = "OD",2:ncol(pl), convert = F)
  data.pl$time <- as.numeric(data.pl$time)
  data.pl$OD <- as.numeric(data.pl$OD)
  data.pl %<>% left_join(ind, by = "Well")
  
  # Subtract the blank for this experiment and multiply by 5 to make it OD / mL
  # (experiment was carried out in 0.2 mL)
  data.pl %<>% mutate(OD = (OD - 0.089) * 5)
  
  d <-  select(data.pl, 1:3)
  d %<>% spread(key = "Well", value = "OD", convert = F) 
  colnames(d)[1] <- "time"
  
  ## Use Growthcurver
  gc_out <- SummarizeGrowthByPlate(d, t_trim =13.5)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
}

data.od2 <- read.my.gc(file.od, plate.ind)

# Put the data together
data_fig4d <- bind_rows(data.od1, data.od2)

# Summarise to have only one data point for AUC for each well
data_fig4d_summary <- data_fig4d %>% ungroup() %>% 
  group_by(Well, Arabinose, TMP, Mutant) %>%
  summarise(auc = mean(auc_e))

# Summarize to show the average of the three replicates
data_fig4d_sum_curves <- data_fig4d %>% ungroup() %>%
  group_by(Arabinose, TMP, Mutant, time) %>%
  summarise(OD = mean(OD))

#### Show the percentage of growth recovery ####

# Pivot the area under the curve to calculate the differences
data_fig4d_wide <- data_fig4d_summary %>% ungroup() %>%
  group_by(Well, Arabinose, Mutant) %>%
  filter(!(is.na(Mutant))) %>%
  pivot_wider(names_from = TMP, values_from = auc, names_prefix = 'AUC_TMP_')

# Calculate the difference between the data with and without TMP
# and then the percentage of growth recovery
data_fig4d_wide %<>% mutate(diff_auc = AUC_TMP_10 - AUC_TMP_0)

# Add a column for the maximum difference (0% arabinose)
max_diff <- data_fig4d_wide %>% filter(Arabinose == 0) %>% ungroup() %>%
  group_by(Arabinose, Mutant) %>%
  summarise(max_diff = median(diff_auc))

# Use a join to add the maximum difference
data_fig4d_final <- left_join(x = data_fig4d_wide, 
                              y = max_diff %>% ungroup() %>% select(-Arabinose), 
                              by = c('Mutant' = 'Mutant'))

data_fig4d_final %<>% 
  mutate(pct_recovery = 100 * (AUC_TMP_10 / AUC_TMP_0))

comps <- compare_means(pct_recovery~Mutant, data = data_fig4d_final %>% ungroup(), 
                       paired = FALSE, group.by = 'Arabinose', method = 't.test') %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
                           ifelse(p >= 0.01, str_c('p = ', round(p, 2), sep = ''), 
                                  sprintf("p = %2.1e", as.numeric(p)))
  ),
  y_pos = c(110)
  )

comps %<>% mutate(Arabinose = as.factor(Arabinose))

data_fig4d_final_summ <- data_fig4d_final %>% ungroup() %>%
  group_by(Arabinose, Mutant) %>%
  summarise(mean_recovery = mean(pct_recovery), 
            sem_recovery = sd(pct_recovery) / sqrt(n()), 
            num_samples = n())


## A function to highlight some of the axis ticks
highlight = function(x, pat, color="black", family="") {
  ifelse(x %in% pat, glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

p_fig4d <- data_fig4d_final %>% 
  ggplot(aes(x = as.factor(Arabinose), y = pct_recovery, shape = Mutant)) +
  geom_point(data = data_fig4d_final_summ, inherit.aes = F, size = 7,
             aes(x = as.factor(Arabinose), shape = Mutant, y = mean_recovery)) +
  scale_shape_manual(values = c(18, 19)) + # 18, 19
  geom_line(aes(group = Mutant, x = as.factor(Arabinose), y = mean_recovery),
            data = data_fig4d_final_summ, inherit.aes = F, size = 0.5) +
  theme(axis.title = element_text(size = 20, face = 'bold'), 
        axis.text = element_text(size = 18), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', 
        legend.justification = 'center', 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18), 
        strip.text = element_text(size = 20, face = 'bold'), 
        strip.background = element_rect(fill = 'white')) +
  xlab('Arabinose (% m/v)') + ylab('Recovered growth (%)') +
  labs(shape = '') +
  geom_text(data = comps, inherit.aes = F,
            aes(x = Arabinose, label = p.format), y = 110, size = 5) +
  scale_y_continuous(breaks = seq(from = 0, to = 120, by =  20), limits = c(0, 120)) +
  scale_x_discrete(labels= function(x) highlight(x, c(0.01, 0.025, 0.05, 0.2, 0.4), "black")) +
  theme(axis.text.x=element_markdown())
p_fig4d


p_fig4_top <- plot_grid(p_fig4a, p_fig4b, labels = c('A', 'B'), nrow = 1,
                       label_size = 20, label_fontface = 'bold')

p_fig4 <- plot_grid(p_fig4_top, p_fig4c, p_fig4d, labels = c('', 'C', 'D'), nrow = 3, 
                   label_size = 20, label_fontface = 'bold', 
                   rel_heights = c(1, 1.5, 1.5))

ggsave(p_fig4, device = cairo_pdf, width = 13, height = 16, dpi = 300,
       filename = 'Figures/Main_figures/Fig4.pdf')

#### Figure 5: Heatmap of differences in fitness effects ####

data_fig_5 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose)

# Separate the data for ara 0.2
data_part_1 <- data_fig_5 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_5 %>%
  filter(Arabinose == 0.01)

## Subtract the scores
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_5_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

summary(data_fig_5_final$diffNormScore)

#### Block of code to save the data for the ML analysis and the chimerax figures ####
aa_changes_block <- all_data_complete %>% filter(TMP == 10, Arabinose == 0.01, Timepoint == 10)

aa_changes <- aa_changes_block[11:ncol(aa_changes_block)]

data_aa_diffNorm_0.2_0.01 <- left_join(x = data_fig_5_final, 
                                       y = aa_changes_block %>% 
                                         select(-mean_sel_coeff, -sd_sel_coeff, -sem_sel_coeff, -num_samples, 
                                                -Timepoint, -TMP), 
                                       by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                                              'Residue' = 'Residue', 
                                              'Arabinose_2' = 'Arabinose')) 

write.table(data_aa_diffNorm_0.2_0.01, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Data/Complete_datasets/dataset_diffNorm_ara0.2_ara0.01_index_differences.txt')

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
            file = 'Figures/Chimerax_figures/max_deltaS_ara0.2_ara0.01.txt')
write.table(deltaS_min, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/Chimerax_figures/min_deltaS_ara0.2_ara0.01.txt')
write.table(deltaS_mean, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
            file = 'Figures/Chimerax_figures/mean_deltaS_ara0.2_ara0.01.txt')

# Check deltaS for mutants with catalytic activity
data_fig_5_final %>% filter(WT_Residue == 'Q', Residue == 'C', Position == 67)
data_fig_5_final %>% filter(WT_Residue == 'I', Residue == 'L', Position == 68)
data_fig_5_final %>% filter(WT_Residue == 'I', Residue == 'M', Position == 68)
data_fig_5_final %>% filter(WT_Residue == 'Y', Residue == 'F', Position == 69)
data_fig_5_final %>% filter(WT_Residue == 'Y', Residue == 'H', Position == 69)


data_fig_5_final %>% filter(WT_Residue == 'S', Residue == 'A', Position == 65)


# Need to convert to wide formatted data
data_fig_5_final_df <- data_fig_5_final %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_5_final <- as.matrix(data_fig_5_final_df %>% select(-Position))

rownames(data_fig_5_final) <- data_fig_5_final_df$Position

# Get a matrix of true/false values for the synonymous codons
fig_5_bool <- data_fig_5 %>%
  filter(Arabinose == 0.2) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_5_bool_final <- as.matrix(fig_5_bool %>% select(-Position))

rownames(fig_5_bool_final) <- fig_5_bool$Position

# Need to reorder the columns in the matrices
data_fig_5_final <- data_fig_5_final[1:nrow(data_fig_5_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_5_bool_final <- fig_5_bool_final[1:nrow(fig_5_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

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
ha2 <- HeatmapAnnotation(
  `Dimerization interface` = data_interfaces_final$`A,C`,
  `Tetramerization interface` = data_interfaces_final$`A,D`,
  `DHF binding` = data_interfaces_final$DHF,
  `NADPH binding` = data_interfaces_final$NADPH,
  `Catalytic residues` = data_interfaces_final$Cat_residues,
  `Disordered region` = data_interfaces_final$Disordered_region,
  `Buried residues` = data_interfaces_final$Buried,
  show_annotation_name = T,
  annotation_name_gp = gpar(fontface = 'bold', fontsize = 16),
  annotation_name_side = 'left',
  show_legend = FALSE,
  col = list(
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


p_fig5a <- Heatmap(
  t(data_fig_5_final), cluster_columns = F, cluster_rows = F, 
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
    if (ara_0.01_bool_final[j,i]){
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-4, -2, 0, 2, 4) / 10,
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

p_fig5a

#### Figure 5B: Load the PDB structures with deltaS values ####
p_fig5b_min <- ggdraw() + 
  draw_image('Figures/Chimerax_figures/2rk1_min_deltaS_ara0.01_ara0.2.png')
p_fig5b_mean <- ggdraw() + 
  draw_image('Figures/Chimerax_figures/2rk1_mean_deltaS_ara0.01_ara0.2.png')
p_fig5b_max <- ggdraw() + 
  draw_image('Figures/Chimerax_figures/2rk1_max_deltaS_ara0.01_ara0.2.png')


#### Figure 5C: Check the distribution of delta(DMS scores) for critical sites ####

data_fig_5 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, Secondary_structure, rSASA)

# Separate the data for ara 0.2
data_part_1 <- data_fig_5 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_5 %>%
  filter(Arabinose == 0.01)

## Subtract the scores
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2",
                           "Arabinose_2", "Secondary_structure", "rSASA")

data_fig_5_final_new <- left_join(x = data_part_1, y = data_part_2, 
                                  by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue',
                                         'Residue' = 'Residue', 'rSASA' = 'rSASA',
                                         'Secondary_structure' = 'Secondary_structure')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)


data_fig_5_final_new <- left_join(x = data_part_1, y = data_part_2, 
                                  by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue', 
                                         'rSASA' = 'rSASA', 'Secondary_structure' = 'Secondary_structure')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Merge the dataframes
critical_sites_deltaDMS <- left_join(x = data_fig_5_final_new, y = data_interfaces_final, 
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
  mutate(Site = ifelse(Site == 'A,C', 'Dimerization interface' , 
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

## Do ANOVA and post hoc comparisons
m <- aov(diffNormScore ~ Site, data=critical_sites_deltaDMS)
anova_test <- anova(m)

tukey_test <- HSD.test(m, trt = 'Site')
tukey_test2 <- TukeyHSD(m)

## Save the table for the Tukey test for the supplementary tables
write.table(tukey_test2$Site, append = F, quote = F, sep = '\t', row.names = T, col.names = T, 
            file = 'Data/Tukey_test_fig5C.tsv')

tukey_annotation <- as.data.frame(tukey_test$groups) %>%
  mutate(y = c(0.4, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25))
tukey_annotation$Site_tukey <- rownames(tukey_annotation)


p_fig5c <- critical_sites_deltaDMS %>% rowwise() %>%
  mutate(mut_id = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  filter(Site != 'Only NADPH') %>%
  mutate(Site = factor(Site, levels = c('Disordered region', 'Catalytic residues',
                                        'DHF binding', 'NADPH binding',
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
  xlab('') +
  labs(y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                            bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                            sep = ''))) +
  ylim(-0.4, 0.4) +  
  annotate('text', x = 2.5, y = 0.35,
           label = expression(italic(s[weak] > s[opt])),
           parse = T, size = 10) +
  annotate('text', x = 2.5, y = -0.35,
           label = expression(italic(s[weak] < s[opt])),
           parse = T, size = 10) +
  geom_text(data = tukey_annotation %>% 
              filter(Site_tukey != 'Only NADPH') %>%
              mutate(Site_tukey = factor(Site_tukey, levels = c('Disordered region', 'Catalytic residues',
                                                                'DHF binding', 'NADPH binding',
                                                                'Buried residues',
                                                                'Tetramerization interface', 'Dimerization interface', 
                                                                'Unannotated'))),
            inherit.aes = F,
            aes(x = Site_tukey, y = y, label = groups), 
            size = 5)
p_fig5c


#### Figure 5D: deltaS vs s (optimal), protein sites ####

## Load data from validations to label them in figure 4C
# Read data for selected mutants
file.od <- 'Data/Growth_curves/14_07_21_bact.xlsx'
plate.ind <- 'Data/Growth_curves/Mutant_id_validation_DMS_growth_curves_IGA_14_07_21_no_empty.csv'

# function to process plates ----------------------------------------------
read.my.gc <- function(file, plate.index){
  pl <- read.xlsx(file,sheetIndex = 1, rowIndex = 4:99, stringsAsfators = FALSE,
                  header = F)
  ind <- read.csv2(plate.index, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
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
  gc_out <- SummarizeGrowthByPlate(d, t_trim =14.4)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
}

data.od1 <- read.my.gc(file.od, plate.ind)


all_data_complete_new <- all_data_complete %>% filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose) %>%
  rowwise()

# Load data from Dam 2000
dam2000_data <- read_delim('Data/Mutants_literature/mutants_Dam2000.csv', delim = '\t', locale = locale(decimal_mark = ','))

# Load the data from Strader 2001
strader2001_data <- read_delim('Data/Mutants_literature/mutants_Strader2001.csv', delim = '\t', locale = locale(decimal_mark = ','))

joined_sets <- inner_join(x = all_data_complete_new, 
                          y = strader2001_data,
                          by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue',
                                 'Residue' = 'Residue')) %>%
  mutate(Arabinose = as.factor(Arabinose))

joined_sets %<>% mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                                           ifelse(Arabinose == 0.025, 'Suboptimal', 
                                                  ifelse(Arabinose == 0.05, 'Near-optimal', 
                                                         ifelse(Arabinose == 0.2, 'Optimal', 
                                                                'Above-optimal')))))

# Mutants to annotate
list_mut <- joined_sets %>% select(Position, WT_Residue, Residue, `kcat/km`) %>% ungroup() %>%
  unique() %>% mutate(Genotype = str_c(WT_Residue, Position, Residue, sep = ''))

# Separate the data for ara 0.2
data_part_1 <- data_fig_5 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_5 %>%
  filter(Arabinose == 0.01)

## Subtract the scores
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2", 
                           "Secondary_structure", "rSASA")


# Join
data_fig_5_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue', 
                                     'Secondary_structure' = 'Secondary_structure', 'rSASA' = 'rSASA')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

mutants_highlight <- inner_join(x = data_fig_5_final %>% rowwise() %>%
                                  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                                                            ifelse(Residue == '*', 'Stop',
                                                                   'Amino acid substitution')),
                                         Genotype = str_c(WT_Residue, Position, Residue, sep = '')), 
                                y = list_mut, 
                                by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue',
                                       'Residue' = 'Residue', 'Genotype' = 'Genotype' 
                                ))

# Rename WT
list_mut %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))
mutants_highlight %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))

# Get data to represent the stop and WT codons
summary_wt_stop <- data_fig_5_final %>% rowwise() %>% ungroup() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop',
                                   'Amino acid substitution'))) %>%
  filter(mut_check != 'Amino acid substitution', Position >= 30, Position <= 70) %>%
  group_by(mut_check) %>%
  summarise(mean_s = mean(mean_sel_coeff), sem_s = sd(mean_sel_coeff) / sqrt(n()), 
            mean_ds = mean(diffNormScore), sem_ds = sd(diffNormScore) / sqrt(n()))

# Add the WT to the mutants to highlight
data_fig5d <- data_fig_5_final %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop',
                                   'Amino acid substitution'))) %>%
  filter(mut_check == 'Amino acid substitution') %>% 
  mutate(sem_s = 0, sem_ds = 0) %>% # These are single points
  select(mut_check, mean_sel_coeff, sem_s, diffNormScore, sem_ds)

colnames(data_fig5d) <- c('mut_check', 'mean_s','sem_s', 'mean_ds', 'sem_ds')

data_fig_5d <- rbind(data_fig5d, summary_wt_stop)

## Rename the column for mutants to highlight
mutants_highlight %<>% mutate(mean_s = mean_sel_coeff, mean_ds = diffNormScore)

mutants_highlight_final <- mutants_highlight %>% 
  select(mut_check, mean_s, mean_ds, `kcat/km`, Genotype) %>%
  add_row(summary_wt_stop %>% 
            select(mut_check, mean_s, mean_ds) %>%
            filter(mut_check == 'WT') %>%
            mutate(`kcat/km` = 0.433, 
                   Genotype = mut_check)
  )

annotation_x = 0.05

# Draw the figure
p_fig5d <- data_fig_5d %>% rowwise() %>%
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
  geom_label_repel(data = mutants_highlight_final,
                   aes(
                     label = Genotype), 
                   box.padding = 0.5, fontface = 'bold',
                   size = 7,
                   direction = 'both', 
  ) +
  scale_size_manual(values = c(3, 5, 5)) +
  guides(alpha = 'none') +
  scale_colour_manual(values = c('grey', 'black', 'blue')) +
  scale_shape_manual(values = c(16, 8, 16)) +
  scale_alpha_manual(values = c(0.5, 1, 1)) +
labs(colour = 'Mutation type', shape = 'Mutation type', size = 'Mutation type',
     y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                          bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                          sep = '')), 
     x =  expression(paste(bolditalic('s '), bold('(optimal promoter activity)'), 
                           sep = ''))
     
) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 20), 
    axis.title.y = element_text(face = 'bold', size = 28),
    axis.text = element_text(size = 18), 
    legend.position = 'top',
    legend.justification = 'center', 
    legend.box = 'top', 
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16), 
    legend.key.width = unit(1, 'cm')
    # legend.spacing.x = unit(1, 'cm')
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = expression(italic(s[weak] > s[opt])),
           parse = T, size = 10) +
  annotate('text', x = -0.7, y = -0.5, 
           label = expression(italic(s[weak] < s[opt])),
           parse = T, size = 10) +
  annotate('text', x = annotation_x, y = -0.1, parse = T, size = 6, hjust = 0, 
           label = expression(paste(italic(k[cat]), '/', italic(K[M])))) +
  annotate('text', x = annotation_x, y = -0.2, parse = F, size = 5, hjust = 0,
           label = 'WT = 0.433') +
  annotate('text', x = annotation_x, y = -0.3, parse = F, size = 5, hjust = 0,
           label = 'S65A = 0.275') +
  annotate('text', x = annotation_x, y = -0.4, parse = F, size = 5, hjust = 0,
           label = 'Q67C = 0.004') +
  annotate('text', x = annotation_x, y = -0.5, parse = F, size = 5, hjust = 0,
           label = 'I68L = 0.012') +
  annotate('text', x = annotation_x, y = -0.6, parse = F, size = 5, hjust = 0,
           label = 'I68M = 0.008') +
  annotate('text', x = annotation_x, y = -0.7, parse = F, size = 5, hjust = 0,
           label = 'Y69F = 0.038') +
  annotate('text', x = annotation_x, y = -0.8, parse = F, size = 5, hjust = 0,
           label = 'Y69H = 0.0001')

p_fig5d

## Figure 5E
data_regions <- read_delim('Data/data_annotation_2.txt', delim = '\t') %>%
  rowwise() %>%
  mutate(Unannotated = ifelse(sum(DHF, NADPH, `A,C`, `A,D`, Disordered_region, Cat_residues) == 0, 1, 0)) %>%
  pivot_longer(cols = c('DHF', 'NADPH', 'A,C', 'A,D', 'Disordered_region', 'Cat_residues', 'Buried', 'Unannotated'), 
               names_to = 'Region', values_to = 'Value') %>%
  filter(Value == 1)

data_regions_new <- left_join(x = data_fig_5_final, 
                              y = data_regions %>% select(Position, Region),
                              by = c('Position' = 'Position')) %>%
  unique()

# Add regions to highlighted mutants
mutants_highlight_regions <- left_join(x = mutants_highlight,
                                       y = data_regions_new %>% select(Position, WT_Residue, Residue, 
                                                                       Region), 
                                       by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                                              'Residue' = 'Residue'))  %>%
  mutate(Region = ifelse(Region == 'A,C', 'Dimerization interface', 
                         ifelse(Region == 'A,D', 'Tetramerization interface', 
                                ifelse(Region == 'Disordered_region', 'Disordered region',
                                       ifelse(Region == 'Cat_residues', 'Catalytic residues',
                                              ifelse(Region == 'Buried', 'Buried residues', 
                                                     Region))))))


data_regions <- read_delim('Data/data_annotation_2.txt', delim = '\t') %>%
  rowwise() %>%
  mutate(Unannotated = ifelse(sum(DHF, NADPH, `A,C`, `A,D`, Disordered_region, Cat_residues) == 0, 1, 0)) %>%
  pivot_longer(cols = c('DHF', 'NADPH', 'A,C', 'A,D', 'Disordered_region', 'Cat_residues', 'Buried', 'Unannotated'), 
               names_to = 'Region', values_to = 'Value') %>%
  filter(Value == 1)

data_regions_new <- left_join(x = data_fig_5_final, 
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
  mutate(Region = ifelse(Region == 'A,C', 'Dimerization interface', 
                         ifelse(Region == 'A,D', 'Tetramerization interface', 
                                ifelse(Region == 'Disordered_region', 'Disordered region',
                                       ifelse(Region == 'Cat_residues', 'DHF binding / catalytic',
                                              ifelse(Region == 'Buried', 'Buried residues', 
                                                     ifelse(Region == 'NADPH', 'NADPH binding',
                                                            Region)))))))

data_regions_new %<>% rowwise() %>%
  mutate(Region = ifelse(Region == 'A,C', 'Dimerization interface', 
                         ifelse(Region == 'A,D', 'Tetramerization interface',
                                ifelse(Region == 'Disordered_region', 'Disordered region',
                                       ifelse(Region == 'Cat_residues', 'DHF binding / catalytic',
                                              ifelse(Region == 'NADPH', 'NADPH binding',
                                                     ifelse(Region == 'Buried', 'Buried residues',
                                                            Region))))))) %>% 
  filter(Region != 'Unannotated', Region != 'DHF') %>%
  mutate(Region_stop = ifelse(Residue == '*', 'Stop', Region)) %>%
  mutate(Region = factor(Region, levels = c('Disordered region',
                                            'DHF binding / catalytic', 'NADPH binding',
                                            'Buried residues',
                                            'Tetramerization interface', 'Dimerization interface'
  ))) %>%
  mutate(Region_stop = factor(Region_stop, levels = c('Disordered region',
                                                      'DHF binding / catalytic', 'NADPH binding',
                                                      'Buried residues',
                                                      'Tetramerization interface', 'Dimerization interface', 
                                                      'Stop'))) %>%
  arrange(Region_stop)

mutants_highlight_regions %<>% 
  filter(Region != 'Unannotated', Genotype != 'WT') %>%
  mutate(Region = factor(Region, levels = c('Disordered region',
                                            'DHF binding / catalytic', 'NADPH binding',
                                            'Buried residues',
                                            'Tetramerization interface', 'Dimerization interface' 
                                            
  )))


p_fig5e <- data_regions_new %>%
  ggplot(aes(x = mean_sel_coeff, y = diffNormScore)) +
  geom_point( 
    aes(colour = Region_stop, size = Region_stop, 
        shape = Region_stop), alpha = 0.7) +
  facet_wrap(~Region, nrow = 3, scales = 'free') +
  labs(
    y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                         bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                         sep = '')), 
    x =  expression(paste(bolditalic('s '), bold('(optimal promoter activity)'), 
                          sep = ''))
  ) +
  scale_colour_manual(values = c('#CC79A7',
                                 '#E69F00', '#D55E00',
                                 '#9933ff',
                                 '#0072B2', '#009E73',
                                 'black'
  )) +
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
p_fig5e

# Put the figures together
p_fig5b <- plot_grid(p_fig5b_min, p_fig5b_mean, p_fig5b_max, ncol = 3, 
                     labels = c('B', '', ''), label_size = 40, label_fontface = 'bold')

p_fig5de <- plot_grid(p_fig5d, p_fig5e, ncol = 2, labels = c('D', 'E'), label_fontface = 'bold', 
                      label_size = 40)

#### Join all the panels of figure 3 ####
p_fig5_new <- plot_grid(grid.grabExpr(draw(p_fig5a)), 
                        p_fig5c,
                        p_fig5b,
                        p_fig5de,
                        nrow = 4, labels = c('A', '', 'C'), 
                        rel_heights = c(1, 1, 1, 1), label_size = 40, 
                        label_fontface = 'bold')

ggsave(p_fig5_new, device = cairo_pdf, width = 17, height = 30, dpi = 300,
       filename = 'Figures/Main_figures/Fig5_buried_2022-06-03_noSecStruc.pdf')

#### Figure 6: Slightly destabilizing mutants are masked by higher expression ####

#### Figure 6B ####
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
                                                         ifelse(Arabinose == 0.2, 'Optimal', 
                                                                'Above-optimal')))))

# Mutants to annotate
list_mut <- joined_sets %>% select(Position, WT_Residue, Residue) %>% ungroup() %>%
  unique() %>% mutate(Genotype = str_c(WT_Residue, Position, Residue, sep = ''))

data_fig_6b <- all_data_complete_new %>% ungroup() %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D)

# Separate the data for ara 0.2
data_part_1 <- data_fig_6b %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_6b %>%
  filter(Arabinose == 0.01)

## Subtract the scores
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2", 
                           "Mean_ddG_stab_HET", "Mean_ddG_int_HM_A_C", "Mean_ddG_int_HM_A_D")

# Join
data_fig_6b_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue',
                                     'Residue' = 'Residue', 'Mean_ddG_stab_HET' = 'Mean_ddG_stab_HET',
                                     'Mean_ddG_int_HM_A_C' = 'Mean_ddG_int_HM_A_C', 
                                     'Mean_ddG_int_HM_A_D' = 'Mean_ddG_int_HM_A_D')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

mutants_highlight <- data_fig_6b_final %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense')), 
         Genotype = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  filter(Genotype %in% list_mut$Genotype)

# Rename WT
list_mut %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))
mutants_highlight %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))

# Get data to represent the stop and WT codons
summary_wt_stop <- data_fig_6b_final %>% rowwise() %>% ungroup() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense'))) %>% 
  filter(mut_check != 'Missense') %>%
  group_by(mut_check) %>%
  summarise(mean_s = mean(mean_sel_coeff), sem_s = sd(mean_sel_coeff) / sqrt(n()), 
            mean_ds = mean(diffNormScore), sem_ds = sd(diffNormScore) / sqrt(n()),
            ddG_stab = mean(Mean_ddG_stab_HET), 
            ddG_dim_int = mean(Mean_ddG_int_HM_A_C), 
            ddG_tet_int = mean(Mean_ddG_int_HM_A_D)) %>%
  mutate(Region = 'Other')

data_fig6b <- data_fig_6b_final %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
  filter(mut_check %in% c('Missense', 'WT')) %>% 
  mutate(sem_s = 0, sem_ds = 0, 
         Region = ifelse(Position <= 20, 'IDR', 'Other')) %>% # These are single points
  select(mut_check, mean_sel_coeff, sem_s, diffNormScore, sem_ds, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D, Region) %>%
  ungroup()

colnames(data_fig6b) <- c('mut_check', 'mean_s','sem_s', 'mean_ds', 'sem_ds',
  'ddG_stab', 'ddG_dim_int', 'ddG_tet_int', 'Region')

data_fig6b <- rbind(data_fig6b, summary_wt_stop)

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
data_fig6b %<>%
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
  mutate(bins_stab = ifelse(bins_stab == '[-40,0]', '<= 0', 
                            ifelse(bins_stab == '(5,100]', '> 5', bins_stab)), 
         bins_dim_int = ifelse(bins_dim_int == '[-40,0]', '<= 0',
                               ifelse(bins_dim_int == '(5,100]', '> 5', bins_dim_int)), 
         bins_tet_int = ifelse(bins_tet_int == '[-40,0]', '<= 0',
                               ifelse(bins_tet_int == '(5,100]', '> 5', bins_tet_int))) %>%
  mutate(bins_stab = ifelse(Region == 'IDR', 'IDR', bins_stab), 
         bins_dim_int = ifelse(Region == 'IDR', 'IDR', bins_dim_int), 
         bins_tet_int = ifelse(Region == 'IDR', 'IDR', bins_tet_int)) %>%
  mutate(bins_stab = factor(bins_stab, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5')),
         bins_dim_int = factor(bins_dim_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5')), 
         bins_tet_int = factor(bins_tet_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5'))
         )

## Get the counts of mutations
table(data_fig6b$bins_stab)
table(data_fig6b$bins_dim_int)
table(data_fig6b$bins_tet_int)

# Check correlation between deltaS and destabilizing effects
corr_ddg_stab_deltaS <- data_fig6b %>% filter(!(is.na(ddG_stab)), mut_check != 'WT') %>%
  select(ddG_stab, ddG_dim_int, ddG_tet_int, mean_ds)

cor.test(corr_ddg_stab_deltaS$ddG_stab, corr_ddg_stab_deltaS$mean_ds, method = 'spearman')
cor.test(corr_ddg_stab_deltaS$ddG_dim_int, corr_ddg_stab_deltaS$mean_ds, method = 'spearman')
cor.test(corr_ddg_stab_deltaS$ddG_tet_int, corr_ddg_stab_deltaS$mean_ds, method = 'spearman')

cor(corr_ddg_stab_deltaS, method = 'spearman')

# Draw the figure for ddG stability
p_fig6b_stab <- data_fig6b %>% rowwise() %>%
  filter(!(is.na(bins_stab))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_stab), size = 3) +
  guides(size = 'none', alpha = 'none') +
  scale_colour_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  labs(colour = expression(
    paste(bold('\u0394\u0394'), bolditalic('G'), bold(' subunit stability (kcal/mol)'), sep = '')), 
    x = expression(paste(bolditalic('s'), bold(' (optimal promoter activity)'), sep = '')), 
    y = expression(paste(bold('\u0394'), bolditalic(s[weak]), ' (', 
      bolditalic(s[weak]), bold(' - '), bolditalic(s[opt]),
      ')', sep = ''))
  ) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 24), 
    axis.title.y = element_text(face = 'bold', size = 30),
    axis.text = element_text(size = 22), 
    legend.position = 'top',
    legend.justification = 'center', 
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22), 
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = expression(paste(italic(s[weak]), ' > ', italic(s[opt]), sep = '')),
           parse = T, size = 10) +
  annotate('text', x = -0.7, y = -0.5, 
           label = expression(paste(italic(s[weak]), ' < ', italic(s[opt]), sep = '')),
           parse = T, size = 10) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5, 
                               nrow = 1))
p_fig6b_stab

# Draw the figure for ddG dim int
p_fig6b_dim_int <- data_fig6b %>% rowwise() %>%
  filter(!(is.na(bins_dim_int))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_dim_int), size = 3) +
  guides(size = 'none', alpha = 'none') +
  scale_colour_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  labs(colour = expression(
    paste(bold('\u0394\u0394'), bolditalic('G'), bold(' dimer binding (kcal/mol)'), sep = '')), 
    x = expression(paste(bolditalic('s'), bold(' (optimal promoter activity)'), sep = '')), 
    y = expression(paste(bold('\u0394'), bolditalic(s[weak]), ' (', 
                         bolditalic(s[weak]), bold(' - '), bolditalic(s[opt]),
                         ')', sep = '')))+
  theme(
    axis.title.x = element_text(face = 'bold', size = 24), 
    axis.title.y = element_text(face = 'bold', size = 30),
    axis.text = element_text(size = 22), 
    legend.position = 'top',
    legend.justification = 'center', 
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           label = expression(paste(italic(s[weak]), ' > ', italic(s[opt]), sep = '')),
           parse = T, size = 10) +
  annotate('text', x = -0.7, y = -0.5, 
           label = expression(paste(italic(s[weak]), ' < ', italic(s[opt]), sep = '')),
           parse = T, size = 10) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5,
                               nrow = 1))

p_fig6b_dim_int

# Draw the figure for ddG tet int
p_fig6b_tet_int <- data_fig6b %>% rowwise() %>%
  filter(!(is.na(bins_tet_int))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_tet_int), size =3) +
  # geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
  #                  box.padding = 0.4, size = 6, fontface = 'bold') +
  guides(size = 'none', alpha = 'none') +
  scale_colour_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  # labs(colour = '\u0394\u0394G tet. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
  labs(colour = expression(
    paste(bold('\u0394\u0394'), bolditalic('G'), bold(' tetramer binding (kcal/mol)'), sep = '')), 
    # x = expression(paste(bolditalic('s'), bold(' (optimal expression)'), sep = '')), 
    x = expression(paste(bolditalic('s'), bold(' (optimal promoter activity)'), sep = '')), 
    y = expression(paste(bold('\u0394'), bolditalic(s[weak]), ' (', 
                         bolditalic(s[weak]), bold(' - '), bolditalic(s[opt]),
                         ')', sep = ''))
  ) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 24), 
    axis.title.y = element_text(face = 'bold', size = 30),
    axis.text = element_text(size = 22), 
    legend.position = 'top',
    legend.justification = 'center', 
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22)
  ) + 
  xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
  annotate('text', x = -0.7, y = 0.25, 
           # label = 's[weak] > s[opt]',
           label = expression(paste(italic(s[weak]), ' > ', italic(s[opt]), sep = '')),
           parse = T, size = 10) +
  annotate('text', x = -0.7, y = -0.5, 
           # label = 's[weak] < s[opt]',
           label = expression(paste(italic(s[weak]), ' < ', italic(s[opt]), sep = '')),
           parse = T, size = 10) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5, 
                               nrow = 1))
p_fig6b_tet_int

#### Figure 6C: Destabilizing mutants ####

ddg_sel_coeff_data <- all_data_complete %>% 
  select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP, mean_sel_coeff, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D) %>%
  filter(# !(is.na(Mean_ddG_stab_HET)),
         Timepoint == 10, TMP == 10)

summary(ddg_sel_coeff_data$Mean_ddG_int_HM_A_C)
summary(ddg_sel_coeff_data$Mean_ddG_int_HM_A_D)
summary(ddg_sel_coeff_data$Mean_ddG_stab_HET)

# Use the same meaningful bins as above
ddg_sel_coeff_data %<>% 
  filter(Residue != '*') %>%
  mutate(bins_stab = cut(Mean_ddG_stab_HET, breaks = c(-40, 0, 1, 2, 5, 100),
                         include.lowest = TRUE, na.rm = TRUE), 
         bins_dim_int = cut(Mean_ddG_int_HM_A_C, breaks = c(-40, 0, 1, 2, 5, 100),
                            include.lowest = TRUE, na.rm = TRUE),
         bins_tet_int = cut(Mean_ddG_int_HM_A_D, breaks = c(-40, 0, 1, 2, 5, 100),
                            include.lowest = TRUE, na.rm = TRUE)
  ) %>% mutate(bins_stab = as.character(bins_stab), 
               bins_dim_int = as.character(bins_dim_int), 
               bins_tet_int = as.character(bins_tet_int)) %>%
  mutate(bins_stab = ifelse(Position <= 20, 'IDR', bins_stab), 
         bins_dim_int = ifelse(Position <= 20, 'IDR', bins_dim_int), 
         bins_tet_int = ifelse(Position <= 20, 'IDR',  bins_tet_int)) %>%
  # Rename bins
  mutate(bins_stab = ifelse(bins_stab == '[-40,0]', '<= 0', 
                            ifelse(bins_stab == '(5,100]', '> 5', bins_stab)), 
         bins_dim_int = ifelse(bins_dim_int == '[-40,0]', '<= 0',
                               ifelse(bins_dim_int == '(5,100]', '> 5', bins_dim_int)), 
         bins_tet_int = ifelse(bins_tet_int == '[-40,0]', '<= 0',
                               ifelse(bins_tet_int == '(5,100]', '> 5', bins_tet_int))) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Above-optimal', 
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                                  'Optimal', 'Above-optimal')))

# Add nonsense mutants for reference
nonsense_mutants <- all_data_complete %>% 
  filter(TMP == 10, Residue == '*', Position >= 30, Position <= 70) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Above-optimal', 
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
    bins_stab = factor(bins_stab, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')),
    bins_dim_int = factor(bins_dim_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    bins_tet_int = factor(bins_tet_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                             'Optimal', 'Above-optimal'))
  )
  
mut_counts_stab <-table(ddg_sel_coeff_data$bins_stab) / 5
mut_counts_stab <- as.data.frame(t(mut_counts_stab))

mut_counts_dim_int <-table(ddg_sel_coeff_data$bins_dim_int) / 5
mut_counts_dim_int <- as.data.frame(t(mut_counts_dim_int))

mut_counts_tet_int <-table(ddg_sel_coeff_data$bins_tet_int) / 5
mut_counts_tet_int <- as.data.frame(t(mut_counts_tet_int))

#### Look for mutations that only affect one of the three ddGs ####

## Prepare the data frame
ddg_sel_coeff_data <- all_data_complete %>%
  select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP, mean_sel_coeff,
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D) %>%
  filter(# !(is.na(Mean_ddG_stab_HET)),
         Timepoint == 10, TMP == 10)

# Use the same meaningful bins
ddg_sel_coeff_data %<>%
  filter(Residue != '*') %>%
  mutate(bins_stab = cut(Mean_ddG_stab_HET, breaks = c(-40, 0, 1, 2, 5, 100),
                         include.lowest = TRUE, na.rm = TRUE),
         bins_dim_int = cut(Mean_ddG_int_HM_A_C, breaks = c(-40, 0, 1, 2, 5, 100),
                            include.lowest = TRUE, na.rm = TRUE),
         bins_tet_int = cut(Mean_ddG_int_HM_A_D, breaks = c(-40, 0, 1, 2, 5, 100),
                            include.lowest = TRUE, na.rm = TRUE)
  ) %>% mutate(bins_stab = as.character(bins_stab),
               bins_dim_int = as.character(bins_dim_int),
               bins_tet_int = as.character(bins_tet_int)) %>%
  mutate(bins_stab = ifelse(Position <= 20, 'IDR', bins_stab), 
         bins_dim_int = ifelse(Position <= 20, 'IDR', bins_dim_int), 
         bins_tet_int = ifelse(Position <= 20, 'IDR', bins_tet_int)) %>%
  # Rename bins
  mutate(bins_stab = ifelse(bins_stab == '[-40,0]', '<= 0',
                            ifelse(bins_stab == '(5,100]', '> 5', bins_stab)),
         bins_dim_int = ifelse(bins_dim_int == '[-40,0]', '<= 0',
                               ifelse(bins_dim_int == '(5,100]', '> 5', bins_dim_int)),
         bins_tet_int = ifelse(bins_tet_int == '[-40,0]', '<= 0',
                               ifelse(bins_tet_int == '(5,100]', '> 5', bins_tet_int))) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak',
                            ifelse(Arabinose == 0.025, 'Suboptimal',
                                   ifelse(Arabinose == 0.05, 'Near-optimal',
                                          ifelse(Arabinose == 0.2, 'Optimal',
                                                 ifelse(Arabinose == 0.4, 'Above-optimal',
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal',
                                                  'Optimal', 'Above-optimal')))

## Use the above data frame to separate mutations with small effects on two of the ddGs
mutations_stab <- ddg_sel_coeff_data %>% 
  filter(or(bins_stab == 'IDR', 
         and(abs(Mean_ddG_int_HM_A_C) < 0.5, abs(Mean_ddG_int_HM_A_D) < 0.5)))
summary(mutations_stab$Mean_ddG_stab_HET)

mutations_dim_int <- ddg_sel_coeff_data %>%
  filter(or(bins_dim_int == 'IDR', 
            and(abs(Mean_ddG_stab_HET) < 0.5, abs(Mean_ddG_int_HM_A_D) < 0.5)))
summary(mutations_dim_int$Mean_ddG_int_HM_A_C)

mutations_tet_int <- ddg_sel_coeff_data %>%
  filter(or(bins_tet_int == 'IDR',
            and(abs(Mean_ddG_stab_HET) < 0.5, abs(Mean_ddG_int_HM_A_C) < 0.5)))
summary(mutations_tet_int$Mean_ddG_int_HM_A_D)

# Add nonsense mutants for reference
nonsense_mutants <- all_data_complete %>% 
  filter(TMP == 10, Residue == '*', Position >= 30, Position <= 70) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Above-optimal', 
                                                        NA))))), 
         bins_stab = 'Stop', bins_dim_int = 'Stop', bins_tet_int = 'Stop') %>% 
  select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP, mean_sel_coeff, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D, 
         bins_stab, bins_dim_int, bins_tet_int, exp_level)

# Final data frame for the figure on ddG_stab bins
mutations_stab_final <- bind_rows(mutations_stab %>% ungroup() %>%
                                  select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C,
                                         -Mean_ddG_int_HM_A_D),
                                nonsense_mutants %>% ungroup() %>%
                                  select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C,
                                         -Mean_ddG_int_HM_A_D)) %>% 
  mutate(
    bins_stab = factor(bins_stab, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')),
    bins_dim_int = factor(bins_dim_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    bins_tet_int = factor(bins_tet_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                             'Optimal', 'Above-optimal'))
  )

# Final data frame for the figure on ddG_dim_int bins
mutations_dim_int_final <- bind_rows(mutations_dim_int %>% ungroup() %>%
                                    select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C,
                                           -Mean_ddG_int_HM_A_D),
                                  nonsense_mutants %>% ungroup() %>%
                                    select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C,
                                           -Mean_ddG_int_HM_A_D)) %>% 
  mutate(
    bins_stab = factor(bins_stab, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')),
    bins_dim_int = factor(bins_dim_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    bins_tet_int = factor(bins_tet_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                             'Optimal', 'Above-optimal'))
  )

# Final data frame for the figure on ddG_dim_int bins
mutations_tet_int_final <- bind_rows(mutations_tet_int %>% ungroup() %>%
                                       select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C,
                                              -Mean_ddG_int_HM_A_D),
                                     nonsense_mutants %>% ungroup() %>%
                                       select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C,
                                              -Mean_ddG_int_HM_A_D)) %>% 
  mutate(
    bins_stab = factor(bins_stab, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')),
    bins_dim_int = factor(bins_dim_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    bins_tet_int = factor(bins_tet_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                             'Optimal', 'Above-optimal'))
  )



### Draw the figures
## ddG_stab
# Check counts of mutations per bin (divide by 5 because each mutant appears 
# once for each arabinose concentration)
mut_counts_stab <- table(mutations_stab_final$bins_stab) / 5
mut_counts_stab <- as.data.frame(t(mut_counts_stab))


heatmap_height <- 14
heatmap_width <- 14
cell_text_size <- 18
header_text_cex <- 1.5

## Heatmap for mutations affecting only stability
data_heatmap_stab <- mutations_stab_final %>% ungroup() %>% 
  select(bins_stab, exp_level, mean_sel_coeff) %>%
  group_by(bins_stab, exp_level) %>%
  summarise(median_s = median(mean_sel_coeff))

data_heatmap_stab_wide <- data_heatmap_stab %>% ungroup() %>%
  pivot_wider(names_from = bins_stab, values_from = median_s)

# Change into a matrix and move the first column to rownames
needed_rownames <- data_heatmap_stab_wide$exp_level
data_heatmap_stab_final <- as.matrix(data_heatmap_stab_wide %>% select(-exp_level))
rownames(data_heatmap_stab_final) <- needed_rownames

## Get and format the counts of mutations
data_heatmap_stab_counts <- mutations_stab_final %>% ungroup() %>% 
  select(bins_stab, exp_level, mean_sel_coeff) %>%
  group_by(bins_stab, exp_level) %>%
  summarise(mut_count = n())

data_heatmap_stab_counts_wide <- data_heatmap_stab_counts %>% ungroup() %>%
  pivot_wider(names_from = bins_stab, values_from = mut_count)

# Change into a matrix and move the first column to rownames
needed_rownames <- data_heatmap_stab_counts_wide$exp_level
data_heatmap_stab_counts_final <- as.matrix(data_heatmap_stab_counts_wide %>% select(-exp_level))
rownames(data_heatmap_stab_counts_final) <- needed_rownames

## Draw the heatmap
# Draw the heatmap
header_text <- str_c(str_c('n = ', data_heatmap_stab_counts_final[1,]), sep = ',')

ha_new <- rowAnnotation(foo = anno_text(header_text, location = 0.5, just = 'center', 
                                           gp = gpar(border = 'white', cex = header_text_cex)))

cell_text_size = 22

p_heatmap_stab_final <- Heatmap(t(data_heatmap_stab_final), cluster_columns = F, cluster_rows = F, 
                          col = colorRamp2(
                            breaks = seq(-0.8, 0.2, length.out = 6),
                            colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#F4A582")),
                          show_column_names = T, row_names_side = 'left',
                          width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
                          border = T,
                          column_title = 'Promoter activity',
                          row_title_gp = gpar(fontsize=24, fontface = 'bold'),
                          row_names_rot = 0, 
                          row_names_centered = T,
                          row_names_gp = gpar(fontsize=18, fontface = 'bold'),
                          row_title =  expression(
                            paste(bold('\u0394\u0394'), bolditalic('G'),
                                  bold(' subunit stability (kcal/mol)'), sep = '')),
                          column_title_gp = gpar(fontsize=24, fontface = 'bold'),
                          column_title_side = 'bottom',
                          column_names_gp = gpar(fontsize=18,fontface='bold'),
                          cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.2f", t(data_heatmap_stab_final)[i, j]), 
                                      x, y,
                                      gp = gpar(fontsize = cell_text_size, fontface = 'bold'))
                          },
                          right_annotation = ha_new,
                          show_heatmap_legend = FALSE,
                          heatmap_legend_param = list(
                            at = round(seq(-0.8, 0.2, length.out = 6),2),
                            title = expression(paste(bold('Median '), bolditalic('s'), sep = '')), 
                            title_gp = gpar(fontsize = 20),
                            legend_height = unit(3.5, "cm"),
                            legend_width = unit(2, "cm"),
                            border='black',
                            lwd=1.7,
                            labels_gp = gpar(fontsize = 18),
                            title_position = "leftcenter-rot"
                          )
)
p_heatmap_stab_final


## Repeat for the bins on the dimer interface
data_heatmap_dim_int <- mutations_dim_int_final %>% ungroup() %>% 
  select(bins_dim_int, exp_level, mean_sel_coeff) %>%
  group_by(bins_dim_int, exp_level) %>%
  summarise(median_s = median(mean_sel_coeff))

data_heatmap_dim_int_wide <- data_heatmap_dim_int %>% ungroup() %>%
  pivot_wider(names_from = bins_dim_int, values_from = median_s)

# Change into a matrix and move the first column to rownames
needed_rownames <- data_heatmap_dim_int_wide$exp_level
data_heatmap_dim_int_final <- as.matrix(data_heatmap_dim_int_wide %>% select(-exp_level))
rownames(data_heatmap_dim_int_final) <- needed_rownames

## Get and format the counts of mutations
data_heatmap_dim_int_counts <- mutations_dim_int_final %>% ungroup() %>% 
  select(bins_dim_int, exp_level, mean_sel_coeff) %>%
  group_by(bins_dim_int, exp_level) %>%
  summarise(mut_count = n())

data_heatmap_dim_int_counts_wide <- data_heatmap_dim_int_counts %>% ungroup() %>%
  pivot_wider(names_from = bins_dim_int, values_from = mut_count)

# Change into a matrix and move the first column to rownames
needed_rownames <- data_heatmap_dim_int_counts_wide$exp_level
data_heatmap_dim_int_counts_final <- as.matrix(data_heatmap_dim_int_counts_wide %>% select(-exp_level))
rownames(data_heatmap_dim_int_counts_final) <- needed_rownames

# Draw the heatmap
header_text <- str_c(str_c('n = ', data_heatmap_dim_int_counts_final[1,]), sep = ',')

ha_new <- rowAnnotation(foo = anno_text(header_text, location = 0.5, just = 'center', 
                                           gp = gpar(border = 'white', cex = header_text_cex)))


## Draw the heatmap
p_heatmap_dim_int_final <- Heatmap(t(data_heatmap_dim_int_final), cluster_columns = F, cluster_rows = F,
                             col = colorRamp2(
                               breaks = seq(-0.8, 0.2, length.out = 6),
                               colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#F4A582")),
                             show_column_names = T, row_names_side = 'left',
                             width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
                             border = T,
                             column_title = 'Promoter activity',
                             row_title_gp = gpar(fontsize=24, fontface = 'bold'),
                             row_names_rot = 0, 
                             row_names_centered = T,
                             row_names_gp = gpar(fontsize=18, fontface = 'bold'),
                             row_title =  expression(
                               paste(bold('\u0394\u0394'), bolditalic('G'),
                                     bold(' dimer binding (kcal/mol)'), sep = '')),
                             column_title_gp = gpar(fontsize=24, fontface = 'bold'),
                             column_title_side = 'bottom',
                             column_names_gp = gpar(fontsize=18,fontface='bold'),
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.text(sprintf("%.2f", t(data_heatmap_dim_int_final)[i, j]), x, y,
                                         gp = gpar(fontsize = cell_text_size, fontface = 'bold'))
                             },
                             show_heatmap_legend = FALSE,
                             right_annotation = ha_new,
                             heatmap_legend_param = list(
                               at = round(seq(-0.8, 0.2, length.out = 6),2),
                               title = expression(paste(bold('Median '), bolditalic('s'), sep = '')), 
                               title_gp = gpar(fontsize = 20),
                               legend_height = unit(3.5, "cm"),
                               legend_width = unit(2, "cm"),
                               border='black',
                               lwd=1.7,
                               labels_gp = gpar(fontsize = 18),
                               title_position = "leftcenter-rot"
                             )
)
p_heatmap_dim_int_final

## Repeat for the bins on the tetramer interface
data_heatmap_tet_int <- mutations_tet_int_final %>% ungroup() %>% 
  select(bins_tet_int, exp_level, mean_sel_coeff) %>%
  group_by(bins_tet_int, exp_level) %>%
  summarise(median_s = median(mean_sel_coeff))

data_heatmap_tet_int_wide <- data_heatmap_tet_int %>% ungroup() %>%
  pivot_wider(names_from = bins_tet_int, values_from = median_s)

# Change into a matrix and move the first column to rownames
needed_rownames <- data_heatmap_tet_int_wide$exp_level
data_heatmap_tet_int_final <- as.matrix(data_heatmap_tet_int_wide %>% select(-exp_level))
rownames(data_heatmap_tet_int_final) <- needed_rownames

## Get and format the counts of mutations
data_heatmap_tet_int_counts <- mutations_tet_int_final %>% ungroup() %>% 
  select(bins_tet_int, exp_level, mean_sel_coeff) %>%
  group_by(bins_tet_int, exp_level) %>%
  summarise(mut_count = n())

data_heatmap_tet_int_counts_wide <- data_heatmap_tet_int_counts %>% ungroup() %>%
  pivot_wider(names_from = bins_tet_int, values_from = mut_count)

# Change into a matrix and move the first column to rownames
needed_rownames <- data_heatmap_tet_int_counts_wide$exp_level
data_heatmap_tet_int_counts_final <- as.matrix(data_heatmap_tet_int_counts_wide %>% select(-exp_level))
rownames(data_heatmap_tet_int_counts_final) <- needed_rownames

# Draw the heatmap
header_text <- str_c(str_c('n = ', data_heatmap_tet_int_counts_final[1,]), sep = ',')

ha_new <- rowAnnotation(foo = anno_text(header_text, location = 0.5, just = 'center', 
                                           gp = gpar(border = 'white', cex = header_text_cex)))

## Draw the heatmap
p_heatmap_tet_int_final <- Heatmap(t(data_heatmap_tet_int_final), cluster_columns = F, cluster_rows = F, 
                             col = colorRamp2(
                               breaks = seq(-0.8, 0.2, length.out = 6),
                               colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#F4A582")),
                             show_column_names = T, row_names_side = 'left',
                             width=unit(heatmap_width, 'cm'), height = unit(heatmap_height, 'cm'),
                             border = T,
                             column_title = 'Promoter activity',
                             row_title_gp = gpar(fontsize=24, fontface = 'bold'),
                             row_names_rot = 0, 
                             row_names_centered = T,
                             row_names_gp = gpar(fontsize=18, fontface = 'bold'),
                             row_title =  expression(
                               paste(bold('\u0394\u0394'), bolditalic('G'), 
                                     bold(' tetramer binding (kcal/mol)'), sep = '')),
                             column_title_gp = gpar(fontsize=24, fontface = 'bold'),
                             column_title_side = 'bottom',
                             column_names_gp = gpar(fontsize=18,fontface='bold'),
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.text(sprintf("%.2f", t(data_heatmap_tet_int_final)[i, j]), x, y,
                                         gp = gpar(fontsize = cell_text_size, fontface = 'bold'))
                             },
                             show_heatmap_legend = TRUE,
                             right_annotation = ha_new,
                             heatmap_legend_param = list(
                               at = round(seq(-0.8, 0.2, length.out = 6),2),
                               title = expression(paste(bold('Median '), bolditalic('s'), sep = '')), 
                               title_gp = gpar(fontsize = 20),
                               legend_height = unit(3.5, "cm"),
                               legend_width = unit(2, "cm"),
                               border='black',
                               lwd=1.7,
                               labels_gp = gpar(fontsize = 18),
                               title_position = "leftcenter-rot"
                             )
)
p_heatmap_tet_int_final


#### Figure 6D: Generate figures with pseudorecovery values ####

df_regions_tmp <- all_data_complete %>% 
  select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP, mean_sel_coeff, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D) %>%
  mutate(ddG_stab = Mean_ddG_stab_HET, 
         ddG_dim_int = Mean_ddG_int_HM_A_C, 
         ddG_tet_int = Mean_ddG_int_HM_A_D) %>%
  select(-Mean_ddG_stab_HET, -Mean_ddG_int_HM_A_C, -Mean_ddG_int_HM_A_D)

## Add the data for the protein sites
df_regions_final <- left_join(x = df_regions_tmp,
                              y = data_interfaces_final %>% 
                                  mutate('Dimer_interface' = `A,C`, 
                                         'Tetramer_interface' = `A,D`) %>%
                                  select(-`A,C`, -`A,D`), 
                                by = c('Position' = 'Position'))


df_regions_final_new <- left_join(x = df_regions_final %>% filter(Residue != '*'), 
                                  y = fig1C_data_sum %>% select(-sem_recovery, -num_samples), 
                                  by = c('Arabinose' = 'Arabinose')) 

df_regions_final_fig6d <- left_join(x = df_regions_final_new %>% filter(TMP == 10), 
                                    y = expression_level, 
                                    by = c('Arabinose' = 'Arabinose', 'exp_level' = 'exp_level')) %>%
  mutate(sel_coeff_factor = mean_sel_coeff + 1) %>%
  mutate(pseudofitness = sel_coeff_factor * mean_recovery)

# Bin the mutations by ddG
df_regions_final_fig6d %<>% ungroup() %>%
  filter(Residue != '*') %>% # Remove stop codons
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
  mutate(bins_stab = ifelse(bins_stab == '[-40,0]', '<= 0', 
                            ifelse(bins_stab == '(5,100]', '> 5', bins_stab)), 
         bins_dim_int = ifelse(bins_dim_int == '[-40,0]', '<= 0',
                               ifelse(bins_dim_int == '(5,100]', '> 5', bins_dim_int)), 
         bins_tet_int = ifelse(bins_tet_int == '[-40,0]', '<= 0',
                               ifelse(bins_tet_int == '(5,100]', '> 5', bins_tet_int))) %>%
  mutate(bins_stab = ifelse(Position <= 20, 'IDR', bins_stab), 
         bins_dim_int = ifelse(Position <= 20, 'IDR', bins_dim_int), 
         bins_tet_int = ifelse(Position <= 20, 'IDR', bins_tet_int)
         )%>%
  mutate(bins_stab = factor(bins_stab, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5')),
         bins_dim_int = factor(bins_dim_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5')), 
         bins_tet_int = factor(bins_tet_int, levels = c('IDR', '<= 0', '(0,1]', '(1,2]', '(2,5]', '> 5')),
         exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                                                'Optimal', 'Above-optimal'))
  )

## Isolate mutations that only seem to have an effect on subunit stability
df_regions_final_fig6d_stab <- df_regions_final_fig6d %>%
  filter(or(bins_stab == 'IDR', 
            and(abs(ddG_dim_int) < 0.5, abs(ddG_tet_int < 0.5))
            )
         ) %>%
  group_by(Arabinose, TMP, exp_level, mean_fluo, bins_stab) %>%
  summarise(median_pseudofitness = median(pseudofitness), 
            sem_pseudofitness = sd(pseudofitness) / sqrt(n()), 
            num_check = n())

p_fig6d_stab <- df_regions_final_fig6d_stab %>%
  ggplot(aes(x = exp_level, y = median_pseudofitness, colour = bins_stab, group = bins_stab,
             fill = bins_stab,
             ymax = median_pseudofitness + sem_pseudofitness, 
             ymin = median_pseudofitness - sem_pseudofitness)) +
  geom_point(size = 3) + geom_line() +
  geom_ribbon(alpha = 0.2) +
  scale_colour_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  scale_fill_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  xlab('Promoter activity') + ylab('Pseudo growth recovery (%)') + 
  theme(axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
        axis.title = element_text(size = 24, face = 'bold'), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top',
        legend.justification = 'center',
        legend.box.just = 'top',
        legend.title = element_text(size = 24, face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 22)) +
  labs(colour = expression(
    paste(bold('\u0394\u0394'), bolditalic('G'), bold(' subunit stability (kcal/mol)'), sep = '')), 
       fill = expression(
         paste(bold('\u0394\u0394'), bolditalic('G'), bold(' subunit stability (kcal/mol)'), sep = ''))
    ) +
  guides(colour = guide_legend(nrow = 1, title.position = 'top'))
p_fig6d_stab

## Isolate mutations that only seem to have an effect on the dimer interface
df_regions_final_fig6d_dim_int <- df_regions_final_fig6d %>%
  filter(or(bins_dim_int == 'IDR', 
            and(abs(ddG_stab) < 0.5, abs(ddG_tet_int < 0.5))
  )
  ) %>%
  group_by(Arabinose, TMP, exp_level, mean_fluo, bins_dim_int) %>%
  summarise(median_pseudofitness = median(pseudofitness), 
            sem_pseudofitness = sd(pseudofitness) / sqrt(n()), 
            num_check = n())

p_fig6d_dim_int <- df_regions_final_fig6d_dim_int %>%
  ggplot(aes(x = exp_level, y = median_pseudofitness, colour = bins_dim_int, group = bins_dim_int,
             fill = bins_dim_int,
             ymax = median_pseudofitness + sem_pseudofitness, 
             ymin = median_pseudofitness - sem_pseudofitness)) +
  geom_point(size = 3) + geom_line() +
  geom_ribbon(alpha = 0.2) +
  scale_colour_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  scale_fill_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  xlab('Promoter activity') + ylab('Pseudo growth recovery (%)') + 
  theme(axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
        axis.title = element_text(size = 24, face = 'bold'), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top',
        legend.justification = 'center',
        legend.title = element_text(size = 24, face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 22)) +
  labs(colour = expression(
    paste(bold('\u0394\u0394'), bolditalic('G'), bold(' dimer binding (kcal/mol)'), sep = '')), 
    fill = expression(
      paste(bold('\u0394\u0394'), bolditalic('G'), bold(' dimer binding (kcal/mol)'), sep = ''))
  ) +
  guides(colour = guide_legend(nrow = 1, title.position = 'top'))
p_fig6d_dim_int

## Isolate mutations that only seem to have an effect on the tetramer interface
df_regions_final_fig6d_tet_int <- df_regions_final_fig6d %>%
  filter(or(bins_tet_int == 'IDR', 
            and(abs(ddG_dim_int) < 0.5, abs(ddG_stab < 0.5))
  )
  ) %>%
  group_by(Arabinose, TMP, exp_level, mean_fluo, bins_tet_int) %>%
  summarise(median_pseudofitness = median(pseudofitness), 
            sem_pseudofitness = sd(pseudofitness) / sqrt(n()), 
            num_check = n())

p_fig6d_tet_int <- df_regions_final_fig6d_tet_int %>%
  ggplot(aes(x = exp_level, y = median_pseudofitness, colour = bins_tet_int, group = bins_tet_int,
             fill = bins_tet_int,
             ymax = median_pseudofitness + sem_pseudofitness, 
             ymin = median_pseudofitness - sem_pseudofitness)) +
  geom_point(size = 3) + geom_line() +
  geom_ribbon(alpha = 0.2) +
  scale_colour_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  scale_fill_manual(values = c('grey', '#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  xlab('Promoter activity') + ylab('Pseudo growth recovery (%)') + 
  theme(axis.text.y = element_text(size = 22),
        axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
        axis.title = element_text(size = 24, face = 'bold'), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top',
        legend.justification = 'center',
        legend.title = element_text(size = 24, face = 'bold', hjust = 0.5),
        legend.text = element_text(size = 22),) +
  labs(colour = expression(
    paste(bold('\u0394\u0394'), bolditalic('G'), bold(' tetramer binding (kcal/mol)'), sep = '')), 
    fill = expression(
      paste(bold('\u0394\u0394'), bolditalic('G'), bold(' tetramer binding (kcal/mol)'), sep = ''))
  ) +
  guides(colour = guide_legend(nrow = 1, title.position = 'top'))
p_fig6d_tet_int

## Put the panels together
p_fig6b <- plot_grid(p_fig6b_stab, p_fig6b_dim_int, p_fig6b_tet_int, ncol = 3)

ht_list = p_heatmap_stab_final + p_heatmap_dim_int_final + p_heatmap_tet_int_final
p_fig6c <- grid.grabExpr(
  draw(ht_list,
       ht_gap = unit(4, "cm"), ), 
)

p_fig6c <- plot_grid(
  grid.grabExpr(draw(p_heatmap_stab_final)), NULL,
  grid.grabExpr(draw(p_heatmap_dim_int_final)), NULL,
  grid.grabExpr(draw(p_heatmap_tet_int_final)),
  ncol = 5, align = 'hv', rel_widths = c(1, 0.05, 1, 0.05, 1)
)

p_fig6d <- plot_grid(p_fig6d_stab, p_fig6d_dim_int, p_fig6d_tet_int, ncol = 3)

p_fig6 <- plot_grid(p_fig6b,
                    p_fig6c,
                    p_fig6d, nrow = 3, labels = c('B', 'C', 'D'), 
                    label_size = 40, label_fontface = 'bold', 
                    rel_heights = c(1, 1.3, 1))
p_fig6
ggsave(p_fig6, device = cairo_pdf, width = 24, height = 24, dpi = 300,
       filename = 'Figures/Main_figures/Fig6.pdf')


ggsave(p_fig6, device = 'png', width = 24, height = 24, dpi = 300, 
       filename = 'Figures/Main_figures/Fig6.png')
