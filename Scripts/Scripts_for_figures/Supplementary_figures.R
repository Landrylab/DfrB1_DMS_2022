###########################################################
####        DMS_suppFigures_paper_2022-05-09           ####
#### This script will organize the data and prepare    ####
#### the supplementary figures for the paper.          ####
#### This is an updated version to work with the new   ####
#### data, including the data with 0.4% arabinose.     ####
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
library(FDRestimation)
library(ARTool)
library(lemon)
library(nlme)
library(lme4)
library(rlme)
library(shadowtext)
theme_set(theme_cowplot()  + 
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

# Load the complete dataset
all_data_complete <- read_delim(
  '../R67_DMS_February2022/Complete_datasets_MiSeq_NovaSeq/complete_dataset_avgBothSequencers.txt', delim = '\t')

#### Supp. figure 1B: Distributions of GFP cytometry data ####

# Load the main table with the sample IDs
# annotation_df <- read.xlsx('Data/Data_cytometry/Cytometry_file_names_2021_06_02_matched.xlsx',
#                            sheetIndex = 1, rowIndex = 1:43, header = T)

# annotation_df <- read_delim(
#   'Data/Data_cytometry/May_2022/Stats_run2_05_05_22_matched.csv',
#   locale = locale(decimal_mark = ','), delim = ',')

## Use a loop to load the rest of the data
# path_files <- 'Data/Data_cytometry/All_raw_files'
# path_files <- 'Data/Data_cytometry/May_2022/All_raw_files/'
# 
# list_files <- list.files(path_files)
# 
# all_data <- c()
# 
# for(infile in list_files){
#   # Load the data
#   new_data <- read_delim(file.path(path_files, infile), delim = ',', col_names = T) %>%
#     select(TIME, 'GRN-B-HLin', 'FSC-HLin', 'SSC-HLin', 'FSC-HLog', 'SSC-HLog')
#   
#   
#   # Add the name of the file as an ID (remove the csv)
#   new_data %<>% mutate(ID = substr(x = infile, start = 1, stop = (nchar(infile) - 4)))
#   
#   all_data <- bind_rows(all_data, new_data)
# }
# 
# ## Add the annotation data to identify each sample
# # all_data <- left_join(x = all_data, 
# #                       y = annotation_df,
# #                       by = c('ID' = 'File'))
# 
# # Add the annotation data to identify each sample
# all_data %<>% separate(col = ID, into = c('tmp', 'Well'), sep = 'am.')
# 
# all_data <- inner_join(x = all_data, 
#                        y = annotation_df,
#                        by = c('Well' = 'Well'))
# 
# 
# # Add the logarithms of the green fluorescence
# all_data %<>% mutate(GRNBHLog = ifelse(log10(`GRN-B-HLin`) < 0, 0, log10(`GRN-B-HLin`)),
#                      `FSC-HLog` = ifelse(log10(`FSC-HLin`) < 0, 0, log10(`FSC-HLin`)),
#                      `SSC-HLog` = ifelse(log10(`SSC-HLin`) < 0, 0, log10(`SSC-HLin`))
# )
# 
# all_data_processed <- all_data %>% rowwise() %>%
#   mutate(circle_test = (`FSC-HLog` - 1.25)^2 + (`SSC-HLog` - 1.25)^2) %>%
#   filter(circle_test < 1) %>%
#   mutate(Time = ifelse(Timepoint == 5, 'After 5 generations', NA)) %>%
#   mutate(
#     # Sample = str_c('Sample', toString(Sample), sep = ' '),
#     Arabinose = str_c(toString(Arabinose), '% arabinose', sep = ''),
#     # Replicate = ifelse(Sample <= 6, 1, 
#     #                    ifelse(Sample <= 24, floor((Sample - 1) / 6),
#     #                           ifelse(Sample > 24, floor((Sample -1) / 6) - 3, NA)
#     #                    )
#     # )
#   )
# 
# # all_data_processed %<>% arrange(Sample)
# 
# all_data_processed %<>% 
#   mutate(Arabinose = factor(Arabinose, levels = c('0% arabinose', 
#                                                   '0.01% arabinose', 
#                                                   '0.025% arabinose', 
#                                                   '0.05% arabinose',
#                                                   '0.2% arabinose', 
#                                                   '0.4% arabinose')
#   )
#   )
# 
# ## Not using
# # all_data_processed_summary <- all_data_processed %>% 
# #   ungroup() %>% group_by(Arabinose, Time, Well) %>%
# #   summarise(med_fluo = median(GRNBHLog)) %>%
# #   ungroup() %>% group_by(Arabinose, Time) %>%
# #   mutate(Replicate = row_number())
# 
# # Label replicates 1-3
# all_data_processed %<>% ungroup() %>%
#   separate(col = Well, into = c('Well_row', 'tmp'), sep = 1) %>%
#   mutate(Replicate = ifelse(Well_row == 'F', 1, 
#                             ifelse(Well_row == 'G', 2, 
#                                    ifelse(Well_row == 'H', 3, NA))))
# 
# p <- all_data_processed %>%
#   ggplot(aes(x = as.numeric(GRNBHLog))) +
#   # facet_grid(Arabinose ~ Time) +
#   facet_rep_grid(Arabinose ~ Time, scales = 'free', repeat.tick.labels = T) +
#   geom_density(aes(colour = as.factor(Replicate), group = as.factor(Replicate)),
#   size = 2.5, alpha = 0.3) +
#   xlab('GFP fluorescence') + ylab('Density')  +
#   # geom_histogram(aes(colour = as.factor(Replicate), group = as.factor(Replicate))) +
#   xlim(0, 4)  +
#   ylim(0, 4) +
#   theme(legend.position = 'none',
#         strip.background = element_rect(fill = 'white'), 
#         strip.text = element_text(face = 'bold', size = 24),
#         axis.text =  element_text(size = 24), axis.title = element_text(face = 'bold', size = 26),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank())
# p
# ggsave(p, device = cairo_pdf, width = 10, height = 35, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS1B_cytometry.pdf')
# ggsave(p, device = 'png', width = 10, height = 35, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS1B_cytometry.png')

# ## Use histograms to check the unexpected point
# p <- all_data_processed %>% rowwise() %>%
#   filter(Arabinose == '0.01% arabinose', Time == 'After 5 generations') %>%
#   mutate(Replicate = str_c('Replicate ', Replicate, sep = '')) %>%
#   ggplot(aes(x = as.numeric(GRNBHLog))) +
#   geom_histogram() +
#   facet_wrap(~Replicate, nrow = 3) +
#   xlab('GFP fluorescence') + ylab('Count') +
#   ggtitle('0.01% arabinose, time = after 5 generations')
# p
# ggsave(p, device = cairo_pdf, width = 10, height = 21, dpi = 300, 
#        filename = 'Figures/2022-04-07_supp_figures/FigS1B_ara0.01_t5_histogram.pdf')

#### Figure S1: Growth curves (WT with TMP vs WT without TMP) ####

plate.ind <- 'Data/GrowthCurves_February2022/Fitness_TMP_ID_IGA_23_02_22.xlsx'
file.od <- 'Data/GrowthCurves_February2022/Data_23_02_2022_bact.xlsx'

read.my.gc <- function(file, plate.index){
  pl <- read.xlsx(file,sheetIndex = 1, rowIndex = 5:100, stringsAsfators = FALSE,
                  header = F)
  # ind <- read.csv(plate.index, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  # ind <- read.csv2(plate.index, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  ind <- read.xlsx(plate.index, sheetIndex = 1, rowIndex = 1:97, header = T) 
  
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

  ## For the WT controls I will trim at t = 13.5
  gc_out <- SummarizeGrowthByPlate(d, t_trim =13.5)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
  
}


data.od1 <- read.my.gc(file.od, plate.ind)
## Remove wells on the border of the plate (rows A and H, columns 1 and 12)
data.od1 %<>% separate(Well, into = c('Well_row', 'Well_col'), sep = 1) %>%
 filter(!(Well_row %in% c('A', 'H')), !(Well_col %in% c(1, 12)))

# Subtract the blank for this experiment and multiply by 5 to make it OD / mL
# (experiment was carried out in 0.2 mL)
data.od1 %<>% mutate(OD = (OD - 0.088) * 5)

# Show the growth curves
p_figs1 <- data.od1 %>% ungroup() %>% rowwise() %>%
  filter(time <= 13.5) %>%
  mutate(TMP = str_c(TMP, ' \u00B5g/mL', sep = '')) %>%
  mutate(
    TMP = as.factor(TMP),
    Arabinose = str_c(Arabinose, '% arabinose')
    ) %>%
  filter(Mutant == 'WT') %>%
  mutate(Arabinose = factor(Arabinose, levels = c('0% arabinose', 
                                                  '0.01% arabinose',
                                                  '0.025% arabinose', 
                                                  '0.05% arabinose',
                                                  '0.2% arabinose',
                                                  '0.4% arabinose'))) %>%
  ggplot(aes(x = time, y = OD, colour = TMP, group = interaction(Mutant, Replicate, TMP))) + 
  facet_wrap(~Arabinose, ncol = 3, scales = 'free') +
  geom_line() +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(face = 'bold', size = 20),
        axis.text = element_text(size = 18),
        axis.line = element_line(),
        legend.position = 'top', 
        legend.title = element_text(size = 20), 
        strip.text = element_text(size = 20, face = 'bold'),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 18),
        legend.justification = 0.5) +
  guides(size = 'none', linetype = 'none', alpha = 'none') +
  xlab('Time (h)') + ylim(0, 6) +
  labs(y = expression(bold(OD[600])))
# p_figs1

# Edit the legend to make the lines stand out more
p_s1_legend <- get_legend(p_figs1 + geom_line(size = 3) +
                          theme(legend.position = 'top',
                                legend.title = element_text(size = 24),
                                legend.text = element_text(size = 22)))
p_figs1_final <- plot_grid(p_s1_legend, p_figs1 + theme(legend.position = 'none'), 
                           nrow = 2, rel_heights = c(0.1, 1))

p_figs1_final

ggsave(plot = p_figs1_final, device = cairo_pdf, width = 21, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS1_growthRecovery_DfrB1.pdf')
ggsave(plot = p_figs1_final, device = 'png', width = 21, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS1_growthRecovery_DfrB1.png')

#### Figure S2: Quality control of the DMS library ####

# Save the order for codons
codon_order = c('TAA', 'TAG', 'TGA', # *
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
)


# Load data
dms_qc <- read_delim('Data/MiSeqIsa/Results/clean_data_pBAD.csv', delim = '\t')
colnames(dms_qc)[1] <- 'Mutant_ID'
colnames(dms_qc)[4] <- 'Position'

total_by_pos <- dms_qc %>% ungroup() %>%
  group_by(Position) %>%
  summarise(read_total = sum(Read_count))

sum(total_by_pos$read_total)

# Use a join to add the total count per position to the matrix
dms_qc_new <- left_join(x = dms_qc, y = total_by_pos, by = c('Position' = 'Position'))

dms_qc <- dms_qc_new %>% mutate(Ratio_deep = Read_count * 100 / read_total) %>%
  mutate(Codon = factor(Codon, levels = codon_order)) 

# Pivot to accomodate in a table
data_heatmap_ratiodeep <- dms_qc %>% select(Position, Codon, Ratio_deep) %>%
  pivot_wider(id_cols = Position, names_from = Codon, values_from = Ratio_deep)

data_ratiodeep_matrix <- as.matrix(data_heatmap_ratiodeep %>% select(-Position))

rownames(data_ratiodeep_matrix) <- data_heatmap_ratiodeep$Position

# Rearrange according to the codon order
data_ratiodeep_matrix <- data_ratiodeep_matrix[, codon_order]

## Add encoded amino acids
aa_encoded <- dms_qc %>% select(Codon, AA) %>%
  unique()
aa_encoded %<>% arrange(factor(Codon, levels = codon_order))

### Prepare the WT matrix ###
fig_qc_bool <- dms_qc %>%
  mutate(WT_check = (WT_codon == Codon)) %>%
  select(Position, Codon, WT_check) %>%
  pivot_wider(names_from = Codon, values_from = WT_check)

fig_qc_bool_final <- as.matrix(fig_qc_bool %>% select(-Position))
rownames(fig_qc_bool_final) <- fig_qc_bool$Position

# Reorder according to codon order
fig_qc_bool_final <- fig_qc_bool_final[, codon_order]

# Modify the column names to indicate the encoded residue
for(i in 1:length(colnames(data_ratiodeep_matrix))){
  colnames(data_ratiodeep_matrix)[i] <- str_c(colnames(data_ratiodeep_matrix)[i], ' (',
                                              aa_encoded$AA[i], ')', sep = '')
  colnames(fig_qc_bool_final)[i] <- str_c(colnames(fig_qc_bool_final)[i], ' (',
                                          aa_encoded$AA[i], ')', sep = '')
}

column_labels = c('', '', '', '', '5', 
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
                  '', '', '', ''
)

# Draw the heatmap
p_figs2 <- Heatmap(t(data_ratiodeep_matrix), cluster_columns = F, cluster_rows = F, 
                ## Color ramp for the old normalized scores
                col = colorRamp2(
                  breaks = seq(0, 4, length.out = 7),
                  # colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#F4A582", "#D6604D")),
                  colors = viridis::viridis(7)),
                show_column_names = T, row_names_side = 'left',
                width=unit(31, 'cm'), height = unit(31, 'cm'),
                border = T,
                row_title = 'Codon (encoded residue)',
                row_title_gp = gpar(fontsize=22, fontface = 'bold'),
                row_names_rot = 0, 
                row_names_centered = T,
                row_names_gp = gpar(fontsize=16, fontface = 'bold'),
                column_title = 'Position', 
                column_title_side = 'bottom',
                column_names_gp = gpar(fontsize=16,fontface='bold'),
                column_title_gp = gpar(fontsize=22, fontface = 'bold'),
                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                  if (fig_qc_bool_final[j,i]){
                    # grid.text(expression(bold('s')), x, y)
                    grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
                  }      
                },
                # right_annotation = ha_empty,
                show_heatmap_legend = TRUE,
                heatmap_legend_param = list(
                  at = c(0, 2, 4),
                  title = "Mutant fraction in library (%)", 
                  title_gp = gpar(fontsize = 20),
                  legend_height = unit(3.5, "cm"),
                  legend_width = unit(2, "cm"),
                  border='black',
                  lwd=1.7,
                  labels_gp = gpar(fontsize = 18),
                  title_position = "leftcenter-rot"
                ), 
                column_labels = column_labels
)
p_figs2

ggsave(grid.grabExpr(draw(p_figs2)), width = 16, height = 14, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS2.heatmap_quality_control.pdf')
ggsave(grid.grabExpr(draw(p_figs2)), width = 16, height = 14, dpi = 300, device = 'png', 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS2.heatmap_quality_control.png')

#### Figure S3 was generated using BioRender ####

#### Figure S4: Heatmaps of correlation between replicates (TMP) ####

all_data_all_reps <- read_delim('../R67_DMS_February2022/Complete_datasets_MiSeq_NovaSeq/all_data_all_reps_bothSequencers.txt', 
                                delim = '\t')

all_data_all_reps_corr_TMP10 <- all_data_all_reps %>% rowwise() %>%
  mutate(Genotype = str_c(WT_Residue, Residue, Position)) %>%
  select(Genotype, Timepoint, TMP, Arabinose, sel_coeff, Sequencer, ID) %>%
  group_by(Genotype) %>%
  filter(TMP == 10, Timepoint == 10)

cor_matrix <- all_data_all_reps_corr_TMP10  %>% ungroup() %>%
  select(-Timepoint, -Arabinose, -TMP, -Sequencer) %>%
  pivot_wider(names_from = ID, values_from = sel_coeff) %>%
  select(-Genotype) %>% cor(method = 'spearman')

annotation_df <- all_data_all_reps_corr_TMP10 %>% ungroup() %>%
  select(-Genotype, -sel_coeff) %>%
  unique() %>%
  mutate(Arabinose = factor(Arabinose, levels = c(0.01, 0.025, 0.05, 0.2, 0.4))) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal',
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal',
                                                  'Optimal', 'Overexpressed')))

ha_corr <- HeatmapAnnotation(`Expression level` = annotation_df$exp_level,
                             show_annotation_name = F,
                             annotation_name_gp = gpar(fontface = 'bold', fontsize = 14),
                             annotation_name_side = 'left',
                             show_legend = TRUE,
                             annotation_legend_param = list(
                               `Expression level` = list(
                                 title_gp = gpar(fontsize = 22, fontface = 'bold'), 
                                 labels_gp = gpar(fontsize = 20)
                               )
                             ),
                             col = list(`Expression level` = c("Weak" = "#fed976", "Suboptimal" = "#fd8d3c", 
                                                               "Near-optimal" = "#bd0026", "Optimal" = "#80001a", 
                                                               'Overexpressed' = 'black')),
                             gp = gpar(col = "black")
)

ha_corr2 <- HeatmapAnnotation(`Expression level` = annotation_df$exp_level,
                              which = 'row',
                             show_annotation_name = F,
                             annotation_name_gp = gpar(fontface = 'bold', fontsize = 14),
                             # annotation_name_side = 'left',
                             show_legend = FALSE,
                             annotation_legend_param = list(
                               `Expression level` = list(
                                 title_gp = gpar(fontsize = 22, fontface = 'bold'), 
                                 labels_gp = gpar(fontsize = 20)
                               )
                             ),
                             col = list(`Expression level` = c("Weak" = "#fed976", "Suboptimal" = "#fd8d3c", 
                                                               "Near-optimal" = "#bd0026", "Optimal" = "#80001a", 
                                                               'Overexpressed' = 'black')),
                             gp = gpar(col = "black")
)


# Draw the heatmap
p_figs4 <- Heatmap(cor_matrix, cluster_columns = T, cluster_rows = T, 
             clustering_distance_rows = 'pearson',
             clustering_distance_columns = 'pearson',
             ## Color ramp for the old normalized scores
             col = colorRamp2(
               breaks = seq(0.5, 1, length.out = 7),
               colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#F4A582", "#D6604D")),
             show_column_names = T, row_names_side = 'left',
             width=unit(28, 'cm'), height = unit(28, 'cm'),
             border = T,
             row_title = 'Sample',
             row_title_gp = gpar(fontsize=22, fontface = 'bold'),
             row_names_rot = 0, 
             row_names_centered = T,
             row_names_gp = gpar(fontsize=18, fontface = 'bold'),
             column_title = 'Sample', 
             column_title_gp = gpar(fontsize=22, fontface = 'bold'),
             column_title_side = 'bottom',
             column_names_gp = gpar(fontsize=18,fontface='bold'),
             top_annotation = ha_corr,
             left_annotation = ha_corr2,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y,
                         gp = gpar(fontsize = 14, fontface = 'bold'))
             },
             show_heatmap_legend = TRUE,
             heatmap_legend_param = list(
               at = c(0.5, 0.75, 1),
               title = "Spearman ⍴", 
               title_gp = gpar(fontsize = 20),
               legend_height = unit(3.5, "cm"),
               legend_width = unit(2, "cm"),
               border='black',
               lwd=1.7,
               labels_gp = gpar(fontsize = 18),
               title_position = "leftcenter-rot"
             )
)
p_figs4
ggsave(grid.grabExpr(draw(p_figs4)), width = 18, height = 17, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS4.heatmap_corr_samples_TMP.pdf')
ggsave(grid.grabExpr(draw(p_figs4)), width = 18, height = 17, dpi = 300, device = 'png', 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS4.heatmap_corr_samples_TMP.png')

#### Figure S5: Ranks of Dam and Strader mutants ####

## Preprocess the DMS data
data_figs5 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose) %>%
  mutate(exp_level = ifelse(Arabinose == 0.2, 'Optimal', 
                            ifelse(Arabinose == 0.05, 'Near-optimal',
                                   ifelse(Arabinose == 0.025, 'Suboptimal', 
                                          ifelse(Arabinose == 0.01, 'Weak',
                                                 ifelse(Arabinose == 0.4, 'Overexpressed',
                                                        'Absent'))))))

# Load data from Dam 2000
dam2000_data <- read_delim('Data/mutants_Dam2000.csv', delim = '\t', locale = locale(decimal_mark = ','))


# Load the data from Strader 2001
strader2001_data <- read_delim('Data/mutants_Strader2001.csv', delim = '\t', locale = locale(decimal_mark = ','))

dam_strader_mutants <- bind_rows(dam2000_data %>% select(Position, WT_Residue, Residue), 
                                 strader2001_data %>% select(Position, WT_Residue, Residue) %>%
                                   filter(WT_Residue != 'WT'))

data_figs5_ranked <- data_figs5 %>% group_by(Arabinose) %>%
  mutate(fit_rank = rank(mean_sel_coeff))

data_figs5ab_ara0.01 <- left_join(x = data_figs5_ranked %>% filter(Arabinose == 0.01), 
                                  y = dam_strader_mutants %>% mutate(sel_mut = TRUE), 
                                  by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(sel_mut = ifelse(is.na(sel_mut), FALSE, TRUE))

labels_dam_strader_ara0.01 <- data_figs5ab_ara0.01 %>% filter(sel_mut, Arabinose == 0.01) %>% rowwise() %>%
  mutate(ID = str_c(WT_Residue, Position, Residue, sep = ''))

p_figs5a <- data_figs5ab_ara0.01 %>% 
  ggplot(aes(x = fit_rank, y = mean_sel_coeff)) + 
  geom_point(aes(colour = sel_mut, alpha = sel_mut)) +
  geom_label_repel(data = labels_dam_strader_ara0.01,
                   aes(label = ID), 
                   box.padding = 0.4,
                   direction = 'both', 
                   position = position_nudge_repel(x = 300, y = -0)
  ) +
  scale_colour_manual(values = c('black', 'red')) +
  scale_alpha_manual(values = c(0.2, 1)) +
  xlab('Rank') + # ylab('s (weak expression)') +
  labs(y = expression(bolditalic(s[weak]))) +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(face = 'bold', size = 20), 
        axis.title.y = element_text(face = 'bold', size = 22), 
        axis.text = element_text(size = 18), 
        legend.position = 'none'
  ) +
  ylim(-1, 0.4)
p_figs5a

data_figs5ab_ara0.2 <- left_join(x = data_figs5_ranked %>% filter(Arabinose == 0.2), 
                                 y = dam_strader_mutants %>% mutate(sel_mut = TRUE), 
                                 by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(sel_mut = ifelse(is.na(sel_mut), FALSE, TRUE))

labels_dam_strader <- data_figs5ab_ara0.2 %>% filter(sel_mut, Arabinose == 0.2) %>% rowwise() %>%
  mutate(ID = str_c(WT_Residue, Position, Residue, sep = ''))

p_figs5b <- data_figs5ab_ara0.2 %>% 
  ggplot(aes(x = fit_rank, y = mean_sel_coeff)) + 
  geom_point(aes(colour = sel_mut, alpha = sel_mut)) +
  geom_label_repel(
    data = labels_dam_strader, 
    aes(label = ID),
    box.padding = 0.4,
    direction = 'both', 
    position = position_nudge_repel(x = 300, y = -0)
  ) +
  scale_colour_manual(values = c('black', 'red')) +
  scale_alpha_manual(values = c(0.2, 1)) +
  xlab('Rank') + # ylab('s (optimal expression)') +
  labs(y = expression(bolditalic(s[opt]))) +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(face = 'bold', size = 20), 
        axis.title.y = element_text(face = 'bold', size = 22), 
        axis.text = element_text(size = 18), 
        legend.position = 'none'
  ) +
  ylim(-1, 0.4)
p_figs5b

p_figs5 <- plot_grid(p_figs5a, p_figs5b, nrow = 2, labels = c('A', 'B'), 
                     label_size = 20, label_fontface =  'bold')

ggsave(plot = p_figs5, device = cairo_pdf, width = 10, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS5_mutants_literature_ranks.pdf')
ggsave(plot = p_figs5, device = 'png', width = 10, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS5_mutants_literature_ranks.png')

#### Figure S6: Changes in fitness effects of deleterious mutants from the literature ####

# Use an inner join with our DMS data
data_fig_dam2000 <- inner_join(x = data_figs5 %>% ungroup(), 
                               y = dam2000_data,
                               by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue'))


# Use an inner join with our DMS data
data_fig_strader2001 <- inner_join(x = data_figs5 %>% ungroup(), 
                                   y = strader2001_data,
                                   by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue'))

# ## Draw the figures for panels C and D
# p_dam2000_newfig <- data_fig_dam2000 %>%rowwise() %>%
#   filter(exp_level %in% c('Weak', 'Optimal')) %>%
#   mutate(ID = str_c(WT_Residue, Position, Residue, sep = ''), 
#          exp_level = factor(exp_level, levels = c('Weak', 'Optimal'))
#   ) %>%
#   ggplot(aes(x = Pct_activity, y = mean_sel_coeff, label = ID)) + 
#   geom_point(aes(
#     colour = exp_level
#   )) +
#   geom_label_repel(aes(colour = exp_level), 
#                    fontface = 'bold', max.overlaps = Inf, show.legend = F, 
#                    box.padding = 0.4, direction = 'x', fill = 'black') +
#   xlab('% activity relative to WT [Dam, 2000]') + ylab('s') +
#   scale_colour_manual(values = c('#fed976', '#bd0026')) +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         axis.title = element_text(face = 'bold', size = 20), 
#         axis.text = element_text(size = 18), 
#         legend.position = 'top', 
#         legend.justification = 0.5, 
#         legend.title = element_text(size = 20), 
#         legend.text = element_text(size = 18)
#   ) + ylim(-1, 0) +
#   labs(colour = 'Expression level')
# p_dam2000_newfig
# 
# p_strader2001_newfig <- data_fig_strader2001 %>%  rowwise() %>%
#   filter(exp_level %in% c('Weak', 'Optimal')) %>%
#   mutate(ID = str_c(WT_Residue, Position, Residue, sep = ''), 
#          exp_level = factor(exp_level, levels = c('Weak', 'Optimal'))) %>%
#   ggplot(aes(x = `kcat/km`, y = mean_sel_coeff, label = ID)) + 
#   geom_point(aes(colour = exp_level)) +
#   geom_label_repel(aes(colour = exp_level), 
#                    fontface = 'bold', max.overlaps = Inf, show.legend = F, 
#                    box.padding = 0.4, direction = 'x', 
#                    fill = 'black') +
#   xlab('Catalytic efficiency (kcat / km) [Strader, 2001]') + ylab('s') +
#   scale_colour_manual(values = c('#fed976', '#bd0026')) +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         axis.title = element_text(face = 'bold', size = 20), 
#         axis.text = element_text(size = 18), 
#         legend.position = 'none'
#   ) + ylim(-1, 0) +
#   labs(colour = 'Expression level')
# p_strader2001_newfig

# figs6_legend <- get_legend(p_dam2000_newfig)
# figs6_plots <- plot_grid(p_dam2000_newfig + theme(legend.position = 'none'), 
#                          p_strader2001_newfig, nrow = 2, labels = c('A', 'B'), 
#                          label_size = 20, label_fontface = 'bold')
# 
# p_figs6 <- plot_grid(figs6_legend, figs6_plots, nrow = 2, 
#                      label_size = 20, label_fontface =  'bold', 
#                      rel_heights = c(0.1, 1))

## Try directly showing the deltaS, might be clearer ##
data_fig_dam2000_new <- data_fig_dam2000 %>% ungroup() %>% group_by(Position, WT_Residue, Residue) %>%
  select(-Arabinose) %>%
  filter(exp_level %in% c('Weak', 'Optimal')) %>%
  pivot_wider(names_from = exp_level, values_from = mean_sel_coeff) %>%
  mutate(diffNormScore = Weak - Optimal)

p_dam2000_newfig <- data_fig_dam2000_new %>% rowwise() %>%
  mutate(ID = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  ggplot(aes(x = Pct_activity, y = diffNormScore, label = ID)) + 
  geom_point() +
  geom_label_repel(fontface = 'bold', max.overlaps = Inf, show.legend = F, 
                   box.padding = 0.6, direction = 'both', size = 4) +
  xlab('% activity relative to WT [Dam, 2000]') + 
  labs(y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                            bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                            sep = ''))
       ) +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(face = 'bold', size = 20), 
        axis.title.y = element_text(face = 'bold', size = 22),
        axis.text = element_text(size = 18), 
        legend.position = 'top', 
        legend.justification = 0.5, 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)
  ) +
  ylim(-0.4, 0.4)
p_dam2000_newfig

data_fig_strader2001_new <- data_fig_strader2001 %>% ungroup() %>% group_by(Position, WT_Residue, Residue) %>%
  select(-Arabinose) %>%
  filter(exp_level %in% c('Weak', 'Optimal')) %>%
  pivot_wider(names_from = exp_level, values_from = mean_sel_coeff) %>%
  mutate(diffNormScore = Weak - Optimal)

p_strader2001_newfig <- data_fig_strader2001_new %>% rowwise() %>%
  mutate(ID = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  ggplot(aes(x = `kcat/km`, y = diffNormScore, label = ID)) + 
  geom_point() +
  geom_label_repel(fontface = 'bold', max.overlaps = Inf, show.legend = F, 
                   box.padding = 0.4, direction = 'both', size = 4) +
  # xlab('Catalytic efficiency (kcat / km) [Strader, 2001]') +
  # labs(y = expression(s[weak] - s[opt])) +
  labs(y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                            bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                            sep = '')), 
       x = expression(paste(bold('Catalytic efficiency ('), 
                            bolditalic(k[cat] / K[m]), 
                            bold(') [Strader, 2001]'), sep = ''))
  ) +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_text(face = 'bold', size = 20), 
        axis.title.y = element_text(face = 'bold', size = 22),
        axis.text = element_text(size = 18), 
        legend.position = 'top', 
        legend.justification = 0.5, 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)
  ) + ylim(-0.4, 0.4)
p_strader2001_newfig

# p_figs6 <- plot_grid(p_dam2000_newfig + theme(legend.position = 'none'), 
#                      p_strader2001_newfig, ncol = 2, labels = c('A', 'B'), 
#                      label_size = 20, label_fontface = 'bold')
# p_figs6

# ggsave(plot = p_figs6, device = cairo_pdf, width = 17, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS6_mutants_literature_exp_level_newAxis.pdf')
# ggsave(plot = p_figs6, device = 'png', width = 17, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS6_mutants_literature_exp_level_newAxis.png')

p_figs6 <- plot_grid(p_dam2000_newfig + theme(legend.position = 'none'), 
                     p_strader2001_newfig, nrow = 2, labels = c('A', 'B'), 
                     label_size = 20, label_fontface = 'bold')
ggsave(plot = p_figs6, device = cairo_pdf, width = 8, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS6_mutants_literature_exp_level_newAxis.pdf')
ggsave(plot = p_figs6, device = 'png', width = 8, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS6_mutants_literature_exp_level_newAxis.png')


#### Figure S7: Validations and the DMS scores are correlated ####

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
  
  # Subtract the blank for this experiment and multiply by 5 to make it OD / mL
  # (experiment was carried out in 0.2 mL)
  data.pl %<>% mutate(OD = (OD - 0.088) * 5)
  
  d <-  select(data.pl, 1:3)
  d %<>% spread(key = "Well", value = "OD", convert = F) 
  colnames(d)[1] <- "time"
  
  ## Use Growthcurver as follows
  ## I stop at t = 14.4 because that is the time it takes the last WT culture to get to its max OD
  gc_out <- SummarizeGrowthByPlate(d, t_trim =14.4)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
  
  
  # print(data.pl)
}

data.od1 <- read.my.gc(file.od, plate.ind)

# Subtract the blank for this experiment and multiply by 5 to make it OD / mL
# (experiment was carried out in 0.2 mL)
# data.od1 %<>% mutate(OD = (OD - 0.087) * 5)

data.od1$ID <- factor(data.od1$ID, levels = unique(data.od1$ID))
lbls <- select(data.od1, c(4,11)) %>% distinct()

data.od1 %<>% filter(!is.na(Replicate)) %>%
  unite(col = ID_new, Mutant, Arabinose, TMP, sep='_', remove = FALSE)
data.od1$ID_new <- factor(data.od1$ID_new, levels = unique(data.od1$ID_new))

# # Check how much time it takes the WT curves to reach their maximum
# data.od1.max <- data.od1 %>% 
#   group_by(ID, TMP, Replicate) %>%
#   filter(OD == max(OD)) %>%
#   filter(Mutant == 'WT')
# ## The WT curves reach their maximum between t = 12.3 and 14.4. I will modify the code above.

# Adjust the IDs of the growthcurver output and summarise
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

all_data_complete_new <- all_data_complete %>% filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose) %>%
  rowwise() %>%
  mutate(ID = str_c(Position, WT_Residue, Residue, '_', Arabinose, sep = '')) %>%
  mutate(ID = ifelse(ID == '2EE_0.01', 'WT_0.01', ID)) %>%
  mutate(ID = ifelse(ID == '2EE_0.2', 'WT_0.2', ID))

# Join the data
joined_sets <- inner_join(x = all_data_complete_new, 
                          # y = data.od.all %>% filter(TMP == 10),
                          y = data.od1.new %>% filter(TMP == 10),
                          by = c('ID' = 'ID')) %>%
  mutate(Arabinose = as.factor(Arabinose))

joined_sets %<>% mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                                           ifelse(Arabinose == 0.025, 'Suboptimal', 
                                                  ifelse(Arabinose == 0.05, 'Near-optimal', 
                                                         ifelse(Arabinose == 0.2, 'Optimal', NA)))))

# Repeat with the area under the curve (empiric)
p_figS7a <- joined_sets %>% ungroup() %>% rowwise() %>%
  mutate(Position = toString(Position)) %>%
  mutate(mut_id = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  mutate(mut_id = ifelse(mut_id == 'E2E', 'WT', mut_id), 
         exp_level = factor(exp_level, levels = c('Weak', 'Optimal'))) %>%
  ggplot(aes(x = mean_sel_coeff, y = auc_e)) +
  geom_smooth(method = 'lm', alpha = 0.5) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.7, size = 6,
           label.y.npc = 0.15, method = 'spearman', cor.coef.name = 'rho'
  ) +
  geom_point(aes(colour = exp_level, label = mut_id), size = 2) +
  geom_errorbarh(aes(xmax = mean_sel_coeff + sem_sel_coeff, xmin = mean_sel_coeff - sem_sel_coeff, 
                     colour = exp_level)) +
  geom_label_repel(aes(colour = exp_level, label = mut_id), 
                   fontface = 'bold', 
                   fill = '#999999',
                   # fill = '#737373',
                   max.overlaps = Inf, show.legend = F, alpha = 1,
                   size = 5
                   # fill = '#000000'
  ) +
  scale_colour_manual(values = c('#fed976', '#80001A')) + 
  # scale_fill_manual(values = c('#000000', '#ffffff')) +
  # xlab('s') + 
  # ylab('Growth in liquid culture (AUC)') +
  labs(y = expression(paste(bold('Growth in liquid culture ('), 
                            bolditalic('AUC'), 
                            bold(')'), sep = '')), 
       x = expression(bolditalic(s))) +
  labs(colour = 'Expression level') +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20),
        legend.position = 'top', 
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'))
p_figS7a

#### Figure S7B: Boxplots showing differences in fitness effects ####

joined_sets %<>% rowwise() %>% mutate(mut_id = str_c(Position, Residue, sep = '')) %>%
  mutate(WT_check = (WT_Residue == Residue))

comps <- compare_means(auc_e~exp_level, data = joined_sets, paired = TRUE) %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16', 
                           sprintf("p = %2.1e", as.numeric(p))
  ),
  y_pos = c(57)
  )

p_figS7b <- joined_sets %>% ungroup() %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Optimal'))) %>%
  ggplot(aes(x = exp_level, y = auc_e, fill = exp_level)) +
  geom_boxplot() + 
  geom_point(aes(colour = WT_check, size = WT_check)) + 
  geom_line(aes(group = mut_id, colour = WT_check)) +
  scale_fill_manual(values = c('#fed976', '#80001A')) +
  scale_colour_manual(values = c('black', 'blue')) +
  scale_size_manual(values = c(1, 2)) +
  labs(fill = '') +
  # stat_compare_means(method = 'wilcox.test', paired = TRUE, 
  #                    comparisons = list(c('Weak', 'Optimal')), 
  #                    aes(
  #                      # label = ifelse(p < 0.05,sprintf("p = %2.1e", as.numeric(..p.format..)), ..p.format..))
  #                      label = paste('p = ', ..p.format.., sep = ''))
  # ) +
  geom_signif(data = as.data.frame(comps),
              inherit.aes = FALSE, aes(xmin = group1, xmax = group2,
                                       annotations=p.format, y_position = y_pos), 
              manual = TRUE, textsize = 7) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20),
        legend.position = 'none', 
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) + 
  # ylab('Growth in liquid culture (AUC)') +
  labs(y = expression(paste(bold('Growth in liquid culture ('), 
                            bolditalic('AUC'), 
                            bold(')'), sep = ''))) +
  xlab('Expression level')
p_figS7b

# Plot them together
p_figS7 <- plot_grid(p_figS7a, p_figS7b, labels = c('A', 'B'), 
                     label_size = 20, label_fontface = 'bold', ncol = 2)

ggsave(plot = p_figS7, device = cairo_pdf, width = 16, height = 8, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS7_validations.pdf')
ggsave(plot = p_figS7, device = png, width = 16, height = 8, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS7_validations.png')


#### Figure S8: Volcano plots from ANOVAs, mean deltaS per position ####

#### Figure S8A: Mean deltaS per position ####

data_fig_2 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose)

data_fig2c <- data_fig_2 %>% ungroup() %>%
  group_by(Position, Arabinose) %>%
  summarise(meanSelCoeff = mean(mean_sel_coeff))

data_fig2c_exp <- data_fig2c %>% 
  mutate(Expression_level = ifelse(Arabinose == 0.01, 'Weak', 
                                   ifelse(Arabinose == 0.025, 'Suboptimal', 
                                          ifelse(Arabinose == 0.05, 'Near-optimal', 
                                                 ifelse(Arabinose == 0.2, 'Optimal', 
                                                        ifelse(Arabinose == 0.4, 'Overexpressed', NA)))))) %>%
  mutate(Expression_level = factor(Expression_level,
                                   levels = c('Weak', 'Suboptimal',
                                              'Near-optimal', 'Optimal', 'Overexpressed')))

## Draw the figure
p_figs8a <- data_fig2c_exp %>% 
  ggplot(aes(x = Expression_level, y = meanSelCoeff, fill = Expression_level)) +
  geom_violin(alpha = 0.8) +
  stat_summary(fun="median", geom="point", size = 5) +
  geom_point() +
  geom_line(aes(group = Position)) +
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
  # ylab('Mean s per position')
  labs(y = expression(paste(bold('Mean '), bolditalic('s'), bold(' per position'), sep = '')))
p_figs8a

#### Figure S8B: Volcano plots ####

data_fig_2 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose)

# Separate the data for optimal expression from the rest
data_part_1 <- data_fig_2 %>%
  # filter(Arabinose == 0.01)
  filter(Arabinose == 0.2)

data_part_2 <- data_fig_2 %>%
  # filter(Arabinose != 0.01)
  filter(Arabinose != 0.2)

# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "sem_sel_coeff_2", "Arabinose_2")

# Join
data_fig_2_final <- inner_join(x = data_part_1, y = data_part_2, 
                               by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
)


data_fig2b <- data_fig_2_final %>% ungroup() %>%
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff,
         Arabinose_2 = str_c(Arabinose_2, '% arabinose', sep = ''))


# Relevel the data set
data_fig2b %<>% mutate(Arabinose_2 = factor(Arabinose_2, 
                                            # levels = c('0.025% arabinose', '0.05% arabinose', '0.2% arabinose'))
                                            levels = c('0.01% arabinose', '0.025% arabinose', '0.05% arabinose', '0.4% arabinose'))
)

data_fig2b_exp <- data_fig2b %>% separate(col = Arabinose_2, into = c('Arabinose_num', 'tmp'), 
                                          sep = '% arabinose') %>%
  mutate(Arabinose_num = as.numeric(Arabinose_num), 
         exp_level = ifelse(Arabinose_num == 0.01, 'Weak', 
                            ifelse(Arabinose_num == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose_num == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose_num == 0.2, 'Optimal',
                                                 ifelse(Arabinose_num == 0.4, 'Overexpressed', NA))))))

# Load data for all replicates from both sequencers
all_data_all_reps <- read_delim('../R67_DMS_February2022/Complete_datasets_MiSeq_NovaSeq/all_data_all_reps_bothSequencers.txt', 
                                delim = '\t')

all_data_all_reps_TMP10 <- all_data_all_reps %>% rowwise() %>%
  filter(TMP == 10, Timepoint == 10, WT_Residue != Residue, Residue != '*') %>%
  arrange(Position, Residue, WT_Residue) %>%
  mutate(Genotype = str_c(WT_Residue, Position, Residue)) %>%
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
  
  # Save values from the ANOVA
  df_arabinose <- anova_test$Df[1]
  df_residuals <- anova_test$Df[2]
  
  sum_sq_arabinose <- anova_test$`Sum Sq`[1]
  sum_sq_residuals <- anova_test$`Sum Sq`[2]
  
  mean_sq_arabinose <- anova_test$`Mean Sq`[1]
  mean_sq_residuals <- anova_test$`Mean Sq`[2]
  
  f_value <- anova_test$`F value`[1]
  p_value <- anova_test$`Pr(>F)`[1]
  
  new_row <- data.frame(Genotype = c(genotype), df_arabinose = c(df_arabinose), df_residuals = c(df_residuals),
                        sum_sq_arabinose = c(sum_sq_arabinose), sum_sq_residuals = c(sum_sq_residuals),
                        f_value = c(f_value), p_val = c(p_value))
  all_anovas <- bind_rows(all_anovas, new_row)
}

## Apply the Benjamini-Hochberg correction for multiple hypotheses
fdr_test <-  p.fdr(pvalues = all_anovas$p_val, adjust.method = 'BH', threshold = 0.05)
all_anovas$p.adj <- fdr_test$fdrs

summary(all_anovas$p.adj)
sum(all_anovas$p.adj < 0.05)

# Prepare the data
temp_data <- data_fig_2 %>% rowwise() %>%
  mutate(ID = str_c(Position, Residue, sep = ''))

# Let's add a column with the difference in scores at low and high expression
score_diff <- data_fig2b %>% 
  filter(Arabinose_2 == '0.01% arabinose') %>%
  # mutate(mut_check_diff = ifelse(abs(diffNormScore) > 0.3, TRUE, FALSE))
  mutate(mut_check_diff = ifelse(abs(diffNormScore) > 0.1, TRUE, FALSE))

score_significant <- left_join(x = data_fig2b_exp %>% filter(Arabinose_num == 0.01), 
                               y = all_anovas %>% rowwise() %>%
                                 # separate(Genotype_2, into = c('WT_Residue', 'Residue', 'Position'), 
                                  #         sep = c(1, 2)) %>%
                                 tidyr::extract(col = Genotype, into = c('WT_Residue', 'Position', 'Residue'), 
                                         regex = '([A-Z])([0-9]+)([A-Z])') %>%
                                 mutate(Position = as.numeric(Position)), 
                               by = c('WT_Residue' = 'WT_Residue', 'Position' = 'Position', 
                                      'Residue' = 'Residue')) %>%
  mutate(mut_check = p.adj < 0.05)

table(score_diff$mut_check_diff)
table(score_significant$mut_check)
table(score_diff$mut_check_diff, score_significant$mut_check)

# Join the data to add these marks
data_fig2d <- left_join(x = temp_data, 
                        y= score_significant %>% select(Position, Residue, p.adj, mut_check),
                        by = c('Position' = 'Position', 'Residue' = 'Residue'))


# Add the expression level
data_fig2d_exp <- data_fig2d %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal',
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA))))))

# Add the data about the minimum change in expression
data_fig2d_exp <- left_join(x = data_fig2d_exp, 
                            y = score_diff %>% select(Position, Residue, diffNormScore, mut_check_diff), 
                            by = c('Position' = 'Position', 'Residue' = 'Residue')) %>%
  rowwise() %>%
  mutate(mut_check_final = and(mut_check, mut_check_diff))

## Try a volcano plot 
p_figs8b <- data_fig2d_exp %>%
  filter(Arabinose == 0.01) %>%
  ggplot(aes(x = diffNormScore, y = -log10(p.adj), colour = mut_check_final)) +
  geom_point() +
  # geom_vline(xintercept = -0.3, linetype = 'dashed') +
  # geom_vline(xintercept = 0.3, linetype = 'dashed') +
  geom_vline(xintercept = -0.1, linetype = 'dashed') +
  geom_vline(xintercept = 0.1, linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  # xlab('\u0394s') + 
  # ylab('-log10(p.adj)') +
  labs(
    # x = expression(s[weak] - s[opt]), 
    x = expression(paste(bold('\u0394'), bolditalic(s[weak]),
                     bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
                     sep = '')),
       y = expression(paste(bold('-log10('), 
                            bolditalic('p.adj'), 
                            bold(')'), sep = ''))) +
  theme(axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.title.x = element_text(size = 22, face = 'bold'),
        axis.text = element_text(size = 18), 
        legend.position = 'none')
p_figs8b  

write.table(x = all_anovas, append = F, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Figures/2022-05-09_Supp_figures_paper/TableSX.ANOVA_individual_mutants.csv')

p_figs8 <- plot_grid(p_figs8a, p_figs8b, ncol = 2, labels = c('A', 'B'), label_size = 20, 
          label_fontface = 'bold')
p_figs8
ggsave(p_figs8, device = cairo_pdf, width = 20, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS8.volcano_plot_signif_mutants.pdf')
ggsave(p_figs8, device = 'png', width = 20, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS8.volcano_plot_signif_mutants.png')

### Do the ANOVA on ranks for all the mutants and all replicates
data_anova <- all_data_all_reps_TMP10 %>% 
  mutate(ID = as.factor(ID),
         Arabinose = as.factor(Arabinose), 
         Genotype = as.factor(Genotype))

## Takes a long time to run
m <- art(sel_coeff ~ Arabinose*Genotype + Error(ID), data=data_anova)
summary(m)
anova_general <- anova(m)

anova_general_table <- data.frame(anova_general) %>%
  mutate(Pr..F. = ifelse(Pr..F. == 0, '<2.2e-16', Pr..F.))
write.table(x = anova_general_table, append = F, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Figures/2022-05-09_Supp_figures_paper/TableSX.ANOVA_general.csv')

#### Fig. S9: k-means clustering ####

# Join the data to add these marks
data_fig1f <- left_join(x = temp_data, 
                        y= score_significant %>% select(Position, Residue, p.adj, mut_check),
                        by = c('Position' = 'Position', 'Residue' = 'Residue'))


# Add the expression level
data_fig1f_exp <- data_fig1f %>% 
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA))))))

# Add the data about the minimum change in expression
data_fig1f_exp <- left_join(x = data_fig1f_exp, 
                            y = score_diff %>% select(Position, Residue, diffNormScore, mut_check_diff), 
                            by = c('Position' = 'Position', 'Residue' = 'Residue')) %>%
  rowwise() %>%
  mutate(mut_check_final = and(mut_check, mut_check_diff))

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

p_figS9A <- kmeans_ss %>% 
  ggplot(aes(x = k, y = tot_within_ss)) +
  geom_point(size = 3) + 
  geom_line() +
  scale_x_continuous(breaks = 1:10) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  # xlab('Number of clusters (k)') + 
  labs(x = expression(paste(bold('Number of clusters ('), 
                        bolditalic('k'), bold(')'), sep = ''))) +
  ylab('Sum of squared errors') +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = 'top',
        legend.justification = 0.5,
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank())
p_figS9A

# Run the k-means clustering with the seed
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

#### Draw Fig. S9B ####

cluster_relabel_new <- cluster_relabel %>% select(-old_cluster) %>%
  pivot_longer(cols = c('Weak', 'Suboptimal', 'Near-optimal', 'Optimal', 'Overexpressed'), 
               names_to = 'exp_level', values_to = 'mean_sel_coeff')

p_figS9B <- cluster_relabel_new %>% 
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                                  'Optimal', 'Overexpressed'))) %>%
  ggplot(aes(x = exp_level, y = mean_sel_coeff, colour = as.factor(new_cluster))) +
  geom_point(size = 3) +
  geom_line(aes(group = as.factor(new_cluster)), size = 2) +
  scale_colour_manual(values = c('#1b9e77', '#ffbf80', '#7570b3', '#e7298a')) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = 'top',
        legend.justification = 0.5,
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  xlab('Expression level') +
  # ylab('s') +
  labs(colour = 'Cluster', y = expression(bolditalic(s)))
p_figS9B

## Prepare the data for Fig. S9C

# Load the table with protein sites
data_interfaces_final <- read_delim('Data/Complete_datasets_TMP0_TMP10/data_annotation_2.txt', delim = '\t')

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

data_fig1f_new_summ <- data_fig1f_new %>% ungroup() %>% 
  group_by(Position, WT_Residue, Residue) %>%
  summarise(Cluster = mean(cluster))

# Add the data for protein regions
data_fig1f_new_summ <- left_join(x = data_fig1f_new_summ, 
                                 y = data_interfaces_final %>% rowwise() %>%
                                   mutate(Unannotated = ifelse(sum(`A,C`, `A,D`, DHF, NADPH,
                                                                   Cat_residues, Disordered_region,
                                                                   Buried) == 0, 
                                                               1, 0)
                                   ), 
                                 by = c('Position' = 'Position'))

#### Run chi square tests with contingency tables for one protein site at a time ####

## Only for the dimerization interface
data_fig1f_new_dim_int <- data_fig1f_new_summ  %>% 
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
                        'Cat_residues', 'Disordered_region', 'Buried', 
                        'Unannotated'), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'A,C')

cont_table <- table(data_fig1f_new_dim_int$Site_check, data_fig1f_new_dim_int$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
dim_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the tetramerization interface
data_fig1f_new_tet <- data_fig1f_new_summ  %>% 
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
                        'Cat_residues', 'Disordered_region', 'Buried', 
                        'Unannotated'), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'A,D')

cont_table <- table(data_fig1f_new_tet$Site_check, data_fig1f_new_tet$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
tet_log2fold <- log2((chisq$observed + 1)/ chisq$expected)[2,]

## Only for DHF binding
data_fig1f_dhf <- data_fig1f_new_summ  %>% 
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
                        'Cat_residues', 'Disordered_region', 'Buried', 
                        'Unannotated'), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'DHF')

cont_table <- table(data_fig1f_dhf$Site_check, data_fig1f_dhf$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
dhf_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the catalytic residues
data_fig1f_new_cat <- data_fig1f_new_summ  %>% 
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
                        'Cat_residues', 'Disordered_region', 'Buried', 
                        'Unannotated'), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'Cat_residues')

cont_table <- table(data_fig1f_new_cat$Site_check, data_fig1f_new_cat$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
cat_log2fold <- log2((chisq$observed + 1)/ chisq$expected)[2,]

## Only for the NADPH binding residues
data_fig1f_new_nadph <- data_fig1f_new_summ  %>% 
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
                        'Cat_residues', 'Disordered_region', 'Buried', 
                        'Unannotated'), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'NADPH')

cont_table <- table(data_fig1f_new_nadph$Site_check, data_fig1f_new_nadph$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
nadph_log2fold <- log2((chisq$observed + 1)/ chisq$expected)[2,]

## Only for the disordered region
data_fig1f_new_dis <- data_fig1f_new_summ  %>% 
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
                        'Cat_residues', 'Disordered_region', 'Buried', 
                        'Unannotated'), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'Disordered_region')

cont_table <- table(data_fig1f_new_dis$Site_check, data_fig1f_new_dis$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
dis_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the buried residues
data_fig1f_new_bur <- data_fig1f_new_summ  %>% 
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
                        'Cat_residues', 'Disordered_region', 'Buried', 
                        'Unannotated'), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'Buried')

cont_table <- table(data_fig1f_new_bur$Site_check, data_fig1f_new_bur$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
bur_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the unannotated residues
data_fig1f_new_unannot <- data_fig1f_new_summ  %>% 
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH',
                        'Cat_residues', 'Disordered_region', 'Buried', 
                        'Unannotated'), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'Unannotated')

cont_table <- table(data_fig1f_new_unannot$Site_check, data_fig1f_new_unannot$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
unannot_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

#### Add the ddG data ####
data_fig1f_new_summ_ddg <- left_join(x = data_fig1f_new_summ, 
                                     y = all_data_complete %>% select(Position, WT_Residue, Residue,
                                                                      Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, 
                                                                      Mean_ddG_int_HM_A_D) %>%
                                       unique(),
                                     by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                                            'Residue' = 'Residue'))

data_fig1f_new_summ_ddg %<>% filter(!(is.na(Mean_ddG_stab_HET))) %>%
  mutate(high_ddG_stab = ifelse(Mean_ddG_stab_HET >= 2, 1, 0), 
         low_ddG_stab = ifelse(Mean_ddG_stab_HET < 2, 1, 0), 
         
         high_ddG_dim = ifelse(Mean_ddG_int_HM_A_C >= 2, 1, 0), 
         low_ddG_dim = ifelse(Mean_ddG_int_HM_A_C < 2, 1, 0), 
         
         high_ddG_tet = ifelse(Mean_ddG_int_HM_A_D >= 2, 1, 0), 
         low_ddG_tet = ifelse(Mean_ddG_int_HM_A_D < 2, 1, 0))

## Only for the mutants with high ddG stability
data_fig1f_high_ddg_stab <- data_fig1f_new_summ_ddg  %>% 
  pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
                        high_ddG_dim, low_ddG_dim, 
                        high_ddG_tet, low_ddG_tet), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'high_ddG_stab')

cont_table <- table(data_fig1f_high_ddg_stab$Site_check, data_fig1f_high_ddg_stab$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
high_ddg_stab_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the mutants with low ddG stability
data_fig1f_low_ddg_stab <- data_fig1f_new_summ_ddg  %>% 
  pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
                        high_ddG_dim, low_ddG_dim, 
                        high_ddG_tet, low_ddG_tet), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'low_ddG_stab')

cont_table <- table(data_fig1f_low_ddg_stab$Site_check, data_fig1f_low_ddg_stab$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
low_ddg_stab_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the mutants with high ddG dimerization interface
data_fig1f_high_ddg_dim <- data_fig1f_new_summ_ddg  %>% 
  pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
                        high_ddG_dim, low_ddG_dim, 
                        high_ddG_tet, low_ddG_tet), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'high_ddG_dim')

cont_table <- table(data_fig1f_high_ddg_dim$Site_check, data_fig1f_high_ddg_dim$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
high_ddg_dim_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the mutants with low ddG stability
data_fig1f_low_ddg_dim <- data_fig1f_new_summ_ddg  %>% 
  pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
                        high_ddG_dim, low_ddG_dim, 
                        high_ddG_tet, low_ddG_tet), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'low_ddG_dim')

cont_table <- table(data_fig1f_low_ddg_dim$Site_check, data_fig1f_low_ddg_dim$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
low_ddg_dim_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the mutants with high ddG tetramerization interface
data_fig1f_high_ddg_tet <- data_fig1f_new_summ_ddg  %>% 
  pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
                        high_ddG_dim, low_ddG_dim, 
                        high_ddG_tet, low_ddG_tet), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'high_ddG_tet')

cont_table <- table(data_fig1f_high_ddg_tet$Site_check, data_fig1f_high_ddg_tet$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
high_ddg_tet_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

## Only for the mutants with low ddG stability
data_fig1f_low_ddg_tet <- data_fig1f_new_summ_ddg  %>% 
  pivot_longer(cols = c(high_ddG_stab, low_ddG_stab, 
                        high_ddG_dim, low_ddG_dim, 
                        high_ddG_tet, low_ddG_tet), names_to = 'Site',
               values_to = 'Site_check') %>%
  filter(Site == 'low_ddG_tet')

cont_table <- table(data_fig1f_low_ddg_tet$Site_check, data_fig1f_low_ddg_tet$Cluster)
chisq <- chisq.test(cont_table)
chisq$observed
chisq$expected
chisq$p.value

# Calculate log2(obs / exp)
low_ddg_tet_log2fold <- log2((chisq$observed + 1) / chisq$expected)[2,]

#### Concatenate the data for all sites and put it in a heatmap ####

# c('A,C', 'A,D', 'DHF', 'NADPH',
#   'Cat_residues', 'Disordered_region', 'Buried', 
#   'Unannotated')
data_enrichment <- rbind(dim_log2fold, tet_log2fold, dhf_log2fold, cat_log2fold, 
                         nadph_log2fold, dis_log2fold, bur_log2fold, unannot_log2fold, 
                         high_ddg_dim_log2fold, high_ddg_tet_log2fold, high_ddg_stab_log2fold, 
                         low_ddg_dim_log2fold, low_ddg_tet_log2fold, low_ddg_stab_log2fold)

rownames(data_enrichment) <- c('Dimerization interface', 'Tetramerization interface', 
                               'DHF binding', 'Catalytic residues', 'NADPH binding', 
                               'Disordered region',
                               'Buried residues', 'Unannotated residues', 
                               '\u0394\u0394G dim >= 2', '\u0394\u0394G stab >= 2',
                               '\u0394\u0394G tet >= 2', 
                               '\u0394\u0394G dim < 2', '\u0394\u0394G stab < 2',
                               '\u0394\u0394G tet < 2')

seq1 <- seq(-8, 0, length.out = 4)
seq2 <- seq(0, 2, length.out = 4)

p_figS9C <- Heatmap(data_enrichment, cluster_columns = F, cluster_rows = T,
                        clustering_distance_rows = 'pearson',
                        ## Color ramp for the old normalized scores
                        col = colorRamp2(
                          breaks = c(seq1[1:(length(seq1)-1)], 0, seq2[2:length(seq2)]),
                          colors = rev(brewer.pal(n = 7, name = 'RdBu'))),
                        show_column_names = T, row_names_side = 'left',
                        # width=unit(31, 'cm'), height = unit(31, 'cm'),
                        width = unit(17, 'cm'), height = unit(10, 'cm'),
                        border = T,
                        row_title = '',
                        row_title_gp = gpar(fontsize=20, fontface = 'bold'),
                        row_names_rot = 0,
                        row_names_centered = T,
                        row_names_gp = gpar(fontsize=14, fontface = 'bold'),
                        column_title = 'Cluster',
                        column_title_side = 'bottom',
                        column_names_rot = 0,
                        column_names_gp = gpar(fontsize=14,fontface='bold'),
                        column_title_gp = gpar(fontsize=20, fontface = 'bold'),
                        # right_annotation = ha_empty,
                        show_heatmap_legend = TRUE,
                        heatmap_legend_param = list(
                          at = c(-8, -4, 0, 1, 2),
                          title = "log2 fold enrichment",
                          title_gp = gpar(fontsize = 20),
                          legend_height = unit(5, "cm"),
                          legend_width = unit(2, "cm"),
                          border='black',
                          lwd=1.7,
                          labels_gp = gpar(fontsize = 18),
                          title_position = "leftcenter-rot"
                        )
)
p_figS9C

p_figS9 <- plot_grid(
  p_figS9A + theme(plot.margin = margin(t = 0, b = 1, r = 6, l = 3, unit = 'cm')),
  p_figS9B + theme(plot.margin = margin(t = 0, b = 1, r = 6, l = 3, unit = 'cm')),
  grid.grabExpr(draw(p_figS9C)), nrow = 3,
  labels = c('A', 'B', 'C'), label_size = 20, label_fontface = 'bold')

ggsave(p_figS9, width = 14, height = 17, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS9_kmeans_clusters.pdf')
ggsave(p_figS9, width = 14, height = 17, dpi = 300, device = 'png', 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS9_kmeans_clusters.png')

#### Figure S10: Faceted view of figure 1E with errorbars ####

data_fig_2 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose)

# Separate the data for optimal expression from the rest
data_part_1 <- data_fig_2 %>%
  # filter(Arabinose == 0.01)
  filter(Arabinose == 0.2)

data_part_2 <- data_fig_2 %>%
  # filter(Arabinose != 0.01)
  filter(Arabinose != 0.2)

# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "sem_sel_coeff_2", "Arabinose_2")

# Join
data_fig_2_final <- inner_join(x = data_part_1, y = data_part_2, 
                               by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
)

data_fig_2_final_exp <- data_fig_2_final %>% 
  mutate(exp_level = ifelse(Arabinose_2 == 0.01, 'Weak expression', 
                            ifelse(Arabinose_2 == 0.025, 'Suboptimal expression', 
                                   ifelse(Arabinose_2 == 0.05, 'Near-optimal expression',
                                          ifelse(Arabinose_2 == 0.4, 'Overexpressed', NA))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak expression', 'Suboptimal expression', 
                                                  'Near-optimal expression', 'Overexpressed')))


## Make sure WT appears only once
data_fig_2_final_exp_nowt <- data_fig_2_final_exp %>% filter(WT_Residue != Residue)
data_fig_2_final_exp_wt <- data_fig_2_final_exp %>% filter(WT_Residue == Residue, Position == 2)

data_fig_2_final_exp_plot <- bind_rows(data_fig_2_final_exp_nowt, data_fig_2_final_exp_wt)

# Define a function for npc annotations
annotation_box <- roundrectGrob(x = unit(0.225, 'npc'), y = unit(0.85, 'npc'), 
                                width = unit(0.3, 'npc'), height = unit(0.125, 'npc'), 
                                gp = gpar(lwd = 2, col = 'black',
                                          # fill = '#999999', 
                                          fill = '#737373',
                                          lty = 2))

p_figs10 <- data_fig_2_final_exp_plot %>%
  ggplot(aes(x = mean_sel_coeff, y = mean_sel_coeff_2, 
             colour = as.factor(exp_level))
  ) +
  facet_wrap(~exp_level, nrow = 2, scales = 'free') +
  geom_errorbarh(aes(xmax = mean_sel_coeff + sem_sel_coeff, xmin = mean_sel_coeff - sem_sel_coeff), 
                 alpha = 0.5) +
  geom_errorbar(aes(ymax = mean_sel_coeff_2 + sem_sel_coeff_2, ymin = mean_sel_coeff_2 - sem_sel_coeff_2), 
                alpha = 0.5) +
  geom_point(alpha = 1, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = 'black') +
  # xlab('s (optimal expression)') + ylab('s (non-optimal expression)') + 
  labs(colour = 'Expression level in y-axis', 
       x = expression(paste(bolditalic('s'), 
                            bold(' (optimal expression)'), 
                            sep = '')), 
       y = expression(paste(bolditalic('s'), 
                            bold(' (non-optimal expression)')))) +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', 'black')) +
  annotation_custom(annotation_box) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
           # label.x.npc = 'left',
           # label.y.npc = 'top',
           label.x.npc = 0.05,
           label.y.npc = 0.9,
           size = 8,
           vjust = 1,
           method = 'spearman', show.legend = F, cor.coef.name = 'rho'# , 
           # geom = 'label', fill = '#999999'
  ) +
  geom_smooth(method = 'loess', show.legend = F) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = 'top', 
        legend.justification = 0.5, 
        strip.text = element_text(size = 18, face = 'bold'),
        strip.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()
  ) +
  xlim(-1, 0.5) + ylim(-1, 0.5)
p_figs10
ggsave(p_figs10, device = cairo_pdf, width = 21, height = 17, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS10_facet_errorbars.pdf')
ggsave(p_figs10, device = 'png', width = 21, height = 17, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS10_facet_errorbars.png')




##### Fig. S11: deltaAUC correlates with deltaS ####

## Have a look at the growth curves ##
plate.ind <- 'Data/2022-04-25_exp_increase/Growth_curves/Growth_curves_20_04_2022_sampleSheet.csv'
file.od <- 'Data/2022-04-25_exp_increase/Growth_curves/Growth_curves_20_04_2022_data.xlsx'

## Define function to read plate data
# function to process plates ----------------------------------------------
read.my.gc <- function(file, plate.index){
  pl <- read.xlsx(file,sheetIndex = 1, rowIndex = 3:62, stringsAsfactors = FALSE,
                  header = F)
  ind <- read_delim(plate.ind, delim = '\t', locale = locale(decimal_mark = ',')) 
  
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
  data.pl %<>% mutate(OD = (OD - 0.087) * 5)
  
  d <-  select(data.pl, 1:3)
  d %<>% spread(key = "Well", value = "OD", convert = F) 
  colnames(d)[1] <- "time"
  
  ## Use Growthcurver
  gc_out <- SummarizeGrowthByPlate(d, t_trim =20)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
}

data.od1 <- read.my.gc(file.od, plate.ind)

# Let's see how the area under the curve correlates with the DMS selection coefficients
# Use a left join to add the selection coefficients
data.od1.new <- left_join(x = data.od1, 
                          y = all_data_complete %>% rowwise() %>% 
                            mutate(ID = str_c(WT_Residue, Position, Residue)) %>%
                            mutate(ID = ifelse(ID == 'E2E', 'WT', ID)) %>%
                            filter(Timepoint == 10, TMP == 10) %>%
                            select(ID, Arabinose, mean_sel_coeff),
                          by = c('ID' = 'ID', 'Arabinose concentration' = 'Arabinose')
)
# Get the mean and standard deviation of replicates
data.od1.summary <- data.od1.new %>% ungroup() %>%
  filter(ID != 'Blank') %>%
  group_by(ID, `Arabinose concentration`, Well) %>%
  summarise(mean_auc = mean(auc_e), mean_sel_coeff = mean(mean_sel_coeff)) %>%
  ungroup() %>% group_by(`Arabinose concentration`, ID) %>%
  summarise(mean_auc_final = mean(mean_auc), sem_auc = sd(mean_auc)/sqrt(n()),
            num_samples = n(), mean_sel_coeff = mean(mean_sel_coeff))
# 
# p <- data.od1.summary %>%
#   ggplot(aes(x = mean_sel_coeff, y = mean_auc_final, 
#              ymax = mean_auc_final + sem_auc, ymin = mean_auc_final - sem_auc, 
#              colour = as.factor(`Arabinose concentration`))) +
#   geom_point(size = 3) +
#   geom_errorbar() +
#   scale_colour_manual(values = c('#FED976', '#80001A')) +
#   geom_label_repel(aes(label = ID), fontface = 'bold', fill = 'black', alpha = 0.6, 
#                    show.legend = F) + 
#   theme(legend.position = 'top',
#         axis.text.x =  element_text(size = 18), 
#         axis.text.y = element_text(size = 18),
#         axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         axis.line=element_line(), 
#         legend.justification = 'center') +
#   xlab('DMS sel. coeff.') + ylab('Growth in liquid culture (AUC)') +
#   labs(colour = 'Arabinose') + guides(label = 'none', fill = 'none')
# p
# ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_growth_crossing.pdf')

# ### Need to check the growth curves
# p <- data.od1.new %>% filter(ID != 'Blank') %>%
#   ggplot(aes(x = time, y = OD, colour = ID)) +
#   facet_grid(ID~`Arabinose concentration`) +
#   geom_line(aes(group = Well))
# p
# ggsave(p, device = cairo_pdf, width = 10, height = 14, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_growth_crossing_curves.pdf')

#### Show the data as a reaction norm ####

## Add labels for the expression level
data.od1.summary %<>% 
  mutate(Expression_level = ifelse(`Arabinose concentration` == 0.01, 'Weak', 
                                   ifelse(`Arabinose concentration` == 0.2, 'Optimal', NA))) %>%
  mutate(Expression_level = factor(Expression_level, levels = c('Weak', 'Optimal')))

# p1 <- data.od1.summary %>% 
#   # mutate(Arabinose = as.factor(`Arabinose concentration`)) %>%
#   ggplot(aes(x = Expression_level,
#              y = mean_auc_final, 
#              ymax = mean_auc_final + sem_auc, ymin = mean_auc_final - sem_auc)) +
#   geom_point() + 
#   # geom_errorbar(width = 0.025) +
#   geom_label_repel(aes(label = ID, colour = Expression_level), fontface = 'bold', 
#                    fill = '#999999') +
#   geom_line(aes(group = ID)) +
#   scale_colour_manual(values = c('#FED976', '#80001A')) +
#   theme(legend.position = 'none',
#         axis.text.x =  element_text(size = 18), 
#         axis.text.y = element_text(size = 18),
#         axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         axis.line=element_line(), 
#         legend.justification = 'center') +
#   xlab('Expression level') + ylab('Growth in liquid culture (AUC)')
# p1
# 
# p2 <- data.od1.summary %>% 
#   # mutate(Arabinose = as.factor(`Arabinose concentration`)) %>%
#   ggplot(aes(x =Expression_level, y = mean_sel_coeff)) +
#   geom_point() + 
#   geom_label_repel(aes(label = ID, colour = Expression_level), fontface = 'bold', 
#                    fill = '#999999') +
#   scale_colour_manual(values = c('#FED976', '#80001A')) +
#   geom_line(aes(group = ID)) +
#   theme(legend.position = 'none',
#         axis.text.x =  element_text(size = 18), 
#         axis.text.y = element_text(size = 18),
#         axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         axis.line=element_line(), 
#         legend.justification = 'center') +
#   xlab('Expression level') + ylab('s')
# p2
# 
# p_final <- plot_grid(p2, p1, ncol = 2, labels = c('A', 'B'), 
#                      label_size = 20)
# p_final

## Alternatively, plot the delta(AUC) vs the deltaS

# Pivot the AUC values
data.od1.summary.auc <- data.od1.summary %>% ungroup() %>%
  select(-mean_sel_coeff, -sem_auc, -num_samples, -`Arabinose concentration`) %>%
  pivot_wider(names_from = Expression_level, values_from = mean_auc_final,
              names_prefix = 'AUC_') %>%
  mutate(delta_AUC = AUC_Weak - AUC_Optimal)

data.od1.summary.deltaS <- data.od1.summary %>% ungroup() %>%
  select(-mean_auc_final, -sem_auc, -num_samples, -`Arabinose concentration`) %>%
  pivot_wider(names_from = Expression_level, values_from = mean_sel_coeff, 
              names_prefix = 'selCoeff_') %>%
  mutate(diffNormScore = selCoeff_Weak - selCoeff_Optimal)

data.od1.summary.final <- inner_join(x = data.od1.summary.auc, 
                                     y = data.od1.summary.deltaS, 
                                     by = c('ID' = 'ID'))

p_figS11 <- data.od1.summary.final %>% 
  ggplot(aes(x = diffNormScore, y = delta_AUC)) +
  geom_point(size = 2) +
  geom_label_repel(aes(label = ID), fontface = 'bold', size = 6) +
  theme(axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.text = element_text(size = 18), 
        legend.title = element_text(size = 20),
        legend.position = 'none', 
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white'),
        plot.margin = margin(t = 0.5, b = 0.5, l = 0.5, r = 0.5, unit = 'cm')
        ) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 6,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  geom_smooth(method = 'lm', show.legend = F) +
  labs(x = expression(paste(bold('\u0394'),
                            bolditalic(s[weak]), 
                            bold(' ('),
                            bolditalic(s[weak] - s[opt]), 
                            bold(')'), sep = '')),
       y = expression(paste(bold('\u0394'),
                            bolditalic(AUC[weak]), 
                            bold(' ('),
                            bolditalic(AUC[weak] - AUC[opt]),
                            bold(')'), sep = '')))
p_figS11
ggsave(p_figS11, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS11_deltaAUC_deltaS.pdf')
ggsave(p_figS11, device = 'png', width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS11_deltaAUC_deltaS.png')

####  Fig. S12: GEMME vs DMS (with TMP) ####

#### Compare the DMS data and the GEMME data ####

# Load GEMME data
gemme_data <- read_delim('GEMME/DfrB1_normPred_evolCombi_rownames.txt', delim = ' ')

# Use gather to get a long data frame
gemme_data_long <- gemme_data %>% gather(-Residue, key = Position_str, value = Fitness) %>%
  mutate(Position = as.numeric(substr(Position_str, 2, length(Position_str))),
         Residue = toupper(Residue))

# Use an inner join to put the data together with the fitness scores
gemme_vs_dms <- inner_join(x = all_data_complete %>% filter(Timepoint == 10) %>%
                             select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, TMP),
                           y = gemme_data_long %>% select(-Position_str) %>% filter(!(is.na(Fitness))),
                           by = c('Position' = 'Position', 'Residue' = 'Residue'))

## Make sure WT appears only once
gemme_vs_dms_nowt <- gemme_vs_dms %>% filter(WT_Residue != Residue)
gemme_vs_dms_wt <- gemme_vs_dms %>% filter(WT_Residue == Residue, Position == 2)

gemme_vs_dms_plot <- bind_rows(gemme_vs_dms_nowt, gemme_vs_dms_wt)

## Prepare plots
p_figS12 <- gemme_vs_dms_plot %>% rowwise() %>% 
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak expression', 
                            ifelse(Arabinose == 0.025, 'Suboptimal expression', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal expression', 
                                          ifelse(Arabinose == 0.2, 'Optimal expression',
                                                 ifelse(Arabinose == 0.4, 'Overexpression', NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak expression', 'Suboptimal expression', 
                                                  'Near-optimal expression', 'Optimal expression',
                                                  'Overexpression'))) %>%
  filter(TMP == 10) %>%
  mutate(Arabinose = str_c(toString(Arabinose), '% arabinose', sep = '')) %>%
  ggplot(aes(x = Fitness, y = mean_sel_coeff)) + 
  geom_point() + 
  xlab('GEMME score') + # ylab('s') +
  labs(y = expression(bolditalic('s'))) +
  # facet_wrap(~Arabinose) +
  facet_wrap(~exp_level, scales = 'free', nrow = 3) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(face = 'bold', size = 20), axis.text = element_text(size = 18),
        strip.text = element_text(size = 20, face = 'bold'), 
        strip.background =  element_rect(fill = 'white'), 
        axis.line = element_line()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 6,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  geom_smooth(method = 'lm', show.legend = F) +
  ylim(-1, 0.7)
# ylim(-1, 2)
p_figS12
ggsave(plot = p_figS12, device = cairo_pdf, width = 14, height = 17, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS12_GEMME.pdf')
ggsave(plot = p_figS12, device = 'png', width = 14, height = 17, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS12_GEMME.png')


#### Figure S13: Expresion-dependent differences in fitness effects become weaker ####
#### as expression approaches the optimum ####

data_fig_4 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 10, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose)

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
annotation_name_gp = gpar(fontface = 'bold', fontsize = 18),
annotation_name_side = 'left',
annotation_name_rot = 90,
annotation_name_offset = c(Entropy = '0.15cm')
)

### 0.01 arabinose
# Separate the data for ara 0.2
data_part_1 <- data_fig_4 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_4 %>%
  filter(Arabinose == 0.01)


## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)


# Need to convert to wide formatted data
data_fig_4_final_df <- data_fig_4_final %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_4_final <- as.matrix(data_fig_4_final_df %>% select(-Position))

rownames(data_fig_4_final) <- data_fig_4_final_df$Position

# Get a matrix of true/false values for the synonymous codons (I will reuse the)
fig_4_bool <- data_fig_4 %>%
  filter(Arabinose == 0.025) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_4_bool_final <- as.matrix(fig_4_bool %>% select(-Position))

rownames(fig_4_bool_final) <- fig_4_bool$Position

# Need to reorder the columns in the matrices
data_fig_4_final <- data_fig_4_final[1:nrow(data_fig_4_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_4_bool_final <- fig_4_bool_final[1:nrow(fig_4_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

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

### Repeat figure 4a with the annotation about expression level
p_figs13_ara0.01 <- Heatmap(
  t(data_fig_4_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-4, 4, length.out = 7) / 10, 
                   colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  show_heatmap_legend = F, 
  row_title = "Residue",
  row_title_gp = gpar(fontsize=22, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  top_annotation = ha1,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_4_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-4, -2, 0, 2, 4) / 10,
    title = "\u0394s", 
    title_gp = gpar(fontsize = 20),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 18),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)
p_figs13_ara0.01

### 0.025 arabinose
# Separate the data for ara 0.2
data_part_1 <- data_fig_4 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.025
data_part_2 <- data_fig_4 %>%
  filter(Arabinose == 0.025)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Need to convert to wide formatted data
data_fig_4_final_df <- data_fig_4_final %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_4_final <- as.matrix(data_fig_4_final_df %>% select(-Position))

rownames(data_fig_4_final) <- data_fig_4_final_df$Position

# Get a matrix of true/false values for the synonymous codons (I will reuse the)
fig_4_bool <- data_fig_4 %>%
  filter(Arabinose == 0.025) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_4_bool_final <- as.matrix(fig_4_bool %>% select(-Position))

rownames(fig_4_bool_final) <- fig_4_bool$Position

# Need to reorder the columns in the matrices
data_fig_4_final <- data_fig_4_final[1:nrow(data_fig_4_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_4_bool_final <- fig_4_bool_final[1:nrow(fig_4_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p_figs13_ara0.025 <- Heatmap(
  t(data_fig_4_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-4, 4, length.out = 7) / 10, 
                   colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  show_heatmap_legend = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=22, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_4_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-4, -2, 0, 2, 4) / 10,
    # title = "\u0394s", 
    title = expression(paste(bold('\u0394'), 
                             bolditalic('s'), sep = '')),
    title_gp = gpar(fontsize = 20, fontface = 'bold'),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 18),
    title_position = "leftcenter-rot"
  ), 
  column_labels= column_labels
)

p_figs13_ara0.025

### 0.05 arabinose
# Separate the data for ara 0.2
data_part_1 <- data_fig_4 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.05
data_part_2 <- data_fig_4 %>%
  filter(Arabinose == 0.05)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Need to convert to wide formatted data
data_fig_4_final_df <- data_fig_4_final %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_4_final <- as.matrix(data_fig_4_final_df %>% select(-Position))

rownames(data_fig_4_final) <- data_fig_4_final_df$Position

# Get a matrix of true/false values for the synonymous codons (I will reuse the)
fig_4_bool <- data_fig_4 %>%
  filter(Arabinose == 0.05) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_4_bool_final <- as.matrix(fig_4_bool %>% select(-Position))

rownames(fig_4_bool_final) <- fig_4_bool$Position

# Need to reorder the columns in the matrices
data_fig_4_final <- data_fig_4_final[1:nrow(data_fig_4_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_4_bool_final <- fig_4_bool_final[1:nrow(fig_4_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p_figs13_ara0.05 <- Heatmap(
  t(data_fig_4_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-4, 4, length.out = 7) / 10, 
                   colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  show_heatmap_legend = F,
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=22, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  # bottom_annotation = ha2,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_4_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-4, -2, 0, 2, 4) / 10,
    title = "\u0394s", 
    title_gp = gpar(fontsize = 20, fontface = 'bold'),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 18),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)

p_figs13_ara0.05

### 0.4 arabinose
# Separate the data for ara 0.2
data_part_1 <- data_fig_4 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.4
data_part_2 <- data_fig_4 %>%
  filter(Arabinose == 0.4)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Need to convert to wide formatted data
data_fig_4_final_df <- data_fig_4_final %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_4_final <- as.matrix(data_fig_4_final_df %>% select(-Position))

rownames(data_fig_4_final) <- data_fig_4_final_df$Position

# Get a matrix of true/false values for the synonymous codons (I will reuse the)
fig_4_bool <- data_fig_4 %>%
  filter(Arabinose == 0.4) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_4_bool_final <- as.matrix(fig_4_bool %>% select(-Position))

rownames(fig_4_bool_final) <- fig_4_bool$Position

# Need to reorder the columns in the matrices
data_fig_4_final <- data_fig_4_final[1:nrow(data_fig_4_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_4_bool_final <- fig_4_bool_final[1:nrow(fig_4_bool_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

data_interfaces_final <- read_delim('Data/Complete_datasets_TMP0_TMP10/data_annotation_2.txt', delim = '\t')

# col_fun = colorRamp2(c(0, 1), c("white", "#595959"))

# # Prepare the annotation on interfaces
# # Interface A-C will be interface 1
# # Interface A-D will be interface 2
# ha2 <- HeatmapAnnotation(`Interface 1` = data_interfaces_final$`A,C`,
#                          `Interface 2` = data_interfaces_final$`A,D`,
#                          DHF = data_interfaces_final$DHF,
#                          NADPH = data_interfaces_final$NADPH,
#                          `Catalytic residues` = data_interfaces_final$Cat_residues,
#                          `Disordered region` = data_interfaces_final$Disordered_region,
#                          show_annotation_name = T,
#                          annotation_name_gp = gpar(fontface = 'bold', fontsize = 14),
#                          annotation_name_side = 'left',
#                          show_legend = FALSE,
#                          col = list(`Interface 1` = colorRamp2(c(0, 1), c("white", "#009E73")), 
#                                     `Interface 2` = colorRamp2(c(0, 1), c("white", "#0072B2")),
#                                     DHF = colorRamp2(c(0, 1), c("white", "#E69F00")),
#                                     NADPH = colorRamp2(c(0, 1), c("white", "#D55E00")),
#                                     `Catalytic residues` = colorRamp2(c(0, 1), c("white", "#56B4E9")),
#                                     `Disordered region` = colorRamp2(c(0, 1), c("white", "#CC79A7"))
#                          ),
#                          gp = gpar(col = "black")
# )

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


# Color order
# c('Disordered region', 'Catalytic residues', 'DHF', 'Interface 1', 'Interface 2', 'NADPH', 'Unannotated')
# c('#CC79A7', '#56B4E9', '#E69F00', '#009E73', '#0072B2', '#D55E00', '#000000')

p_figs13_ara0.4 <- Heatmap(
  t(data_fig_4_final), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-4, 4, length.out = 7) / 10, 
                   colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  show_heatmap_legend = F,
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=22, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  bottom_annotation = ha2,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_4_bool_final[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-4, -2, 0, 2, 4) / 10,
    title = "\u0394s", 
    title_gp = gpar(fontsize = 20, fontface = 'bold'),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 18),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)
p_figs13_ara0.4


## Put the three figures together
ht_list = p_figs13_ara0.01 %v% p_figs13_ara0.025 %v% p_figs13_ara0.05 %v% p_figs13_ara0.4 
p_figs13_heatmaps <- grid.grabExpr(
  draw(ht_list,
       row_title_gp = gpar(fontsize=20, fontface = 'bold'),
       ht_gap = unit(1, "cm"))
)


### Add text labels
## Arabinose 0.01
text_fig_ara0.01 <- ggplot() + 
  draw_label(# expression(paste('s'[weak], ' - s'[opt], sep = '')),
    expression(atop(paste(bold('\u0394'), bolditalic(s[weak]), sep = ''),
                    paste(bold('('), 
                     bolditalic(s[weak] - s[opt]), bold(')'), sep = ''))),
             # x = 0.7, y = 0.45, 
    x = 0.7, y = 0.3,
             fontface = 'bold', size = 35, angle = 90, colour = '#fed976') +
  theme(axis.line = element_blank())
text_fig_ara0.01

# Arabinose 0.025
text_fig_ara0.025 <- ggplot() + draw_label(# expression(paste('s'[subopt], ' - s'[opt], sep = '')),
  expression(atop(paste(bold('\u0394'), bolditalic(s[subopt]), sep = ''),
                  paste(bold('('), 
                        bolditalic(s[subopt] - s[opt]), bold(')'), sep = ''))),
                                           # x = 0.7, y = 0.55, 
  x = 0.7, y = 0.4,
                                           fontface = 'bold', size = 35, angle = 90, colour = '#fd8d3c') +
  theme(axis.line = element_blank())
text_fig_ara0.025

# Arabinose 0.05
text_fig_ara0.05 <- ggplot() + draw_label(# expression(paste('s'[near-opt], ' - s'[opt], sep = '')),
  expression(atop(paste(bold('\u0394'), bolditalic(s[near-opt]), sep = ''),
                  paste(bold('('), 
                        bolditalic(s[near-opt] - s[opt]), bold(')'), sep = ''))),
                                          # x = 0.7, y = 0.7, 
  x = 0.7, y = 0.4,
                                          fontface = 'bold', size = 35, angle = 90, colour = '#bd0026') +
  theme(axis.line = element_blank())
text_fig_ara0.05

# Arabinose 0.4
text_fig_ara0.4 <- ggplot() + draw_label(# expression(paste('s'[over], ' - s'[opt], sep = '')),
  expression(atop(paste(bold('\u0394'), bolditalic(s[over]), sep = ''),
                  paste(bold('('), 
                        bolditalic(s[over] - s[opt]), bold(')'), sep = ''))),
                                          # x = 0.7, y = 0.8, 
  x = 0.7, y = 0.35,
                                          fontface = 'bold', size = 35, angle = 90, colour = 'black') +
  theme(axis.line = element_blank())
text_fig_ara0.4

## Add a panel with the distributions of deltaS
# Separate the data for ara 0.2
data_part_1 <- data_fig_4 %>%
  filter(Arabinose == 0.2)

# Separate the rest of the data
data_part_2 <- data_fig_4 %>%
  filter(Arabinose != 0.2)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_deltaS_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Add expression level
data_deltaS_final %<>%
  mutate(Expression_level = ifelse(Arabinose_2 == 0.01, 'Weak', 
                                   ifelse(Arabinose_2 == 0.025, 'Suboptimal', 
                                          ifelse(Arabinose_2 == 0.05, 'Near-optimal',
                                                 ifelse(Arabinose_2 == 0.4, 'Overexpressed', NA))))) %>%
  mutate(Expression_level = factor(Expression_level,
                                   levels = c('Weak', 'Suboptimal', 'Near-optimal', 'Overexpressed')))

comps <- compare_means(diffNormScore~Expression_level, data = data_deltaS_final,
                       paired = TRUE) %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16', 
                           sprintf("p = %2.1e", as.numeric(p))
  ),
  y_pos = c(0.6, 0.75, 0.9, 1.05, 1.20, 1.35)
  )

p_figs13b <- data_deltaS_final %>% ungroup() %>% 
  ggplot(aes(x = Expression_level, y = diffNormScore, fill = Expression_level)) +
  geom_point(aes(colour = Expression_level),alpha = 1, 
             position = position_jitterdodge(jitter.width = 0.25)) +
  geom_violin(aes(fill = Expression_level),alpha = 0.5) +
  geom_signif(data = as.data.frame(comps), inherit.aes = FALSE, aes(xmin = group1, xmax = group2,
                                                                    annotations=p.format, y_position = y_pos), 
              manual = TRUE, textsize = 6) +
  stat_summary(fun="median", geom="point", size = 5, colour = 'grey') +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', 'black')) +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', 'black')) +
  # ylab('\u0394s') +  # Delta fitness
  labs(y = expression(paste(bold('\u0394'), bolditalic('s'), sep = ''))) +
  xlab('Expression level') +
  theme(axis.title = element_text(face = 'bold', size = 24), 
        axis.text = element_text(size = 22), 
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = 'none', 
        legend.justification = 0.5, 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
  plot.margin = margin(t = 1, r = 9, b = 0, l = 15, 'cm'))
p_figs13b

### Put everything together
# p_figs13_text <- plot_grid(text_fig_ara0.01, text_fig_ara0.025, text_fig_ara0.05, text_fig_ara0.4,
#                            nrow = 4, rel_heights = c(1.1, 1, 1.1, 1))
# 
# p_figs13a <- plot_grid(p_figs13_text + theme(plot.margin = margin(t = 0, r = -3, b = 0, l = 8, 'cm')),
#                        p_figs13_heatmaps, ncol = 2, rel_widths = c(0.2, 1), 
#                        labels = c('', 'A'), label_size = 20, label_fontface = 'bold')
# 
# p_figs13b_null <- plot_grid(NULL, p_figs13b, rel_widths = c(0.2, 1), ncol = 2, 
#                             labels = c('', 'B'), label_size = 20, label_fontface = 'bold')
# 
# p_figs13 <- plot_grid(p_figs13a, p_figs13b_null, nrow = 2, rel_heights = c(1, 0.3))
# 
# ggsave(p_figs13, width = 23, height = 31, dpi = 300, device = cairo_pdf, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS13_supp_allDiff_noSecStruc_buried.pdf')

## Save text and figures separately
p_figs13_text <- plot_grid(text_fig_ara0.01, text_fig_ara0.025, text_fig_ara0.05, text_fig_ara0.4,
                           nrow = 4, rel_heights = c(1.1, 1, 1.1, 1))

p_figs13 <- plot_grid(p_figs13_heatmaps, p_figs13b, nrow = 2, rel_heights = c(1, 0.3), 
                      labels = c('A', 'B'), label_size = 20, label_fontface = 'bold')

ggsave(p_figs13, width = 23, height = 31, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS13_supp_allDiff_noSecStruc_buried.pdf')
ggsave(p_figs13_text, width = 23, height = 31, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS13_supp_allDiff_noSecStruc_buried_text.pdf')

# Plot correlation of (s_nearopt - s_opt) vs (s_over - s_opt)
data_corr_check <- data_deltaS_final %>%
  ungroup() %>% group_by(Position, WT_Residue, Residue,
                         Arabinose, Arabinose_2) %>%
  select(-Expression_level, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Arabinose_2, names_prefix = 'ara', values_from = diffNormScore)

## Comparing the deltaS from overexpression to the others
cor.test(data_corr_check$ara0.01, data_corr_check$ara0.4, method = 'spearman')
cor.test(data_corr_check$ara0.025, data_corr_check$ara0.4, method = 'spearman')
cor.test(data_corr_check$ara0.05, data_corr_check$ara0.4, method = 'spearman')

cor(data_corr_check %>% ungroup() %>% select(ara0.01, ara0.025, ara0.05, ara0.4),
    method = 'spearman')

## Near-optimal vs overexpressed (deltaS with Sopt comparison)
# p <- data_corr_check %>% 
#   ggplot(aes(y = ara0.05, x = ara0.4)) +
#   geom_point() +
#   stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
#            label.x.npc = 0.6, label.y.npc = 0.2, size = 6, 
#            cor.coef.name = 'r', method = 'pearson') +
#   geom_smooth(method = 'lm') +
#   theme(axis.text = element_text(size = 18), 
#         axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank()
#         ) +
#   labs(y = expression(S[near-opt] - S[opt]), x = expression(S[over] - S[opt]))
# p
# ggsave(p, width = 10, height = 7, dpi = 300, device = cairo_pdf, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigSX_deltaS_over_nearopt_anticorrelation.pdf')
# 
# ## deltaS comparison (ara0.01 vs ara0.4)
# p <- data_corr_check %>% 
#   ggplot(aes(y = ara0.01, x = ara0.4)) +
#   geom_point() +
#   stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
#            label.x.npc = 0.6, label.y.npc = 0.2, size = 6, 
#            cor.coef.name = 'r', method = 'pearson') +
#   geom_smooth(method = 'lm') +
#   theme(axis.text = element_text(size = 18), 
#         axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank()
#   ) +
#   labs(y = expression(S[weak] - S[opt]), x = expression(S[over] - S[opt]))
# p
# ggsave(p, width = 10, height = 7, dpi = 300, device = cairo_pdf, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigSX_deltaS_over_weak_anticorrelation.pdf')
# 
# ## deltaS comparison (ara0.025 vs ara0.4)
# p <- data_corr_check %>% 
#   ggplot(aes(y = ara0.025, x = ara0.4)) +
#   geom_point() +
#   stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
#            label.x.npc = 0.6, label.y.npc = 0.2, size = 6, 
#            cor.coef.name = 'r', method = 'pearson') +
#   geom_smooth(method = 'lm') +
#   theme(axis.text = element_text(size = 18), 
#         axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank()
#   ) +
#   labs(y = expression(S[subopt] - S[opt]), x = expression(S[over] - S[opt]))
# p
# ggsave(p, width = 10, height = 7, dpi = 300, device = cairo_pdf, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigSX_deltaS_over_subopt_anticorrelation.pdf')


#### Figure S14: E2R and E2V increase expression ####

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

# Draw the figures: columns will be arabinose concentrations, rows will be groups
# One figure for each timepoint, replicates will appear on the same panel
# p <- all_data_processed %>%
#   filter(Timepoint == 'Stationary phase') %>%
#   ggplot(aes(x = as.numeric(GRNBHLog))) +
#   facet_rep_grid(Name ~ Arabinose, scales = 'free', repeat.tick.labels = T) +
#   geom_density(aes(colour = as.factor(Replicate), group = as.factor(Replicate)),
#                size = 1, alpha = 0.3) +
#   xlab('GFP fluorescence') + ylab('Density')  +
#   # geom_histogram(aes(colour = as.factor(Replicate), group = as.factor(Replicate))) +
#   xlim(0, 3)  +
#   ylim(0, 6) +
#   theme(legend.position = 'none',
#         strip.background = element_rect(fill = 'white'), strip.text = element_text(face = 'bold', size = 16),
#         axis.text =  element_text(size = 16), axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         axis.line=element_line()) +
#   ggtitle('Stationary phase')
# p
# ggsave(p, device = cairo_pdf, width = 17, height = 17, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_2ER_cytometry_statPhase.pdf')
# ggsave(p, device = 'png', width = 17, height = 17, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_2ER_cytometry_statPhase.png')
# 
# p <- all_data_processed %>%
#   filter(Timepoint == 'Exponential phase') %>%
#   ggplot(aes(x = as.numeric(GRNBHLog))) +
#   facet_rep_grid(Name ~ Arabinose, scales = 'free', repeat.tick.labels = T) +
#   geom_density(aes(colour = as.factor(Replicate), group = as.factor(Replicate)),
#                size = 1, alpha = 0.3) +
#   xlab('GFP fluorescence') + ylab('Density')  +
#   # geom_histogram(aes(colour = as.factor(Replicate), group = as.factor(Replicate))) +
#   xlim(0, 3)  +
#   ylim(0, 8.5) +
#   theme(legend.position = 'none',
#         strip.background = element_rect(fill = 'white'), strip.text = element_text(face = 'bold', size = 16),
#         axis.text =  element_text(size = 16), axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         axis.line=element_line()) +
#   ggtitle('Exponential phase')
# p
# ggsave(p, device = cairo_pdf, width = 17, height = 17, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_2ER_cytometry_expPhase.pdf')
# ggsave(p, device = 'png', width = 17, height = 17, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_2ER_cytometry_expPhase.png')

#### New version, one curve for each group, multiple groups per panel ####

# p <- all_data_processed %>%
#   mutate(Name = factor(Name, levels = c('DfrB1', 'DfrB1-GFP', 'DfrB1-M26-GFP',
#                                         'DfrB1-2ER-GFP', 'DfrB1-2ER-M26-GFP','GFP'))) %>%
#   ggplot(aes(x = as.numeric(GRNBHLog))) +
#   facet_rep_grid(Arabinose ~ Timepoint, scales = 'free', repeat.tick.labels = T) +
#   geom_density(aes(colour = as.factor(Name), group = as.factor(Name)),
#                size = 1, alpha = 0.3) +
#   xlab('GFP fluorescence') + ylab('Density')  +
#   # geom_histogram(aes(colour = as.factor(Replicate), group = as.factor(Replicate))) +
#   xlim(0, 3)  +
#   ylim(0, 8.5) +
#   scale_colour_manual(values = c('#000000', # DfrB1
#                                  '#D55E00', # DfrB1-GFP
#                                  '#CC79A7', # DfrB1-M26-GFP
#                                  '#56B4E9', # DfrB1-2ER-GFP
#                                  '#0072B2', # DfrB1-2ER-M26-GFP
#                                  '#009E73' # GFP
#   )) +
#   theme(legend.position = 'top',
#         strip.background = element_rect(fill = 'white'), strip.text = element_text(face = 'bold', size = 16),
#         axis.text =  element_text(size = 16), axis.title = element_text(face = 'bold', size = 20),
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         plot.title = element_text(hjust = 0.5), 
#         axis.line=element_line(), 
#         legend.justification = 'center') +
#   labs(colour = 'Construct')
# p
# ggsave(p, device = cairo_pdf, width = 17, height = 17, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_2ER_cytometry_both.pdf')
# ggsave(p, device = 'png', width = 17, height = 17, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_2ER_cytometry_both.png')

## Draw boxplots ##

# p_figS14 <- all_data_processed %>%
#   filter(Timepoint == 'Exponential phase') %>%
#   mutate(Name = ifelse(Name == 'DfrB1-2ER-GFP', 'DfrB1-E2R-GFP', 
#                        ifelse(Name == 'DfrB1-2ER-M26-GFP', 'DfrB1-E2R-M26-GFP', Name))) %>%
#   mutate(Name = factor(Name, levels = c('DfrB1', 'DfrB1-GFP', 'DfrB1-M26-GFP',
#                                         'DfrB1-E2R-GFP', 'DfrB1-E2R-M26-GFP','GFP')), 
#          Arabinose = str_c(Arabinose, '% arabinose', sep = '')) %>%
#   mutate(Arabinose = factor(Arabinose, 
#                             levels = c('0% arabinose', '0.01% arabinose', '0.2% arabinose'))) %>%
#   ggplot(aes(x = Name, y = as.numeric(GRNBHLog))) +
#   # geom_jitter(alpha = 0.3, width = 0.2) +
#   # geom_violin(alpha = 0.6) + 
#   geom_boxplot(alpha = 1, outlier.shape = NA, colour = 'black', fill = '#bfbfbf') + 
#   # facet_wrap(~Expression_level, scales = 'free') +
#   facet_wrap(~Arabinose, scales = 'free') +
#   # scale_colour_manual(values = c('#000000', # DfrB1
#   #                                '#D55E00', # DfrB1-GFP
#   #                                '#CC79A7', # DfrB1-M26-GFP
#   #                                '#56B4E9', # DfrB1-2ER-GFP
#   #                                '#0072B2', # DfrB1-2ER-M26-GFP
#   #                                '#009E73' # GFP
#   # )) +
#   # scale_fill_manual(values = c('#000000', # DfrB1
#   #                                '#D55E00', # DfrB1-GFP
#   #                                '#CC79A7', # DfrB1-M26-GFP
#   #                                '#56B4E9', # DfrB1-2ER-GFP
# #                                '#0072B2', # DfrB1-2ER-M26-GFP
# #                                '#009E73' # GFP
# # )) +
# theme(legend.position = 'none',
#       strip.background = element_rect(fill = 'white'), strip.text = element_text(face = 'bold', size = 16),
#       axis.text.x =  element_text(size = 16, angle = 45, hjust = 1), 
#       axis.text.y = element_text(size = 16),
#       axis.title = element_text(face = 'bold', size = 20),
#       panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#       panel.grid.minor = element_blank(), 
#       plot.title = element_text(hjust = 0.5), 
#       axis.line=element_line(), 
#       legend.justification = 'center') +
#   xlab('Construct') + ylab('GFP fluorescence (log10)') +
#   labs(fill = 'Construct', colour = 'Construct')
# p_figS14

# Prepare the panel for the constructs
p_figS14_constructs <- ggdraw() + 
  # draw_image('Figures/2022-04-19_extra_figures/Figure_S12_expressionIncrease_2ER_constructs.png')
  # draw_image('Figures/2022-04-19_extra_figures/Figure_S12_expressionIncrease_2ER_nocolors.png')
  draw_image('Figures/2022-05-09_Supp_figures_paper/Figure_S14_expressionIncrease_2ER.png')
p_figS14_constructs

## Load the data on enzymatic activity
infile = 'Data/Activity_assays_Kiana_2022-05-30/WT,E2R,E2V activity.xlsx'
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

p_figs14a <- activity_data %>% 
  mutate(Mutant = factor(Mutant, levels = c('WT', 'E2R', 'E2V'))) %>%
  ggplot(aes(x = Mutant, y = Avg_rel_activity)) +
  geom_bar(fill = '#737373', stat = 'identity') +
  geom_errorbar(aes(ymax = Avg_rel_activity + sem_rel_activity, 
                    ymin = Avg_rel_activity - sem_rel_activity), 
                width = 0.2) +
  theme(axis.title = element_text(size = 20, face = 'bold'), 
        axis.text = element_text(size = 18), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  geom_signif(data = as.data.frame(comps), inherit.aes = FALSE,
              aes(xmin = group1, xmax = group2, annotations=p.format, y_position = y_pos), 
              manual = TRUE, textsize = 6) +
  # stat_compare_means(method = 't.test', comparisons = list(c('WT', 'E2R'), 
  #                                                          c('WT', 'E2V'),
  #                                                          c('E2R', 'E2V'))) +
  geom_hline(yintercept = 100, linetype = 'dashed') +
  xlab('Mutant') + ylab('Relative activity (%)')
p_figs14a

## An alternative version of the figure on protein abundance
p_figS14 <- all_data_processed %>%
  filter(Timepoint == 'Exponential phase') %>%
  ## Adjust names
  mutate(Name = ifelse(Name == 'DfrB1-2ER-GFP', 'pBAD-dfrB1(E2R)-sfGFP', 
                      ifelse(Name == 'DfrB1-2ER-M26-GFP', 'pBAD-dfrB1[1-25](E2R)-sfGFP',
                             ifelse(Name == 'DfrB1', 'pBAD-dfrB1',
                                    ifelse(Name == 'DfrB1-GFP', 'pBAD-dfrB1-sfGFP',
                                           ifelse(Name == 'DfrB1-M26-GFP', 'pBAD-dfrB1[1-25]-sfGFP',
                                                  ifelse(Name == 'GFP', 'pBAD-sfGFP', Name))))))) %>%
  # mutate(Name = factor(Name, levels = c('DfrB1', 'DfrB1-GFP', 'DfrB1-M26-GFP',
  #                                       'DfrB1-E2R-GFP', 'DfrB1-E2R-M26-GFP','GFP')),
  mutate(Name = factor(Name, levels = c('pBAD-dfrB1', 'pBAD-dfrB1-sfGFP', 
                                        'pBAD-dfrB1[1-25]-sfGFP',
                                        'pBAD-dfrB1(E2R)-sfGFP',
                                        'pBAD-dfrB1[1-25](E2R)-sfGFP',
                                        'pBAD-sfGFP'
                                        )),
         Arabinose = factor(Arabinose, levels = c(0, 0.01, 0.2))) %>%
  ggplot(aes(x = Name, y = as.numeric(GRNBHLog), colour = Arabinose, fill = Arabinose)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) + 
  scale_colour_manual(values = c('grey', '#FED976', '#80001A')) +
  scale_fill_manual(values = c('grey', '#FED976', '#80001A')) +
  theme(legend.position = 'top',
      strip.background = element_rect(fill = 'white'), strip.text = element_text(face = 'bold', size = 16),
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
  xlab('Construct') + ylab('GFP fluorescence (log10)') +
  labs(fill = 'Arabinose (% m/v)', colour = 'Arabinose (% m/v)')
p_figS14

# p_figS13_final <- plot_grid(p_figS13_constructs, p_figS13, nrow = 2, rel_heights = c(0.5, 1))

p_figS14_abundance <- plot_grid(p_figS14_constructs, p_figS14, nrow = 2,
                                rel_heights = c(0.5, 1))
p_figS14_final <- plot_grid(p_figs14a, p_figS14_abundance, labels = c('A', 'B'), nrow = 2, 
                            label_size = 20, label_fontface = 'bold', 
                            rel_heights = c(0.65, 1))

# ggsave(p_figS13_final, device = cairo_pdf, width = 17, height = 10, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS13_2ER_cytometry_boxplot_new_nocolors.pdf')
# ggsave(p_figS13_final, device = 'png', width = 17, height = 10, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS13_2ER_cytometry_boxplot_new_nocolors.png')

ggsave(p_figS14_final, device = cairo_pdf, width = 10, height = 15, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS14_2ER_activity_abundance.pdf')
# ggsave(p_figS14_final, device = 'png', width = 10, height = 14, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS14_2ER_activity_abundance.png')

#### Supp. figure 15: F18 and P19 mutations ####

## Load panel A
p_figs15a <- ggdraw() + draw_image('Figures/2022-05-09_Supp_figures_paper/FigS15A_pymol_W45.png')
p_figs15a

## Panel B ##
# Load data
data_plddt <- read_delim('Figures/chimerax/AF2_Results/model1_pLDDT.txt', delim = '\t')

data_plddt$Residue <- factor(data_plddt$Residue, levels = c(data_plddt$Residue))

data_plddt %<>% mutate(color_check = (Residue %in% c('F18', 'P19')))

color_final <- ifelse(data_plddt$color_check == TRUE, 'red', 'black')

p_figs15b <- data_plddt %>% 
  ggplot(aes(x = Residue, y = pLDDT)) +
  theme(axis.title = element_text(face = 'bold', size = 20),
        axis.text.x = element_text(size = 18, angle = 90, vjust = 0.5,
                                   colour = color_final),
        axis.text.y = element_text(size = 18),
        panel.grid.major.y = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'none') +
  scale_y_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) +
  annotate(geom = 'rect', xmax = Inf, xmin = -Inf, ymax = 100, ymin = 90,
                     alpha = 0.4, fill = 'green') +
  annotate(geom = 'rect', xmax = Inf, xmin = -Inf, ymax = 90, ymin = 70,
           alpha = 0.4, fill = 'yellow') +
  annotate(geom = 'rect', xmax = Inf, xmin = -Inf, ymax = 70, ymin = 0,
           alpha = 0.4, fill = 'red') +
  geom_point(aes(colour = color_check), size = 3) +
  scale_colour_manual(values = c('black', 'red')) +
  geom_line(aes(group = '')) +
  geom_vline(xintercept = 21, linetype = 'dashed') +
  xlab('Residue') + 
  # ylab('AlphaFold2 pLDDT') +
  labs(y = expression(paste(bold('AlphaFold2 '), bolditalic(pLDDT), sep = ' ')))
p_figs15b

#### Figure S15C: Boxplots of effects per position (residues 16-26) ####

p_figs15c <- all_data_complete %>% 
  filter(TMP == 10, Position %in% seq(from = 16, to = 26, by = 1)) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Overexpressed',
                                                        NA)))))) %>%
  mutate(Position = str_c(WT_Residue, Position, sep = ''), 
         exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                                  'Optimal', 'Overexpressed'))) %>%
  mutate(Position = factor(Position, 
                           levels = c('F16', 'V17', 'F18', 'P19', 'S20', 'D21', 
                                      'A22', 'T23', 'F24', 'G25', 'M26'))) %>%
  ggplot(aes(x = Position, y = mean_sel_coeff,
             colour = exp_level, fill = exp_level)) +
  geom_boxplot(alpha = 0.4) +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'))+
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme(axis.title = element_text(face = 'bold', size = 20),
        axis.text = element_text(size = 18),
        panel.grid.major.y = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', 
        legend.justification = 'center', 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18)) +
  labs(colour = 'Expression level', fill = 'Expression level') +
  xlab('Position') +
  labs(y = expression(bolditalic(s)))
p_figs15c

# Put the panels together
p_figs15 <- plot_grid(p_figs15a, p_figs15b, p_figs15c, nrow = 3,
                      label_size = 20, labels = c('A', 'B', 'C'), label_fontface = 'bold')
ggsave(p_figs15, device = cairo_pdf, width = 21, height = 24, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS15_F18_P19_AF2.pdf')
ggsave(p_figs15, device = 'png', width = 21, height = 24, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS15_F18_P19_AF2.png')

#### Fig. S16 ####

#### Prepare figures for RSA and stability ####

# Separate the data for ara 0.2
data_part_1 <- data_fig_4 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_4 %>%
  filter(Arabinose == 0.01)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

summary(data_fig_4_final$diffNormScore)

# Need to convert to wide formatted data
data_fig_4_final_df <- data_fig_4_final %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_4_final <- as.matrix(data_fig_4_final_df %>% select(-Position))

rownames(data_fig_4_final) <- data_fig_4_final_df$Position
data_fig_4_final <- data_fig_4_final[1:nrow(data_fig_4_final), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

## Panel on solvent accessibility
## Now need a figure relating differences of effects (ara0.2 - ara0.01) to rSASA ##
data_ara0.2_ara0.01 <- all_data_complete %>% ungroup() %>%
  filter(Solvent_accessibility != 'Missing', Arabinose %in% c(0.01, 0.2), TMP == 10, Timepoint == 10) %>%
  select(Residue, Position, Arabinose, mean_sel_coeff, rSASA) %>%
  group_by(Residue, Position, rSASA) %>%
  pivot_wider(names_from = Arabinose, values_from = mean_sel_coeff) %>%
  mutate(diff_ara0.2_ara0.01 = `0.01` - `0.2`) %>%
  ungroup() %>% group_by(Position) %>%
  summarise(med_effect = median(diff_ara0.2_ara0.01, na.rm = T),
            low_pct = quantile(diff_ara0.2_ara0.01, probs = 0.25, na.rm = T),
            high_pct = quantile(diff_ara0.2_ara0.01, probs = 0.75, na.rm = T),
            rSASA = mean(rSASA))

# Draw the figure
p_diff_ara0.2_ara0.01_rsa <- data_ara0.2_ara0.01 %>%
  ggplot(aes(x = rSASA, y = med_effect, ymax = high_pct, ymin = low_pct)) +
  geom_point(size = 2) + geom_errorbar(alpha = 0.5) +
  # geom_text_repel(aes(label = Position)) +
  xlab('Relative solvent accessibility (with TMP)') + ylab('\u0394s') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
           label.x.npc = 0.7, label.y.npc = 0.2, size = 6, 
           cor.coef.name = 'rho', method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
           label.x.npc = 0.7, label.y.npc = 0.1, size = 6, 
           cor.coef.name = 'r', method = 'pearson') +
  geom_smooth(method = 'lm') +
  theme(axis.title = element_text(face = 'bold', size = 20),
        axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank())
p_diff_ara0.2_ara0.01_rsa

## Repeat the RSA figure for the data without TMP ##

## Panel on solvent accessibility
## Now need a figure relating differences of effects (ara0.2 - ara0.01) to rSASA ##
data_ara0.2_ara0.01 <- all_data_complete %>% ungroup() %>%
  filter(Solvent_accessibility != 'Missing', Arabinose %in% c(0.01, 0.2),
         TMP == 0, Timepoint == 10) %>%
  select(Residue, Position, Arabinose, mean_sel_coeff, rSASA) %>%
  group_by(Residue, Position, rSASA) %>%
  pivot_wider(names_from = Arabinose, values_from = mean_sel_coeff) %>%
  mutate(diff_ara0.2_ara0.01 = `0.01` - `0.2`) %>%
  ungroup() %>% group_by(Position) %>%
  summarise(med_effect = median(diff_ara0.2_ara0.01, na.rm = T),
            low_pct = quantile(diff_ara0.2_ara0.01, probs = 0.25, na.rm = T),
            high_pct = quantile(diff_ara0.2_ara0.01, probs = 0.75, na.rm = T),
            rSASA = mean(rSASA))

# Draw the figure
p_diff_ara0.2_ara0.01_rsa_notmp <- data_ara0.2_ara0.01 %>%
  ggplot(aes(x = rSASA, y = med_effect, ymax = high_pct, ymin = low_pct)) +
  geom_point(size = 2) + geom_errorbar(alpha = 0.5) +
  # geom_text_repel(aes(label = Position)) +
  xlab('Relative solvent accessibility (no TMP)') + ylab('\u0394s') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
           label.x.npc = 0.7, label.y.npc = 0.2, size = 6, 
           cor.coef.name = 'rho', method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, 
           label.x.npc = 0.7, label.y.npc = 0.1, size = 6, 
           cor.coef.name = 'r', method = 'pearson') +
  geom_smooth(method = 'lm') +
  theme(axis.title = element_text(face = 'bold', size = 20),
        axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank())
p_diff_ara0.2_ara0.01_rsa_notmp

p_diff_vs_rsa <- plot_grid(p_diff_ara0.2_ara0.01_rsa, p_diff_ara0.2_ara0.01_rsa_notmp, ncol = 2)
# ggsave(p_diff_vs_rsa, device = cairo_pdf, width = 17, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/Fig_RSA_deltaS.pdf')

## Panel on subunit stability
data_fig_4_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

data_fig_5 <- left_join(x = data_fig_4_final, 
                        y = all_data_complete %>% 
                          filter(TMP == 10, Arabinose == 0.01, Timepoint == 10) %>% # To avoid repetition
                          select(Position, WT_Residue, Residue, 
                                 Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, 
                                 Mean_ddG_int_HM_A_D, Mean_ddG_int_HM_A_B), 
                        by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                               'Residue' = 'Residue')) 

## Apply a log-modulus transformation to the data on stability
log_modulus <- function(x){
  return(sign(x) * log2(abs(x) + 1))
}

data_fig_5 %<>% 
  mutate(logm_diffScore = lapply(diffNormScore, FUN = log_modulus),
         logm_ddG = lapply(Mean_ddG_stab_HET, FUN = log_modulus))

data_fig_5 %<>%
  mutate(logm_diffScore = log_modulus(diffNormScore), 
         logm_ddG = log_modulus(Mean_ddG_stab_HET))

p_diff_vs_stab <- data_fig_5 %>%
  filter(!(is.na(Mean_ddG_stab_HET)), WT_Residue != Residue) %>%
  ggplot(aes(
    # x = logm_ddG, y = diffNormScore
    x = Mean_ddG_stab_HET, y = diffNormScore
    )
    ) +
  geom_point(size = 2) + 
  xlab('ddG on subunit stability [kcal / mol]') + ylab('\u0394s') +
  theme(axis.title = element_text(size = 20, face = 'bold'), axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.2, size = 6, 
           cor.coef.name = 'rho', method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.1, size = 6, 
           cor.coef.name = 'r', method = 'pearson') +
  geom_smooth(method = 'lm')
p_diff_vs_stab

## Repeat for the dimerization interface
p_diff_vs_dimInt <- data_fig_5 %>%
  filter(!(is.na(Mean_ddG_int_HM_A_C)), WT_Residue != Residue) %>%
  ggplot(aes(
    # x = logm_ddG, y = diffNormScore
    x = Mean_ddG_int_HM_A_C, y = diffNormScore
  )
  ) +
  geom_point(size = 2) + 
  xlab('ddG on dim. interface [kcal / mol]') + ylab('\u0394s') +
  theme(axis.title = element_text(size = 20, face = 'bold'), axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.2, size = 6, 
           cor.coef.name = 'rho', method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.1, size = 6, 
           cor.coef.name = 'r', method = 'pearson') +
  geom_smooth(method = 'lm')
p_diff_vs_dimInt

## Repeat for the dimerization interface
p_diff_vs_tetInt <- data_fig_5 %>%
  filter(!(is.na(Mean_ddG_int_HM_A_D)), WT_Residue != Residue) %>%
  ggplot(aes(
    # x = logm_ddG, y = diffNormScore
    x = Mean_ddG_int_HM_A_D, y = diffNormScore
  )
  ) +
  geom_point(size = 2) + 
  xlab('ddG on tet. interface (kcal / mol)') + ylab('\u0394s') +
  theme(axis.title = element_text(size = 20, face = 'bold'), axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.2, size = 6, 
           cor.coef.name = 'rho', method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.1, size = 6, 
           cor.coef.name = 'r', method = 'pearson') +
  geom_smooth(method = 'lm')
p_diff_vs_tetInt

p_ddGs <- plot_grid(p_diff_vs_stab, p_diff_vs_dimInt, p_diff_vs_tetInt, ncol = 3)
# ggsave(p_ddGs, device = cairo_pdf, width = 27, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/Fig_ddGs_TMP.pdf')

#### Repeat without TMP ####

data_fig_4 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 0, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose)

# Separate the data for ara 0.2
data_part_1 <- data_fig_4 %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2 <- data_fig_4 %>%
  filter(Arabinose == 0.01)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final <- left_join(x = data_part_1, y = data_part_2, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

data_fig_5 <- left_join(x = data_fig_4_final, 
                        y = all_data_complete %>% 
                          filter(TMP == 0, Arabinose == 0.01, Timepoint == 10) %>% # To avoid repetition
                          select(Position, WT_Residue, Residue, 
                                 Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, 
                                 Mean_ddG_int_HM_A_D, Mean_ddG_int_HM_A_B), 
                        by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
                               'Residue' = 'Residue')) 

## Apply a log-modulus transformation to the data on stability
log_modulus <- function(x){
  return(sign(x) * log2(abs(x) + 1))
}

data_fig_5 %<>% 
  mutate(logm_diffScore = lapply(diffNormScore, FUN = log_modulus),
         logm_ddG = lapply(Mean_ddG_stab_HET, FUN = log_modulus))

data_fig_5 %<>%
  mutate(logm_diffScore = log_modulus(diffNormScore), 
         logm_ddG = log_modulus(Mean_ddG_stab_HET))

p_diff_vs_stab <- data_fig_5 %>%
  filter(!(is.na(Mean_ddG_stab_HET))) %>%
  ggplot(aes(
    # x = logm_ddG, y = diffNormScore
    x = Mean_ddG_stab_HET, y = diffNormScore
  )
  ) +
  geom_point(size = 2) + 
  xlab('ddG on subunit stability [kcal / mol]') + ylab('\u0394s') +
  theme(axis.title = element_text(size = 20, face = 'bold'), axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.2, size = 6, 
           cor.coef.name = 'rho', method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.1, size = 6, 
           cor.coef.name = 'r', method = 'pearson') +
  geom_smooth(method = 'lm')
p_diff_vs_stab

## Repeat for the dimerization interface
p_diff_vs_dimInt <- data_fig_5 %>%
  filter(!(is.na(Mean_ddG_int_HM_A_C))) %>%
  ggplot(aes(
    # x = logm_ddG, y = diffNormScore
    x = Mean_ddG_int_HM_A_C, y = diffNormScore
  )
  ) +
  geom_point(size = 2) + 
  xlab('ddG on dim. interface [kcal / mol]') + ylab('\u0394s') +
  theme(axis.title = element_text(size = 20, face = 'bold'), axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.2, size = 6, 
           cor.coef.name = 'rho', method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.1, size = 6, 
           cor.coef.name = 'r', method = 'pearson') +
  geom_smooth(method = 'lm')
p_diff_vs_dimInt

## Repeat for the dimerization interface
p_diff_vs_tetInt <- data_fig_5 %>%
  filter(!(is.na(Mean_ddG_int_HM_A_D))) %>%
  ggplot(aes(
    # x = logm_ddG, y = diffNormScore
    x = Mean_ddG_int_HM_A_D, y = diffNormScore
  )
  ) +
  geom_point(size = 2) + 
  xlab('ddG on tet. interface (kcal / mol)') + ylab('\u0394s') +
  theme(axis.title = element_text(size = 20, face = 'bold'), axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.2, size = 6, 
           cor.coef.name = 'rho', method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01,
           label.x.npc = 0.7, label.y.npc = 0.1, size = 6, 
           cor.coef.name = 'r', method = 'pearson') +
  geom_smooth(method = 'lm')
p_diff_vs_tetInt

p_ddGs <- plot_grid(p_diff_vs_stab, p_diff_vs_dimInt, p_diff_vs_tetInt, ncol = 3)
# ggsave(p_ddGs, device = cairo_pdf, width = 27, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/Fig_ddGs_noTMP.pdf')

# Predictions for the model with all variables 
pred_rf_all <- read_delim('Results_ML/diffNorm_ara0.2_ara0.01/pred_rf_allVariables.txt',
                          delim =  '\t')

# Draw the figures
p_pred_rf_all <- pred_rf_all %>% ggplot(aes(x = pred_data, y = test_data)) +
  geom_point(size = 3) +
  theme(axis.title = element_text(size = 28, face = 'bold'), 
        axis.text = element_text(size = 26), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  geom_smooth(method = 'lm') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.05, size = 9,
           label.y.npc = 0.9, method = 'pearson', cor.coef.name = 'r'
  ) +
  labs(x = expression(paste(bold('Predicted \u0394'), bolditalic('s'), 
                            bold(' (RF, all variables)'), sep = '')), 
       y = expression(paste(bold('Observed \u0394'), bolditalic('s'), 
                            bold(' (test set)'), sep = ''))
  ) 
  # xlab('Predicted \u0394s (RF, all variables)') + ylab('Observed \u0394s (test set)')
p_pred_rf_all

# Relative importances for the model with all variables (permutations)
rel_importance_perm_all <- read_delim('Results_ML/diffNorm_ara0.2_ara0.01/model_diffFit_permImportances_allVariables.txt',
                                      delim =  '\t') %>% arrange(desc(Importance)) %>%
  mutate(Feature = ifelse(Feature == 'random_var', 'Random variable', 
                          ifelse(Feature == 'Mean_ddG_int_HM_A_D', '\u0394\u0394G tet. interface', 
                                 ifelse(Feature == 'Mean_ddG_int_HM_A_C', '\u0394\u0394G dim. interface',
                                        ifelse(Feature == 'parallel_beta_strand', 'Parallel beta strand', 
                                               ifelse(Feature == 'hydrophilicity_hopp', 'Hydrophilicity Hopp', 
                                                      ifelse(Feature == 'Mean_ddG_stab_HET', '\u0394\u0394G subunit stability', 
                                                             ifelse(Feature == 'rSASA', 'Rel. solvent accessibility', 
                                                                    ifelse(Feature == 'polarity_zimmerman', 'Polarity Zimmerman', 
                                                                           Feature)))))))))

rel_importance_perm_all %<>% mutate(Feature = factor(Feature, levels = Feature))


p_rel_importance <- rel_importance_perm_all %>% 
  mutate(color_check = (Feature == 'Random variable')) %>%
  ggplot(aes(y = Importance, x = Feature, fill = color_check)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c('grey', 'red')) +
  theme(axis.title = element_text(size = 24, face = 'bold'), 
        axis.text.x = element_text(size = 30, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 28),
        legend.position = 'none', 
        aspect.ratio = 1/5) +
  xlab('')
p_rel_importance

# # Show the figures together
# p_suppl_ml_bottom <- plot_grid(p_pred_rf_all, p_diff_ara0.2_ara0.01, p_diff_vs_stab, 
#                                labels = c('B', 'C', 'D'), ncol = 3, label_size = 20, label_fontface = 'bold')
# p_suppl_ml_bottom
# 
# p_supp_ml <- plot_grid(p_rel_importance, p_suppl_ml_bottom, labels = c('A', ''), nrow = 2,
#                        label_size = 20, label_fontface = 'bold')
# p_supp_ml

# Figure for predictions with the best variables
pred_rf_best <- read_delim('Results_ML/diffNorm_ara0.2_ara0.01/pred_rf_bestVariables.txt',
                          delim =  '\t')

p_pred_rf_best <- pred_rf_best %>% ggplot(aes(x = pred_data, y = test_data)) +
  geom_point(size = 3) +
  theme(axis.title = element_text(size = 28, face = 'bold'), 
        axis.text = element_text(size = 26), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  geom_smooth(method = 'lm') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.05, size = 9,
           label.y.npc = 0.9, method = 'pearson', cor.coef.name = 'r'
  ) +
  labs(x = expression(paste(bold('Predicted \u0394'), bolditalic('s'), 
                            bold(' (RF, best variables)'), sep = '')), 
       y = expression(paste(bold('Observed \u0394'), bolditalic('s'), 
                            bold(' (test set)'), sep = ''))
         ) 
  # xlab('Predicted \u0394s (RF, best variables)') + ylab('Observed \u0394s (test set)')
p_pred_rf_best

# Relative importance with the drop column method in the best model
rel_importance_dropCol_best <- read_delim('Results_ML/diffNorm_ara0.2_ara0.01/model_diffFit_dropCol_bestVariables.txt',
                                      delim =  '\t') %>% arrange(Importance) %>%
  mutate(Feature = ifelse(Feature == 'random_var', 'Random variable', 
                          ifelse(Feature == 'Mean_ddG_int_HM_A_D', '\u0394\u0394G tet. interface', 
                                 ifelse(Feature == 'Mean_ddG_int_HM_A_C', '\u0394\u0394G dim. interface',
                                        ifelse(Feature == 'parallel_beta_strand', 'Parallel beta strand', 
                                               ifelse(Feature == 'hydrophilicity_hopp', 'Hydrophilicity Hopp', 
                                                      ifelse(Feature == 'Mean_ddG_stab_HET', '\u0394\u0394G subunit stability', 
                                                             ifelse(Feature == 'rSASA', 'Rel. solvent accessibility', 
                                                                    ifelse(Feature == 'polarity_zimmerman', 'Polarity Zimmerman', 
                                                                           Feature)))))))))

rel_importance_dropCol_best %<>% mutate(Feature = factor(Feature, levels = Feature))


p_rel_importance_dropCol <- rel_importance_dropCol_best %>% 
  mutate(color_check = (Feature == 'Random variable')) %>%
  ggplot(aes(x = Importance, y = Feature, fill = color_check)) +
  geom_bar(stat = 'identity') +
  geom_vline(xintercept = 0) +
  scale_fill_manual(values = c('grey', 'red')) +
  theme(axis.title = element_text(size = 28, face = 'bold'), 
        axis.text.x = element_text(size = 26, angle = 0),
        axis.text.y = element_text(size = 26),
        legend.position = 'none')
p_rel_importance_dropCol

p_suppl_ml_bottom <- plot_grid(p_pred_rf_all, p_pred_rf_best,
                               p_rel_importance_dropCol, NULL,
                               ncol = 4,
                            labels = c('B', 'C', 'D', ''), label_size = 30, label_fontface = 'bold', 
                            rel_widths = c(1, 1, 0.8, 0.05))

p_figS16 <- plot_grid(p_rel_importance, p_suppl_ml_bottom, nrow = 2, labels = c('A', ''), 
                        label_size = 30, label_fontface = 'bold', rel_heights = c(1.3, 0.7))
p_figS16

ggsave(p_figS16, device = cairo_pdf, width = 28, height = 24, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS16_RF_training.pdf')
ggsave(p_figS16, device = 'png', width = 28, height = 24, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS16_RF_training.png')

#### Show results of the RF on the data without TMP ####

# Predictions for the model with all variables 
pred_rf_all <- read_delim('Results_ML/diffNorm_ara0.2_ara0.01/pred_rf_allVariables_notmp.txt',
                          delim =  '\t')

# Draw the figures
p_pred_rf_all <- pred_rf_all %>% ggplot(aes(x = pred_data, y = test_data)) +
  geom_point(size = 3) +
  theme(axis.title = element_text(size = 19, face = 'bold'), 
        axis.text = element_text(size = 18), 
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank()) +
  geom_smooth(method = 'lm') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = 0.7, size = 6,
           label.y.npc = 0.15, method = 'pearson', cor.coef.name = 'r'
  ) +
  xlab('Predicted \u0394s (RF, all variables, no TMP)') + ylab('Observed \u0394s (test set, no TMP)')
p_pred_rf_all

# Relative importances for the model with all variables (permutations)
rel_importance_perm_all <- read_delim('Results_ML/diffNorm_ara0.2_ara0.01/model_diffFit_permImportances_allVariables_notmp.txt',
                                      delim =  '\t') %>% arrange(desc(Importance)) %>%
  mutate(Feature = ifelse(Feature == 'random_var', 'Random variable', 
                          ifelse(Feature == 'Mean_ddG_int_HM_A_D', 'ddG tet. interface', 
                                 ifelse(Feature == 'Mean_ddG_int_HM_A_C', 'ddG dim. interface',
                                        ifelse(Feature == 'parallel_beta_strand', 'Parallel beta strand', 
                                               ifelse(Feature == 'hydrophilicity_hopp', 'Hydrophilicity Hopp', 
                                                      ifelse(Feature == 'Mean_ddG_stab_HET', 'ddG subunit stability', 
                                                             ifelse(Feature == 'rSASA', 'Rel. solvent accessibility', 
                                                                    ifelse(Feature == 'polarity_zimmerman', 'Polarity Zimmerman', 
                                                                           Feature)))))))))

rel_importance_perm_all %<>% mutate(Feature = factor(Feature, levels = Feature))


p_rel_importance <- rel_importance_perm_all %>% 
  mutate(color_check = (Feature == 'Random variable')) %>%
  ggplot(aes(y = Importance, x = Feature, fill = color_check)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = c('grey', 'red')) +
  theme(axis.title = element_text(size = 20, face = 'bold'), 
        axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 16),
        legend.position = 'none')
p_rel_importance

p_suppl_ml_notmp <- plot_grid(p_pred_rf_all, p_rel_importance, ncol = 2, labels = c('A', 'B'), 
                              label_size = 20, label_fontface = 'bold', 
                              rel_widths = c(0.5, 1))
ggsave(p_suppl_ml_notmp, device = cairo_pdf, width = 30, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS15_RF_training_notmp.pdf')
ggsave(p_suppl_ml_notmp, device = 'png', width = 30, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS15_RF_training_notmp.png')


# #### New panels on the cumulative distribution of significant mutants based on effects
# #### on stability and binding affinity
# data_compensated <- score_significant %>% rowwise() %>%
#   mutate(signif_check = ifelse(and(diffNormScore > 0.1, p.adj <= 0.05), 'Better at low expression', 
#                                ifelse(and(diffNormScore < -0.1, p.adj <= 0.05), 'Better at high expression', 
#                                       'Not significant')))
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
#   xlab('ddG on subunit stability') + ylab('Cumulative proportion of mutants') +
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
#   xlab('ddG on binding affinity (interface 1)') + ylab('Cumulative proportion of mutants') +
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
#   xlab('ddG on binding affinity (interface 2)') + ylab('Cumulative proportion of mutants') +
#   labs(colour = '')
# p_int2
# 
# # Repeat with a column of the max effect
# p_max <- data_compensated %>% filter(!(is.na(Mean_ddG_int_HM_A_D))) %>% rowwise() %>%
#   mutate(max_ddg = max(Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D)) %>%
#   ggplot(aes(x = max_ddg, colour = signif_check)) + 
#   stat_ecdf(size = 1.5) +
#   theme(axis.text.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18, angle = 90),
#         axis.title = element_text(size = 20, face = 'bold'), 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'none') +
#   xlab('Maximum ddG (stability or binding affinity)') + ylab('Cumulative proportion of mutants') +
#   labs(colour = '')
# p_max
# 
# p_cumulative_dist <- plot_grid(p_stab, p_int1, p_int2, p_max, nrow = 4)
# ggsave(p_cumulative_dist, device = cairo_pdf, width = 10, height = 28, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_ddGs_cumulative_dist.pdf')
# 
# # Similar figure but with the sum of ddGs
# p_sum <- data_compensated %>% filter(!(is.na(Mean_ddG_int_HM_A_D))) %>% rowwise() %>%
#   mutate(sum_ddg = sum(Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D)) %>%
#   ggplot(aes(x = sum_ddg, colour = signif_check)) + 
#   stat_ecdf(size = 1.5) +
#   theme(axis.text.y = element_text(size = 18),
#         axis.text.x = element_text(size = 18, angle = 90),
#         axis.title = element_text(size = 20, face = 'bold'), 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'none') +
#   xlab('Sum of ddGs (stability or binding affinity)') + ylab('Cumulative proportion of mutants') +
#   labs(colour = '')
# p_sum
# ggsave(p_sum, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-04-19_extra_figures/Fig_sum_ddGs_cumulative_dist.pdf')

#### Figure S17: New data with stop codon at first position, effect of arabinose ####

plate.ind <- 'GrowthCurves_Stop_2022-02-23/Fitness_negctl_ID_IGA_03_03_22.xlsx'
file.od <- 'GrowthCurves_Stop_2022-02-23/Growth_curves_03_03_2022.xlsx'

read.my.gc <- function(file, plate.index){
  pl <- read.xlsx(file,sheetIndex = 1, rowIndex = 4:99, stringsAsfators = FALSE,
                  header = F)
  # ind <- read.csv(plate.index, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  # ind <- read.csv2(plate.index, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  ind <- read.xlsx(plate.index, sheetIndex = 1, rowIndex = 1:97, header = T) 
  
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
  
  ## For the WT controls I will trim at t = 13.5
  # gc_out <- SummarizeGrowthByPlate(d, t_trim =13.5)
  gc_out <- SummarizeGrowthByPlate(d)
  colnames(gc_out)[1] <- "Well"
  data.pl %<>% left_join(gc_out, by = "Well")
  
}


data.od1 <- read.my.gc(file.od, plate.ind)

# Subtract the blank for this experiment and multiply by 5 to make it OD / mL
# (experiment was carried out in 0.2 mL)
data.od1 %<>% mutate(OD = (OD - 0.088) * 5)

# Remove wells on the border of the plate (rows A and H, columns 1 and 12)
data.od1 %<>% separate(Well, into = c('Well_row', 'Well_col'), sep = 1) %>%
  filter(!(Well_row %in% c('A', 'H')), !(Well_col %in% c(1, 12)))

# Show the growth curves
p_figs17 <- data.od1 %>% ungroup() %>% rowwise() %>%
  mutate(TMP = str_c(TMP, ' \u00B5g/mL', sep = ''), 
         Mutant = ifelse(Mutant == 'deltaMET', '\u0394MET', Mutant)) %>%
  mutate(Mutant = factor(Mutant, levels = c('\u0394MET', 'WT'))) %>%
  mutate(TMP = as.factor(TMP), Arabinose = str_c(Arabinose, '% arabinose')
  ) %>%
  # filter(Mutant == 'WT') %>%
  mutate(Arabinose = factor(Arabinose, levels = c('0% arabinose', 
                                                  '0.01% arabinose',
                                                  '0.025% arabinose', 
                                                  '0.05% arabinose',
                                                  '0.2% arabinose', 
                                                  '0.4% arabinose'))) %>%
  ggplot(aes(x = time, y = OD, colour = Mutant, group = interaction(Mutant, Replicate, TMP, Arabinose))) + 
  facet_wrap(~Arabinose, nrow = 2, scales = 'free') +
  geom_line() +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(face = 'bold', size = 28),
        axis.text = element_text(size = 24),
        axis.line = element_line(),
        legend.position = 'top', legend.title = element_text(size = 20), 
        strip.text = element_text(size = 20, face = 'bold'),
        strip.background = element_rect(fill = 'white'),
        legend.text = element_text(size = 18), legend.justification = 0.5) +
  guides(size = 'none', linetype = 'none', alpha = 'none') +
  xlab('Time (h)') + ylim(-0.1, 6) +
  labs(y = expression(bolditalic(OD[600]))) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
p_figs17

# # Edit the legend to make the lines stand out more
# p_s17_legend <- get_legend(p_figs17 + geom_line(size = 3) +
#                             theme(legend.position = 'top',
#                                   legend.title = element_text(size = 24),
#                                   legend.text = element_text(size = 22)))
# p_figs17_final <- plot_grid(p_s17_legend, p_figs17 + theme(legend.position = 'none'), 
#                            nrow = 2, rel_heights = c(0.1, 1))


ggsave(# p_figs17_final,
       p_figs17, device = cairo_pdf, width = 21, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS17.Arabinose_effect.pdf')
ggsave(# p_figs17_final,
       p_figs17, device = 'png', width = 21, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS17.Arabinose_effect.png')

## Add a new panel for the difference in growth: the average cost of expression ##

data.od1.summary <- data.od1 %>%
  group_by(TMP, Arabinose, Replicate, Mutant) %>%
  summarise(auc = mean(auc_e))

# Get the cost of the highest expression level
data.od1.summary.wide <- data.od1.summary %>% ungroup() %>%
  filter(Arabinose == 0.2) %>%
  mutate(Replicate = ifelse(Replicate == 5, 4, Replicate)) %>% # Adjust replicate numbering to have matching samples
  pivot_wider(names_from = Mutant, values_from = auc) %>%
  mutate(auc_diff = WT - deltaMET)

max_cost <- median(data.od1.summary.wide$auc_diff)

data.od1.summary %<>% ungroup() %>% group_by(Mutant, Arabinose) %>%
  mutate(Replicate = row_number()) # Renumber replicates

# Pivot to calculate costs of expression
data.od1.summary.pivot <- data.od1.summary %>% ungroup() %>%
  pivot_wider(names_from = Mutant, values_from = auc) %>%
  mutate(auc_diff = WT - deltaMET, max_cost = max_cost)

# Show a figure of the costs of expression
p_cost <- data.od1.summary.pivot %>% 
  ggplot(aes(x = as.factor(Arabinose), y = auc_diff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab('Arabinose concentration (% m/v)') + ylab('Translation cost (AUC)') + 
  theme(# panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
    # panel.grid.minor = element_blank(), 
    axis.title = element_text(face = 'bold', size = 20),
    axis.text = element_text(size = 18),
    legend.position = 'top', legend.title = element_text(size = 20), 
    strip.text = element_text(size = 20, face = 'bold'),
    strip.background = element_rect(fill = 'white'),
    legend.text = element_text(size = 18), legend.justification = 0.5)
p_cost
ggsave(p_cost, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS16.Arabinose_effect_cost.pdf')


### Fig. S18: Heatmaps of correlations between replicates (no TMP) ####
all_data_all_reps_corr_TMP0 <- all_data_all_reps %>% rowwise() %>%
  mutate(Genotype = str_c(WT_Residue, Residue, Position)) %>%
  select(ID, Genotype, Timepoint, TMP, Arabinose, sel_coeff, Sequencer) %>% 
  group_by(Genotype) %>%
  filter(TMP == 0, Timepoint == 10)#  %>%
  # filter(ID != 'E2.BR2.010.NT.SM') # Replicate with low read counts

cor_matrix <- all_data_all_reps_corr_TMP0  %>% ungroup() %>%
  select(-Timepoint, -Arabinose, -TMP, -Sequencer) %>%
  pivot_wider(names_from = ID, values_from = sel_coeff) %>%
  select(-Genotype) %>% cor(method = 'spearman')

# Average correlation between replicates, removing correlations with self
cor_vector <- c()
for(i in cor_matrix){
  if(i != 1){
    cor_vector <- c(cor_vector, i)
  }
}
mean(cor_vector)

annotation_df <- all_data_all_reps_corr_TMP0 %>% ungroup() %>%
  select(-sel_coeff, -Genotype) %>%
  unique() %>%
  mutate(Arabinose = factor(Arabinose, levels = c(0.01, 0.025, 0.05, 0.2, 0.4))) %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Overexpressed',
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal',
                                                  'Optimal', 'Overexpressed')))

ha_corr <- HeatmapAnnotation(`Expression level` = annotation_df$exp_level,
                             show_annotation_name = F,
                             annotation_name_gp = gpar(fontface = 'bold', fontsize = 14),
                             annotation_name_side = 'left',
                             show_legend = TRUE,
                             annotation_legend_param = list(
                               `Expression level` = list(
                                 title_gp = gpar(fontsize = 22, fontface = 'bold'), 
                                 labels_gp = gpar(fontsize = 20)
                               )
                             ),
                             col = list(`Expression level` = c("Weak" = "#fed976", "Suboptimal" = "#fd8d3c", 
                                                               "Near-optimal" = "#bd0026", "Optimal" = "#80001a", 
                                                               'Overexpressed' = 'black')),
                             gp = gpar(col = "black")
)

ha_corr2 <- HeatmapAnnotation(`Expression level` = annotation_df$exp_level,
                              which = 'row',
                              show_annotation_name = F,
                              annotation_name_gp = gpar(fontface = 'bold', fontsize = 14),
                              # annotation_name_side = 'left',
                              show_legend = FALSE,
                              annotation_legend_param = list(
                                `Expression level` = list(
                                  title_gp = gpar(fontsize = 22, fontface = 'bold'), 
                                  labels_gp = gpar(fontsize = 20)
                                )
                              ),
                              col = list(`Expression level` = c("Weak" = "#fed976", "Suboptimal" = "#fd8d3c", 
                                                                "Near-optimal" = "#bd0026", "Optimal" = "#80001a", 
                                                                'Overexpressed' = 'black')),
                              gp = gpar(col = "black")
)

# Draw the heatmap
p_figS18 <- Heatmap(cor_matrix, cluster_columns = T, cluster_rows = T, 
             clustering_distance_rows = 'pearson',
             clustering_distance_columns = 'pearson',
             ## Color ramp for the old normalized scores
             col = colorRamp2(
               # breaks = seq(-1, 1, length.out = 7), 
               breaks = seq(0, 1, length.out = 7),
               colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#F4A582", "#D6604D")),
             show_column_names = T, row_names_side = 'left',
             width=unit(25, 'cm'), height = unit(25, 'cm'),
             border = T,
             row_title = 'Sample',
             row_title_gp = gpar(fontsize=22, fontface = 'bold'),
             row_names_rot = 0, 
             row_names_centered = T,
             row_names_gp = gpar(fontsize=18, fontface = 'bold'),
             column_title = 'Sample', 
             column_title_gp = gpar(fontsize=22, fontface = 'bold'),
             column_title_side = 'bottom',
             # column_title = '0.01% arabinose', column_title_side = 'top',
             column_names_gp = gpar(fontsize=18,fontface='bold'),
             top_annotation = ha_corr,
             left_annotation = ha_corr2,
             cell_fun = function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y,
                         gp = gpar(fontsize = 18, fontface = 'bold'))
             },
             show_heatmap_legend = TRUE,
             heatmap_legend_param = list(
               ## Ticks positions for the old normalized scores
               # at = c(-1, 0, 1, 2),
               ## Ticks positions for the new normalized scores
               # at = c(-4, -2, 0, 2, 4),
               ## Ticks for the scale from -8 to +4
               at = c(0, 0.5, 1),
               #labels = c("low", "zero", "high"),
               # title = "Spearman rho", 
               title = 'Spearman ⍴',
               title_gp = gpar(fontsize = 20),
               legend_height = unit(3.5, "cm"),
               legend_width = unit(2, "cm"),
               border='black',
               lwd=1.7,
               labels_gp = gpar(fontsize = 18),
               title_position = "leftcenter-rot"
             )
)
p_figS18
ggsave(grid.grabExpr(draw(p_figS18)), width = 20, height = 14, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS18.heatmap_corr_samples_noTMP.pdf')
ggsave(grid.grabExpr(draw(p_figS18)), width = 20, height = 14, dpi = 300, device = 'png', 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS18.heatmap_corr_samples_noTMP.png')

#### Fig. S19: General boxplots of DMS s values with and without TMP ####

#### Draw boxplots of the selection coefficient distributions with and without TMP side by side ####
p_figs19 <- all_data_complete %>% filter(Timepoint == 10) %>%
  mutate(Arabinose = factor(Arabinose, levels = c(0.01, 0.025, 0.05, 0.2, 0.4)), 
         TMP = factor(TMP, levels = c(0, 10)),
         exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal',
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal',
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal',
                                                  'Near-optimal', 'Optimal', 
                                                  'Overexpressed'))) %>%
  ggplot(aes(x = TMP, y = mean_sel_coeff, fill = exp_level)) +
  geom_point(alpha = 1, position = position_jitterdodge(jitter.width = 0.2),
             aes(colour = exp_level)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(), 
               alpha = 0.5) +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'))+
  theme(axis.title = element_text(face = 'bold', size = 20),
        axis.text = element_text(size = 18),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'top', 
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18), 
        legend.justification = 'center') +
  # ylab('s') +
  xlab('TMP (\u00B5g/mL)') +
  labs(fill = 'Expression level', colour = 'Expression level',
       y = expression(bolditalic('s')))
p_figs19
ggsave(p_figs19, width = 14, height = 10, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS19.boxplots_whole_dist_TMP_noTMP.pdf')
ggsave(p_figs19, width = 14, height = 10, dpi = 300, device = 'png', 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS19.boxplots_whole_dist_TMP_noTMP.png')

# Similar boxplots but per position
p <- all_data_complete %>% rowwise() %>% filter(Timepoint == 10) %>%
  mutate(TMP = str_c('TMP = ', TMP, ' \u00B5g/mL', sep = ''),
         exp_level = ifelse(Arabinose == 0.01, 'Weak expression', 
                            ifelse(Arabinose == 0.025, 'Suboptimal expression',
                                   ifelse(Arabinose == 0.05, 'Near-optimal expression', 
                                          ifelse(Arabinose == 0.2, 'Optimal expression',
                                                 ifelse(Arabinose == 0.4, 'Overexpression',
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak expression', 'Suboptimal expression',
                                                  'Near-optimal expression', 'Optimal expression', 
                                                  'Overexpression'))) %>%
  ggplot(aes(x = Position, y = mean_sel_coeff, group = Position)) +
  facet_grid(exp_level~TMP) +
  geom_point(alpha = 1, position = position_jitterdodge(jitter.width = 0.2),
             aes(colour = exp_level)) +
  geom_boxplot(aes(fill = exp_level), outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black')) +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'))+
  theme(axis.title = element_text(face = 'bold', size = 20),
        axis.text = element_text(size = 18),
        strip.text = element_text(size = 18, face = 'bold'),
        strip.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        legend.position = 'none') +
  ylab('Sel. coeff.') + xlab('Position')
p
ggsave(p, width = 14, height = 24, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigSX.boxplots_whole_distPosition_TMP_noTMP.pdf')
ggsave(p, width = 14, height = 24, dpi = 300, device = 'png', 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigSX.boxplots_whole_distPosition_TMP_noTMP.png')

#### Fig. S20: stop codons ####

# all_data_stop_check <- all_data_complete %>% 
#   select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP,
#          mean_sel_coeff) %>% 
#   mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
#                             ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
#   filter(Timepoint == 10) %>%
#   mutate(TMP = str_c('TMP = ', TMP, ' \u00B5g/mL', sep = ''), 
#          Timepoint = str_c('Time = ', Timepoint, ' generations'), 
#          exp_level = ifelse(Arabinose == 0.01, 'Weak expression', 
#                             ifelse(Arabinose == 0.025, 'Suboptimal expression', 
#                                    ifelse(Arabinose == 0.05, 'Near-optimal expression', 
#                                           ifelse(Arabinose == 0.2, 'Optimal expression',
#                                                  ifelse(Arabinose == 0.4, 'Overexpression', 
#                                                         NA)))))
#   ) %>%
#   mutate(exp_level = factor(exp_level, levels = c('Weak expression', 'Suboptimal expression',
#                                                   'Near-optimal expression', 'Optimal expression', 
#                                                   'Overexpression'))) %>%
#   mutate(mut_check = factor(mut_check, levels = c('Stop', 'Missense', 'WT')))
# 
# # Calculate p values for comparisons
# comps <- compare_means(mean_sel_coeff~mut_check, 
#                        data = all_data_stop_check %>% 
#                          select(-Position, -WT_Residue, -Residue, -Timepoint),
#                        paired = F, group.by = c('TMP', 'exp_level'), 
#                        method = 'wilcox.test') %>%
#   mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
#                            ifelse(p > 0.01, str_c('p = ', round(as.numeric(p), 2), sep = ''),
#                            sprintf("p = %2.1e", as.numeric(p)))
#   ),
#   y_pos = rep(c(0.6, 0.7, 0.8), 10)
#   )
# 
# # Let's try to do a figure of the distributions of log2FC for stop and WT codons
# p_figs19 <- all_data_stop_check %>% rowwise() %>%
#   ggplot(aes(x = mut_check,y = mean_sel_coeff)) + 
#   geom_point(aes(colour = mut_check), position = position_jitterdodge(), show.legend = F) +
#   geom_boxplot(aes(fill = mut_check), colour = 'black', outlier.shape = NA, alpha = 0.5) + 
#   facet_grid(TMP~exp_level) + 
#   geom_signif(data = as.data.frame(comps), inherit.aes = FALSE,
#               aes(xmin = group1, xmax = group2, annotations=p.format, y_position = y_pos), 
#               manual = TRUE, textsize = 5) +
#   geom_hline(yintercept = 0, linetype = 'dashed') +
#   xlab('Mutation type') +
#   labs(fill = 'Mutant type') +
#   ylab('s') +
#   theme(axis.title = element_text(face = 'bold', size = 20), 
#         axis.text = element_text(size = 18), 
#         strip.text = element_text(size = 20, face = 'bold'), 
#         strip.background = element_rect(fill = 'white'))
# p_figs19
# ggsave(p_figs19, width = 21, height = 14, device = cairo_pdf, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS19_Stop_vs_WT.pdf')
# ggsave(p_figs19, width = 21, height = 14, device = 'png', dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS19_Stop_vs_WT.png')

#### Repeat the figure at the codon level ####

# Load the complete dataset
all_data_complete_codons <- read_delim(
  '../R67_DMS_February2022/Complete_datasets_MiSeq_NovaSeq/complete_dataset_avgBothSequencers_Codons.txt', delim = '\t')

all_data_stop_check_codons <- all_data_complete_codons %>% 
  select(Position, WT_Codon, Codon, WT_Residue, Encoded_residues, Timepoint, Arabinose, TMP,
         mean_sel_coeff) %>% 
  filter(Codon != WT_Codon, Codon != 'TAG') %>%
  mutate(mut_check = ifelse(Encoded_residues == WT_Residue, 'Synonymous',
                            ifelse(Codon %in% c('TAA', 'TGA'), 'Stop',
                                   # 'Amino acid\nsubstitution'))) %>%
                                   'AA\nsubstitution'))) %>%
  filter(Timepoint == 10) %>%
  mutate(TMP = str_c('TMP = ', TMP, ' \u00B5g/mL', sep = ''), 
         Timepoint = str_c('Time = ', Timepoint, ' generations'), 
         exp_level = ifelse(Arabinose == 0.01, 'Weak expression', 
                            ifelse(Arabinose == 0.025, 'Suboptimal expression', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal expression', 
                                          ifelse(Arabinose == 0.2, 'Optimal expression',
                                                 ifelse(Arabinose == 0.4, 'Overexpression', 
                                                        NA)))))
  ) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak expression', 'Suboptimal expression',
                                                  'Near-optimal expression', 'Optimal expression', 
                                                  'Overexpression'))) %>%
  mutate(mut_check = factor(mut_check, levels = c('Stop',
                                                  # 'Amino acid\nsubstitution',
                                                  'AA\nsubstitution',
                                                  'Synonymous')))

# Calculate p values for comparisons
comps <- compare_means(mean_sel_coeff~mut_check, 
                       data = all_data_stop_check_codons %>% 
                         select(-Position, -WT_Codon, -Codon, -Timepoint),
                       paired = F, group.by = c('TMP', 'exp_level'), 
                       method = 'wilcox.test') %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
                           ifelse(p > 0.01, str_c('p = ', round(as.numeric(p), 2), sep = ''),
                                  sprintf("p = %2.1e", as.numeric(p)))
  ),
  y_pos = rep(c(0.55, 0.75, 0.95), 10)
  )

# Let's try to do a figure of the distributions of log2FC for stop and WT codons
p_figs20 <- all_data_stop_check_codons %>% rowwise() %>%
  ggplot(aes(x = mut_check,y = mean_sel_coeff)) + 
  geom_point(aes(colour = mut_check), position = position_jitterdodge(), show.legend = F) +
  geom_boxplot(aes(fill = mut_check), colour = 'black', outlier.shape = NA, alpha = 0.5) + 
  facet_grid(TMP~exp_level, scales = 'free') + 
  geom_signif(data = as.data.frame(comps), inherit.aes = FALSE,
              aes(xmin = group1, xmax = group2, annotations=p.format, y_position = y_pos), 
              manual = TRUE, textsize = 10) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  xlab('Mutation type') +
  labs(fill = 'Mutant type', y = expression(bolditalic(s))) +
  theme(axis.title = element_text(face = 'bold', size = 30), 
        axis.text = element_text(size = 28), 
        strip.text = element_text(size = 30, face = 'bold'), 
        strip.background = element_rect(fill = 'white'), 
        legend.position = 'none',
        axis.line = element_line()) +
  ylim(-0.9, 1.1)
p_figs20
ggsave(p_figs20, width = 32, height = 14, device = cairo_pdf, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS20_Stop_vs_WT_codons.pdf')
ggsave(p_figs20, width = 32, height = 14, device = 'png', dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS20_Stop_vs_WT_codons.png')

#### Look at differences in codon usage for codons synonymous to WT and selection coefficients ####

# Load data on E. coli codon usage
codon_usage_ecoli <- read_delim('Data/Codon_usage/E.coli/Codon_usage_E.coli_clean.txt', delim = '\t') %>%
  rowwise() %>%
  mutate(DNA_codon = str_replace_all(string = Codon, pattern = 'U', replacement = 'T'))

# Add the usage data to the selection coefficients of codons synonymous to WT
all_data_complete_codons_usage <- left_join(x = all_data_complete_codons %>%
                                              filter(WT_Residue == Encoded_residues, WT_Codon != Codon), 
                                            y = codon_usage_ecoli %>% select(-Codon),
                                            by = c('Codon' = 'DNA_codon'))

# Check if codon usage itself correlates to selection coefficient
p <- all_data_complete_codons_usage %>% ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = Percentage_synonymous)) +
  geom_point() +
  facet_wrap(~Arabinose) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/1.Codon_usage_sel_coeff.pdf')

### Check if difference in codon usage between the WT and the synonymous codon correlates to selection coefficient
all_data_complete_codons_usage <- left_join(x = all_data_complete_codons %>%
                                              filter(WT_Residue == Encoded_residues), 
                                            y = codon_usage_ecoli %>% select(-Codon),
                                            by = c('Codon' = 'DNA_codon'))

all_data_complete_codons_usage_wt <-all_data_complete_codons_usage %>%
  filter(Codon == WT_Codon)

all_data_complete_codons_usage_syn <- all_data_complete_codons_usage %>%
  filter(Codon != WT_Codon)

all_data_complete_codons_usage_new <- left_join(x = all_data_complete_codons_usage_syn, 
                                                y = all_data_complete_codons_usage_wt %>% 
                                                  mutate(Percentage_synonymous_wt = Percentage_synonymous) %>%
                                                  select(Position, Timepoint, Arabinose, TMP, Percentage_synonymous_wt), 
                                                by = c('Position' = 'Position', 'Timepoint' = 'Timepoint', 
                                                       'Arabinose' = 'Arabinose', 'TMP' = 'TMP'))

# Calculate the usage difference between the synonymous codon and the WT
all_data_complete_codons_usage_new %<>% mutate(diffUsage = Percentage_synonymous - Percentage_synonymous_wt)

p <- all_data_complete_codons_usage_new %>% ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = diffUsage)) +
  facet_wrap(~Arabinose) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/2.Codon_usage_diff_sel_coeff.pdf')

p <- all_data_complete_codons_usage_new %>% 
  mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt) %>%
  ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = rel_Usage)) +
  facet_wrap(~Arabinose) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/3.Rel_codon_usage_sel_coeff.pdf')

## Repeat the above figure indicating the position of the mismatch
codons_usage_split <- all_data_complete_codons_usage_new %>% rowwise() %>%
  separate(col = WT_Codon, into = c('WT_Pos1', 'WT_Pos2', 'WT_Pos3'), sep = c(1,2), remove = FALSE) %>%
  separate(col = Codon, into = c('Mut_Pos1', 'Mut_Pos2', 'Mut_Pos3'), sep = c(1,2), remove = FALSE) %>%
  mutate(Mismatch1 = (WT_Pos1 != Mut_Pos1), 
         Mismatch2 = (WT_Pos2 != Mut_Pos2), 
         Mismatch3 = (WT_Pos3 != Mut_Pos3)) %>%
  mutate(Mismatch_total = Mismatch1 + Mismatch2 + Mismatch3)

p <- codons_usage_split %>% 
  mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt) %>%
  ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = rel_Usage, colour = as.factor(Mismatch_total))) +
  facet_wrap(~Arabinose) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  # stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
  #          label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/3.Rel_codon_usage_sel_coeff_colorMismatches.pdf')

## Try showing GC content of codons
codons_GC_content_split <- all_data_complete_codons_usage_new %>% rowwise() %>%
  separate(col = WT_Codon, into = c('WT_Pos1', 'WT_Pos2', 'WT_Pos3'), sep = c(1,2), remove = FALSE) %>%
  separate(col = Codon, into = c('Mut_Pos1', 'Mut_Pos2', 'Mut_Pos3'), sep = c(1,2), remove = FALSE) %>%
  rowwise() %>%
  mutate(GC_content_mutant = sum(
    Mut_Pos1 %in% c('C', 'G'), Mut_Pos2 %in% c('C', 'G'), Mut_Pos3 %in% c('C', 'G')
  ), 
  GC_content_WT = sum(
    WT_Pos1 %in% c('C', 'G'), WT_Pos2 %in% c('C', 'G'), WT_Pos3 %in% c('C', 'G')
  )) %>%
  mutate(GC_diff = GC_content_mutant - GC_content_WT)

# Redraw the figure with just GC content of the mutant
p <- codons_GC_content_split %>% 
  mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt) %>%
  ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = rel_Usage, colour = as.factor(GC_content_mutant))) +
  facet_wrap(~Arabinose) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  # stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
  #          label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/3.Rel_codon_usage_sel_coeff_colorGCMutant.pdf')

## Check the effect of the difference in GC content
p <- codons_GC_content_split %>% 
  mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt) %>%
  ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = rel_Usage, colour = as.factor(GC_diff))) +
  facet_wrap(~Arabinose) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  # stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
  #          label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/3.Rel_codon_usage_sel_coeff_colorMismatches.pdf')



# Show the codons with the highest and lowest relative usage per position
max_rel_usage <- all_data_complete_codons_usage_new %>% ungroup() %>%
  group_by(Position, WT_Codon, Encoded_residues, Timepoint, Arabinose, TMP) %>%
  mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt) %>%
  summarise(num_codons = n(), 
            max_rel_Usage = max(rel_Usage)) %>%
  # Filter for residues with at least 2 synonymous codons other than the WT
  filter(num_codons >= 2)

min_rel_usage <- all_data_complete_codons_usage_new %>% ungroup() %>%
  group_by(Position, WT_Codon, Encoded_residues, Timepoint, Arabinose, TMP) %>%
  mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt) %>%
  summarise(num_codons = n(), 
            min_rel_Usage = min(rel_Usage)) %>%
  # Filter for residues with at least 2 synonymous codons other than the WT
  filter(num_codons >= 2)

## Alternative ways to filter for max and min usage ##
# max_rel_usage <- all_data_complete_codons_usage_new %>% ungroup() %>%
#   group_by(Position, WT_Codon, Encoded_residues, Timepoint, Arabinose, TMP) %>%
#   mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt) %>%
#   filter()
#   summarise(num_codons = n(), 
#             max_rel_Usage = max(rel_Usage)) %>%
#   # Filter for residues with at least 2 synonymous codons other than the WT
#   filter(num_codons >= 2)
# 
# min_rel_usage <- all_data_complete_codons_usage_new %>% ungroup() %>%
#   group_by(Position, WT_Codon, Encoded_residues, Timepoint, Arabinose, TMP) %>%
#   mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt) %>%
#   summarise(num_codons = n(), 
#             min_rel_Usage = min(rel_Usage)) %>%
#   # Filter for residues with at least 2 synonymous codons other than the WT
#   filter(num_codons >= 2)

# Concatenate minimum usage and maximum usage
usage_checks <- bind_rows(max_rel_usage %>% rowwise() %>% mutate(usage_check = max_rel_Usage, usage_bool = 'max') %>%
                            select(-max_rel_Usage), 
                          min_rel_usage %>% rowwise() %>% mutate(usage_check = min_rel_Usage, usage_bool = 'min') %>%
                            select(-min_rel_Usage))

# Use an inner join to only keep the codons from the usage checks (maximum and minimum)
filtered_codon_usage <- inner_join(x = all_data_complete_codons_usage_new %>%
                                     mutate(rel_Usage = Percentage_synonymous / Percentage_synonymous_wt), 
                                   y = usage_checks, 
                                   by = c('Position' = 'Position', 'WT_Codon' = 'WT_Codon', 'Encoded_residues' = 'Encoded_residues', 
                                          'Timepoint' = 'Timepoint', 'Arabinose' = 'Arabinose', 'TMP' = 'TMP', 
                                          'rel_Usage' = 'usage_check'))

# Try to get the medians of maximum vs minimum
medians <- filtered_codon_usage %>% filter(TMP == 10) %>%
  group_by(Arabinose, usage_bool) %>% summarise(median_check = median(mean_sel_coeff)) 

p <- filtered_codon_usage %>% filter(TMP == 10) %>%
  ggplot(aes(x = as.factor(Position), colour = usage_bool)) +
  geom_point(aes(y = mean_sel_coeff), size = 2) +
  geom_hline(aes(colour = usage_bool, yintercept = median_check), data = medians, 
             linetype = 'dashed') +
  facet_wrap(~Arabinose, nrow = 5) +
  theme(legend.position = 'top', 
        panel.grid.major.x = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(p, device = cairo_pdf, width = 10, height = 24, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/7.MaxMin_Rel_codon_usage_sel_coeff.pdf')

p <- filtered_codon_usage %>% filter(TMP == 10) %>%
  ggplot(aes(x = as.factor(Position), colour = usage_bool)) +
  geom_point(aes(y = rel_Usage), size = 2) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  theme(legend.position = 'top', 
        panel.grid.major.x = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5))
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/8.MaxMin_Rel_codon_usage_usageComp.pdf')

comps <- compare_means(mean_sel_coeff~usage_bool,
                       # Positions 12, 30, 71 (all valines) have two codons with identical minimum usage
                       # but different selection coefficients. I will remove them here
                       data = filtered_codon_usage %>% filter(!(Position %in% c(12, 30, 71))),
                       group.by = 'Arabinose',
                       paired = TRUE, method = 'wilcox.test')%>%
  mutate(
    # p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16', sprintf("p = %2.1e", as.numeric(p))),
  y_pos = c(0.3)
  )

# Try boxplots to compare the distributions
p <- filtered_codon_usage %>% filter(TMP == 10) %>%
  filter(!(Position %in% c(12, 30, 71))) %>%
  ggplot(aes(x = usage_bool, y = mean_sel_coeff, colour = usage_bool, fill = usage_bool)) +
  geom_point() +
  geom_line(aes(group = interaction(Arabinose, Position)), colour = 'black') +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_signif(data = as.data.frame(comps), inherit.aes = FALSE,
              aes(xmin = group1, xmax = group2, annotations=p.format, y_position = y_pos), 
              manual = TRUE) +
  facet_wrap(~Arabinose, ncol = 5)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/9.MaxMin_Rel_codon_boxplots.pdf')

filtered_codon_usage_new <- filtered_codon_usage %>% ungroup() %>%
  select(Position, mean_sel_coeff, Arabinose, TMP, usage_bool) %>%
  pivot_wider(names_from = usage_bool, values_from = mean_sel_coeff, 
              values_fn = function(in_list) return(in_list[[1]]))

p <- filtered_codon_usage_new %>% filter(TMP == 10) %>%
  ggplot(aes(x = min, y = max)) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F) +
  facet_wrap(~Arabinose) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  xlim(-0.3, 0.3) + ylim(-0.3, 0.3)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/10.MaxMin_Rel_codon_scatterplots.pdf')

# Check if we have more points above or under the diagonal
filtered_codon_usage_new %<>% mutate(diag_check = ifelse(max > min, 'above', 'below'))

test_check <- filtered_codon_usage_new %>%filter(TMP == 10)
table(test_check$diag_check, test_check$Arabinose)

#### Repeat with tAI ####

# Load tAi values
tAI_codon <- read_delim('Data/Codon_usage/tuller_2010_supp_tables_tAI.csv',
                        delim = '\t', locale = locale(decimal_mark = ','))
colnames(tAI_codon) <- c('Codon', 'Anticodon', 'S.cer', 'H.sap', 'E.coliK12', 'C.ele', 'D.mel', 'A.per')

# Add the usage data to the selection coefficients of codons synonymous to WT
all_data_complete_codons_tAI <- left_join(x = all_data_complete_codons %>%
                                              filter(WT_Residue == Encoded_residues, WT_Codon != Codon), 
                                            y = tAI_codon %>% select(Codon, E.coliK12),
                                            by = c('Codon' = 'Codon'))

# Check if codon usage itself correlates to selection coefficient
p <- all_data_complete_codons_tAI %>% ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = E.coliK12)) +
  facet_wrap(~Arabinose) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F) +
  xlab('tAI (E. coli K12)')
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/4.tAI_sel_coeff.pdf')

### Check if difference in codon usage between the WT and the synonymous codon correlates to selection coefficient
all_data_complete_codons_tAI <- left_join(x = all_data_complete_codons %>%
                                              filter(WT_Residue == Encoded_residues), 
                                          y = tAI_codon %>% select(Codon, E.coliK12),
                                          by = c('Codon' = 'Codon'))

all_data_complete_codons_tAI_wt <-all_data_complete_codons_tAI %>%
  filter(Codon == WT_Codon)

all_data_complete_codons_tAI_syn <- all_data_complete_codons_tAI %>%
  filter(Codon != WT_Codon)

all_data_complete_codons_tAI_new <- left_join(x = all_data_complete_codons_tAI_syn, 
                                                y = all_data_complete_codons_tAI_wt %>% 
                                                  mutate(tAI_wt = E.coliK12) %>%
                                                  select(Position, Timepoint, Arabinose, TMP, tAI_wt), 
                                                by = c('Position' = 'Position', 'Timepoint' = 'Timepoint', 
                                                       'Arabinose' = 'Arabinose', 'TMP' = 'TMP'))

# Calculate the usage difference between the synonymous codon and the WT
all_data_complete_codons_tAI_new %<>% mutate(difftAI = E.coliK12 - tAI_wt)

p <- all_data_complete_codons_tAI_new %>% ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = difftAI)) +
  facet_wrap(~Arabinose) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/5.difftAI_sel_coeff.pdf')

p <- all_data_complete_codons_tAI_new %>% 
  filter(TMP == 10) %>%
  mutate(rel_tAI = E.coliK12 / tAI_wt) %>%
  ungroup() %>% filter(TMP == 10) %>%
  ggplot(aes(y = mean_sel_coeff, x = log10(rel_tAI))) +
  facet_wrap(~Arabinose) +
  geom_point() +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'r', size = 4,
           label.x.npc = 0.05, label.y.npc = 0.8, method = 'pearson') +
  geom_smooth(method = 'lm', show.legend = F)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-05-19_Figures_codon_usage/6.Rel_tAI_sel_coeff.pdf')

#### Figure S21: GEMME vs DMS (without TMP) ####

p_figS21 <- gemme_vs_dms_plot %>% rowwise() %>% 
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak expression', 
                            ifelse(Arabinose == 0.025, 'Suboptimal expression', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal expression', 
                                          ifelse(Arabinose == 0.2, 'Optimal expression',
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak expression', 'Suboptimal expression', 
                                                  'Near-optimal expression', 'Optimal expression', 
                                                  'Overexpressed'))) %>%
  filter(TMP == 0) %>%
  mutate(Arabinose = str_c(toString(Arabinose), '% arabinose', sep = '')) %>%
  ggplot(aes(x = Fitness, y = mean_sel_coeff)) + 
  geom_point() + 
  xlab('GEMME score') + # ylab('s') +
  labs(y = expression(bolditalic(s))) +
  # facet_wrap(~Arabinose) +
  facet_wrap(~exp_level, scales = 'free', nrow = 3) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(face = 'bold', size = 20), axis.text = element_text(size = 18),
        strip.text = element_text(size = 20, face = 'bold'), 
        strip.background = element_rect(fill = 'white'), 
        axis.line = element_line()) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = 'rho', size = 8,
           label.x.npc = 0.05, label.y.npc = 0.9, method = 'spearman') +
  geom_smooth(method = 'lm', show.legend = F) +
  ylim(-0.1, 0.2)
p_figS21
ggsave(plot = p_figS21, device = cairo_pdf, width = 14, height = 17, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS21_GEMME_noTMP.pdf')
ggsave(plot = p_figS21, device = 'png', width = 14, height = 17, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS21_GEMME_noTMP.png')


#### Fig. S22: Destabilization of DfrB1 does not result in deleterious effects ####

all_data_complete_new <- all_data_complete %>% filter(TMP == 0, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D) %>%
  rowwise() %>%
  mutate(ID = str_c(Position, WT_Residue, Residue, '_', Arabinose, sep = '')) %>%
  mutate(ID = ifelse(ID == '2EE_0.01', 'WT_0.01', ID)) %>%
  mutate(ID = ifelse(ID == '2EE_0.2', 'WT_0.2', ID))

## Reload the data for the validated mutants
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
                          y = data.od1.new %>% filter(TMP == 0),
                          by = c('ID' = 'ID')) %>%
  mutate(Arabinose = as.factor(Arabinose))

joined_sets %<>% mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                                           ifelse(Arabinose == 0.025, 'Suboptimal', 
                                                  ifelse(Arabinose == 0.05, 'Near-optimal', 
                                                         ifelse(Arabinose == 0.2, 'Optimal',
                                                                ifelse(Arabinose == 0.4, 'Overexpressed', NA))))))

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

# Draw the figure for ddG stability
p_figs22_stab <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_stab))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_stab), size = 3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 6, fontface = 'bold') +
  # xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  # scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  scale_colour_manual(values = c('#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 26), 
    axis.title.y = element_text(face = 'bold', size = 32),
    axis.text = element_text(size = 24), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22)
  ) + 
  xlim(-0.05, 0.1) + ylim(-0.05, 0.05) +
  annotate('text', x = 0.07, y = 0.025, 
           label = expression(italic(s[weak] > s[opt])), parse = T, size = 10) +
  annotate('text', x = 0.07, y = -0.025, 
           label = expression(italic(s[weak] > s[opt])), parse = T, size = 10) +
  labs(colour = expression(paste(bold('\u0394\u0394'), bolditalic('G'),
                                 bold(' subunit stability [kcal / mol]'), sep = '')),
       y = expression(paste(bold('\u0394'), bolditalic(s[weak]), bold(' ('), 
                            bolditalic(s[weak]), bold(' - '), bolditalic(s[opt]),
                            bold(')'), sep = '')), 
       x = expression(paste(bolditalic('s'), bold(' (optimal expression)'), sep =  ''))) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))
p_figs22_stab

# Draw the figure for ddG dim int
p_figs22_dim_int <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_dim_int))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_dim_int), size = 3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 6, fontface = 'bold') +
  # xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  # scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  scale_colour_manual(values = c('#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  # labs(colour = '\u0394\u0394G dim. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 26), 
    axis.title.y = element_text(face = 'bold', size = 32),
    axis.text = element_text(size = 24), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22)
  ) + 
  xlim(-0.05, 0.1) + ylim(-0.05, 0.05) +
  # annotate('text', x = -0.03, y = 0.025, 
  #          label = 's[weak] > s[opt]', parse = T, size = 10) +
  # annotate('text', x = -0.03, y = -0.025, 
  #          label = 's[weak] < s[opt]', parse = T, size = 10) +
  annotate('text', x = 0.07, y = 0.025, 
           label = expression(italic(s[weak] > s[opt])), parse = T, size = 10) +
  annotate('text', x = 0.07, y = -0.025, 
           label = expression(italic(s[weak] > s[opt])), parse = T, size = 10) +
  labs(colour = expression(paste(bold('\u0394\u0394'), bolditalic('G'),
                                 bold(' dim. interface [kcal / mol]'), sep = '')),
       y = expression(paste(bold('\u0394'), bolditalic(s[weak]), bold(' ('), 
                            bolditalic(s[weak]), bold(' - '), bolditalic(s[opt]),
                            bold(')'), sep = '')), 
       x = expression(paste(bolditalic('s'), bold(' (optimal expression)'), sep =  ''))) +
  
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))

p_figs22_dim_int

# Draw the figure for ddG tet int
p_figs22_tet_int <- data_fig_4a %>% rowwise() %>%
  filter(!(is.na(bins_tet_int))) %>%
  ggplot(aes(x = mean_s, y = mean_ds)) +
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
  geom_point(aes(colour = bins_tet_int), size =3) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4, size = 6, fontface = 'bold') +
  # xlab('s (optimal expression)') + 
  # scale_size_manual(values = c(3, 5, 5)) +
  guides(size = 'none', alpha = 'none') +
  # scale_colour_viridis(option = 'A', discrete = TRUE) +
  # scale_colour_manual(values = c('#b3cde3', '#8c96c6', '#8856a7', '#810f7c')) +
  scale_colour_manual(values = c('#dadaeb', '#bcbddc', '#9e9ac8', '#756bb1', '#54278f')) +
  # labs(colour = '\u0394\u0394G tet. interface [kcal/mol]', y = expression(s[weak] - s[opt])) +
  theme(
    axis.title.x = element_text(face = 'bold', size = 26), 
    axis.title.y = element_text(face = 'bold', size = 32),
    axis.text = element_text(size = 24), 
    legend.position = 'top',
    legend.justification = 'center', 
    # legend.key.width = unit(1.15, 'cm'),
    legend.title = element_text(size = 24, face = 'bold'),
    legend.text = element_text(size = 22)
  ) + 
  xlim(-0.05, 0.1) + ylim(-0.05, 0.05) +
  # annotate('text', x = -0.03, y = 0.025, 
  #          label = 's[weak] > s[opt]', parse = T, size = 10) +
  # annotate('text', x = -0.03, y = -0.025, 
  #          label = 's[weak] < s[opt]', parse = T, size = 10) +
  annotate('text', x = 0.07, y = 0.025, 
           label = expression(italic(s[weak] > s[opt])), parse = T, size = 10) +
  annotate('text', x = 0.07, y = -0.025, 
           label = expression(italic(s[weak] > s[opt])), parse = T, size = 10) +
  labs(colour = expression(paste(bold('\u0394\u0394'), bolditalic('G'),
                                 bold(' tet. interface [kcal / mol]'), sep = '')),
       y = expression(paste(bold('\u0394'), bolditalic(s[weak]), bold(' ('), 
                            bolditalic(s[weak]), bold(' - '), bolditalic(s[opt]),
                            bold(')'), sep = '')), 
       x = expression(paste(bolditalic('s'), bold(' (optimal expression)'), sep =  ''))) +
  guides(colour = guide_legend(title.position = 'top', 
                               title.hjust = 0.5))
p_figs22_tet_int

#### Test the association between s and ddG ####

# # Use correlations
# test_corr_stab <- all_data_complete %>% ungroup() %>%
#   filter(!(is.na(Mean_ddG_stab_HET))) %>% rowwise() %>%
#   select(TMP, Arabinose, Position, WT_Residue, Residue, Mean_ddG_stab_HET, mean_sel_coeff) %>%
#   group_by(TMP, Arabinose)  %>%
#   summarise(
#     corr_stab_s = cor.test(Mean_ddG_stab_HET, mean_sel_coeff, method = 'spearman')$estimate, 
#     pval_corr_stab_s = cor.test(Mean_ddG_stab_HET, mean_sel_coeff, method = 'spearman')$p.value
#   )
# write.table(test_corr_stab, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
#             file = 'Figures/2022-05-09_Supp_figures_paper/Tables_corr_ddG/corr_ddG_stab.csv')
# 
# test_corr_dim_int <- all_data_complete %>% ungroup() %>%
#   filter(!(is.na(Mean_ddG_int_HM_A_C))) %>% rowwise() %>%
#   select(TMP, Arabinose, Position, WT_Residue, Residue, Mean_ddG_int_HM_A_C, mean_sel_coeff) %>%
#   group_by(TMP, Arabinose)  %>%
#   summarise(
#     corr_dim_int_s = cor.test(Mean_ddG_int_HM_A_C, mean_sel_coeff, method = 'spearman')$estimate, 
#     pval_corr_dim_int_s = cor.test(Mean_ddG_int_HM_A_C, mean_sel_coeff, method = 'spearman')$p.value
#   )
# write.table(test_corr_dim_int, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
#             file = 'Figures/2022-05-09_Supp_figures_paper/Tables_corr_ddG/corr_ddG_dim_int.csv')
# 
# test_corr_tet_int <- all_data_complete %>% ungroup() %>%
#   filter(!(is.na(Mean_ddG_int_HM_A_D))) %>% rowwise() %>%
#   select(TMP, Arabinose, Position, WT_Residue, Residue, Mean_ddG_int_HM_A_D, mean_sel_coeff) %>%
#   group_by(TMP, Arabinose)  %>%
#   summarise(
#     corr_tet_int_s = cor.test(Mean_ddG_int_HM_A_D, mean_sel_coeff, method = 'spearman')$estimate, 
#     pval_corr_tet_int_s = cor.test(Mean_ddG_int_HM_A_D, mean_sel_coeff, method = 'spearman')$p.value
#   )
# write.table(test_corr_tet_int, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
#             file = 'Figures/2022-05-09_Supp_figures_paper/Tables_corr_ddG/corr_ddG_tet_int.csv')
# 
# ## Use linear fits
# 
# # # Tests
# # a <- lm(formula = mean_sel_coeff ~ Mean_ddG_stab_HET, data = all_data_complete %>% 
# #           filter(!(is.na(Mean_ddG_stab_HET))))
# # b <- glm(formula = mean_sel_coeff ~ Mean_ddG_stab_HET,
# #          data = all_data_complete %>% 
# #           filter(!(is.na(Mean_ddG_stab_HET))))
# 
# ## Define a function to calculate the percentage of variance explained by ddGs
# ## in a linear fit
# var_explained <- function(formula_model, data_model){
#   glm.fit <- glm(formula = formula_model, data = data_model)
#   anova_table <- anova(glm.fit)
#   
#   # Extract the sums of squares
#   anova_ss <- anova_table$`Sum Sq`
#   
#   # Variance explained
# }
# 
# 
# test_lm_stab <- all_data_complete %>% ungroup() %>%
#   filter(!(is.na(Mean_ddG_stab_HET))) %>% rowwise() %>%
#   select(TMP, Arabinose, Position, WT_Residue, Residue, Mean_ddG_stab_HET, mean_sel_coeff) %>%
#   group_by(TMP, Arabinose)  %>%
#   summarise(
#     slope_stab_s = lm(formula = mean_sel_coeff ~ Mean_ddG_stab_HET)$coefficients[[2]],
#     glm_test = glm(formula = mean_sel_coeff ~ Mean_ddG_stab_HET, data = all_data_complete %>% 
#                      filter(!(is.na(Mean_ddG_stab_HET))))$
#   )
# write.table(test_lm_stab, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
#             file = 'Figures/2022-05-09_Supp_figures_paper/Tables_corr_ddG/corr_lm_stab.csv')
# 
# test_lm_dim_int <- all_data_complete %>% ungroup() %>%
#   filter(!(is.na(Mean_ddG_int_HM_A_C))) %>% rowwise() %>%
#   select(TMP, Arabinose, Position, WT_Residue, Residue, Mean_ddG_int_HM_A_C, mean_sel_coeff) %>%
#   group_by(TMP, Arabinose)  %>%
#   summarise(
#     slope_dim_int_s = lm(formula = mean_sel_coeff ~ Mean_ddG_int_HM_A_C)$coefficients[[2]] 
#   )
# write.table(test_lm_dim_int, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
#             file = 'Figures/2022-05-09_Supp_figures_paper/Tables_corr_ddG/corr_lm_dim_int.csv')
# 
# test_lm_tet_int <- all_data_complete %>% ungroup() %>%
#   filter(!(is.na(Mean_ddG_int_HM_A_D))) %>% rowwise() %>%
#   select(TMP, Arabinose, Position, WT_Residue, Residue, Mean_ddG_int_HM_A_D, mean_sel_coeff) %>%
#   group_by(TMP, Arabinose)  %>%
#   summarise(
#     slope_tet_int_s = lm(formula = mean_sel_coeff ~ Mean_ddG_int_HM_A_D)$coefficients[[2]] 
#   )
# write.table(test_lm_tet_int, append = F, quote = F, sep = '\t', row.names = F, col.names = T, 
#             file = 'Figures/2022-05-09_Supp_figures_paper/Tables_corr_ddG/corr_lm_tet_int.csv')

## Use a generalized linear model
data_glm_tmp <- all_data_complete %>% 
  filter(TMP == 10, !(is.na(Mean_ddG_stab_HET)), !(is.na(Mean_ddG_int_HM_A_C)), 
         !(is.na(Mean_ddG_int_HM_A_D)))
glm.fit <- glm(formula = mean_sel_coeff ~ Arabinose + Mean_ddG_stab_HET + 
                 Mean_ddG_int_HM_A_C + Mean_ddG_int_HM_A_D +
                 Arabinose*Mean_ddG_stab_HET +
                 Arabinose*Mean_ddG_int_HM_A_C +
                 Arabinose*Mean_ddG_int_HM_A_D,
               data = data_glm_tmp)
summary(glm.fit)
# anova(glm.fit)

data_glm_notmp <- all_data_complete %>% 
  filter(TMP == 0, !(is.na(Mean_ddG_stab_HET)), !(is.na(Mean_ddG_int_HM_A_C)), 
         !(is.na(Mean_ddG_int_HM_A_D)))
glm.fit_notmp <- glm(formula = mean_sel_coeff ~ Arabinose + Mean_ddG_stab_HET + 
                 Mean_ddG_int_HM_A_C + Mean_ddG_int_HM_A_D +
                 Arabinose*Mean_ddG_stab_HET +
                 Arabinose*Mean_ddG_int_HM_A_C +
                 Arabinose*Mean_ddG_int_HM_A_D,
               data = data_glm_notmp)
summary(glm.fit_notmp)

#### Figure S22B: Like figure 4B but without TMP ####

ddg_sel_coeff_data <- all_data_complete %>% 
  select(Position, WT_Residue, Residue, Timepoint, Arabinose, TMP, mean_sel_coeff, 
         Mean_ddG_stab_HET, Mean_ddG_int_HM_A_C, Mean_ddG_int_HM_A_D) %>%
  filter(!(is.na(Mean_ddG_stab_HET)), Timepoint == 10, TMP == 0)

summary(ddg_sel_coeff_data$Mean_ddG_int_HM_A_C)
summary(ddg_sel_coeff_data$Mean_ddG_int_HM_A_D)
summary(ddg_sel_coeff_data$Mean_ddG_stab_HET)

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

### Add a category of nonsense mutants
nonsense_mutants_notmp <- all_data_complete %>%
  filter(TMP == 0, Residue == '*', Position >= 30, Position <= 70) %>%
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

ddg_sel_coeff_data <- bind_rows(ddg_sel_coeff_data, nonsense_mutants_notmp) %>%
  mutate(
    bins_stab = factor(bins_stab, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')),
    bins_dim_int = factor(bins_dim_int, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    bins_tet_int = factor(bins_tet_int, levels = c('< 0', '(0,1]', '(1,2]', '(2,5]', '> 5', 'Stop')), 
    exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                             'Optimal', 'Overexpressed'))
  )

# Draw the figure for subunit stability
p_figs22B_stab <- ddg_sel_coeff_data %>% 
  ggplot(aes(x = bins_stab, y = mean_sel_coeff, fill = exp_level, colour = exp_level)) + 
  # geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.2)) +
  geom_boxplot(alpha = 0.4) +
  # xlab('\u0394\u0394G subunit stability [kcal /mol]') + ylab('s') +
  labs(y = expression(bolditalic('s')),
       x = expression(paste(bold('\u0394\u0394'), bolditalic('G'), 
                            bold(' subunit stability [kcal / mol]')))) +
  scale_fill_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'),
                    name = 'Expression level') +
  scale_colour_manual(values = c('#fed976', '#fd8d3c', '#bd0026', '#80001a', 'black'), 
                      name = 'Expression level') +
  scale_y_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0), limits = c(-1, 0.2)) +
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
p_figs22B_stab

# Draw the figure for the dimerization interface
p_figs22B_dim_int <- ddg_sel_coeff_data %>% 
  ggplot(aes(x = bins_dim_int, y = mean_sel_coeff, fill = exp_level, colour = exp_level)) + 
  # geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.2)) +
  geom_boxplot(alpha = 0.4) +
  # xlab('\u0394\u0394G dim. interface [kcal /mol]') + ylab('s') +
  labs(y = expression(bolditalic('s')),
       x = expression(paste(bold('\u0394\u0394'), bolditalic('G'), 
                            bold(' dim. interface [kcal / mol]')))) +
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
  ) +
  scale_y_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0), limits = c(-1, 0.2))
p_figs22B_dim_int


# Draw the figure for the dimerization interface
p_figs22B_tet_int <- ddg_sel_coeff_data %>% 
  ggplot(aes(x = bins_tet_int, y = mean_sel_coeff, fill = exp_level, colour = exp_level)) + 
  # geom_point(alpha = 0.6, position = position_jitterdodge(jitter.width = 0.2)) +
  geom_boxplot(alpha = 0.4) +
  # xlab('\u0394\u0394G tet. interface [kcal /mol]') + ylab('s') +
  labs(y = expression(bolditalic('s')),
       x = expression(paste(bold('\u0394\u0394'), bolditalic('G'), 
                            bold(' tet. interface [kcal / mol]')))) +
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
  ) +
  scale_y_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0), limits = c(-1, 0.2))
p_figs22B_tet_int

## Put the panels together
p_figs22a <- plot_grid(p_figs22_stab, p_figs22_dim_int, p_figs22_tet_int, ncol = 3)

legend_figs22b <- get_legend(p_figs22B_stab +
                               theme(legend.position = 'top', 
                                     legend.title = element_text(size = 26), 
                                     legend.text = element_text(size = 24)))
p_figs22b_panels <- plot_grid(p_figs22B_stab + theme(legend.position = 'none'),
                            p_figs22B_dim_int, p_figs22B_tet_int, ncol = 3)

p_figs22b <- plot_grid(legend_figs22b, p_figs22b_panels, nrow = 2, rel_heights = c(0.2, 1))

p_figs22 <- plot_grid(p_figs22a, p_figs22b, nrow = 2, labels = c('A', 'B'), 
                    label_size = 20, label_fontface = 'bold')
p_figs22
ggsave(p_figs22, device = cairo_pdf, width = 24, height = 14, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS22_noTMP_stab.pdf')
# ggsave(p_figs22, device = 'png', width = 24, height = 14, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS22_noTMP_stab.png')

#### Fig. S23: Similar to figure 3 but with the data without TMP ####

data_fig_4_notmp <- all_data_complete %>% ungroup() %>%
  filter(TMP == 0, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose)

entropy <- all_data_complete %>% ungroup() %>%
  filter(TMP == 0, Arabinose == 0.01) %>%
  group_by(Position) %>%
  summarise(Entropy = mean(Entropy))

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
annotation_name_gp = gpar(fontface = 'bold', fontsize = 14),
annotation_name_side = 'left',
annotation_name_rot = 90,
annotation_name_offset = c(Entropy = '0.15cm')
)

### 0.01 arabinose
# Separate the data for ara 0.2
data_part_1_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.01
data_part_2_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.01)


## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2_notmp) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final_notmp <- left_join(x = data_part_1_notmp, y = data_part_2_notmp, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)


# Need to convert to wide formatted data
data_fig_4_final_df_notmp <- data_fig_4_final_notmp %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_4_final_notmp <- as.matrix(data_fig_4_final_df_notmp %>% select(-Position))

rownames(data_fig_4_final_notmp) <- data_fig_4_final_df_notmp$Position

# Get a matrix of true/false values for the synonymous codons (I will reuse the)
fig_4_bool_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.01) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_4_bool_final_notmp <- as.matrix(fig_4_bool_notmp %>% select(-Position))

rownames(fig_4_bool_final_notmp) <- fig_4_bool_notmp$Position

# Need to reorder the columns in the matrices
data_fig_4_final_notmp <- data_fig_4_final_notmp[1:nrow(data_fig_4_final_notmp), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_4_bool_final_notmp <- fig_4_bool_final_notmp[1:nrow(fig_4_bool_final_notmp), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

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

### Repeat figure 4a with the annotation about expression level
p_figs23_ara0.01 <- Heatmap(
  t(data_fig_4_final_notmp), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-8, 8, length.out = 7) / 100, 
                   # breaks = seq(-4, 4, length.out = 7) / 100, 
                   colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  show_heatmap_legend = F, 
  row_title = "Residue",
  row_title_gp = gpar(fontsize=24, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  top_annotation = ha1,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_4_bool_final_notmp[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 4, 8) / 100,
    # at = c(-4, -2, 0, 2, 4) / 100,
    title = "\u0394s", 
    title_gp = gpar(fontsize = 16),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 14),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)
p_figs23_ara0.01

# A version of the same heatmap but with the scale I used for the data with TMP
# p_figs22_ara0.01_bigscale <- Heatmap(
#   t(data_fig_4_final_notmp), cluster_columns = F, cluster_rows = F, 
#   col = colorRamp2(# breaks = seq(-8, 8, length.out = 7) / 100, 
#     breaks = seq(-4, 4, length.out = 7) / 10, 
#     colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
#   show_column_names = T, row_names_side = 'left',
#   width=unit(31, 'cm'), height = unit(11.5, 'cm'),
#   border = T,
#   show_heatmap_legend = F, 
#   row_title = "Residue",
#   row_title_gp = gpar(fontsize=18, fontface = 'bold'),
#   row_names_rot = 90, 
#   row_names_centered = T,
#   row_names_gp = gpar(fontsize=13,fontface='bold'),
#   column_names_gp = gpar(fontsize=13,fontface='bold'),
#   top_annotation = ha1,
#   cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#     if (fig_4_bool_final_notmp[j,i]){
#       # grid.text('s', x, y)
#       grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
#     }      
#   },
#   heatmap_legend_param = list(
#     # at = c(-8, -4, 0, 4, 8) / 100,
#     at = c(-4, -2, 0, 2, 4) / 10,
#     title = "\u0394s", 
#     title_gp = gpar(fontsize = 16),
#     legend_height = unit(3.5, "cm"),
#     legend_width = unit(2, "cm"),
#     border='black',
#     lwd=1.7,
#     labels_gp = gpar(fontsize = 14),
#     title_position = "leftcenter-rot"
#   )
# )
# p_figs22_ara0.01_bigscale


### 0.025 arabinose
# Separate the data for ara 0.2
data_part_1_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.025
data_part_2_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.025)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2_notmp) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final_notmp <- left_join(x = data_part_1_notmp, y = data_part_2_notmp, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Need to convert to wide formatted data
data_fig_4_final_df_notmp <- data_fig_4_final_notmp %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_4_final_notmp <- as.matrix(data_fig_4_final_df_notmp %>% select(-Position))

rownames(data_fig_4_final_notmp) <- data_fig_4_final_df_notmp$Position

# Get a matrix of true/false values for the synonymous codons (I will reuse the)
fig_4_bool_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.025) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_4_bool_final_notmp <- as.matrix(fig_4_bool_notmp %>% select(-Position))

rownames(fig_4_bool_final_notmp) <- fig_4_bool_notmp$Position

# Need to reorder the columns in the matrices
data_fig_4_final_notmp <- data_fig_4_final_notmp[1:nrow(data_fig_4_final_notmp), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_4_bool_final_notmp <- fig_4_bool_final_notmp[1:nrow(fig_4_bool_final_notmp), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p_figs23_ara0.025 <- Heatmap(
  t(data_fig_4_final_notmp), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-8, 8, length.out = 7) / 100, 
    # breaks = seq(-4, 4, length.out = 7) / 100, 
                   colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  show_heatmap_legend = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=24, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_4_bool_final_notmp[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 4, 8) / 100,
    # at = c(-4, -2, 0, 2, 4) / 100,
    # title = "\u0394s", 
    title = expression(paste(bold('\u0394'), bolditalic('s'), sep = '')),
    title_gp = gpar(fontsize = 24, fontface = 'bold'),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 20),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)
p_figs23_ara0.025

# A version with the other scale
# p_figs22_ara0.025_bigscale <- Heatmap(
#   t(data_fig_4_final_notmp), cluster_columns = F, cluster_rows = F, 
#   col = colorRamp2(# breaks = seq(-8, 8, length.out = 7) / 100, 
#     breaks = seq(-4, 4, length.out = 7) / 10, 
#     colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
#   show_column_names = T, row_names_side = 'left',
#   width=unit(31, 'cm'), height = unit(11.5, 'cm'),
#   border = T,
#   show_heatmap_legend = T,
#   row_title = "Residue",
#   row_title_gp = gpar(fontsize=18, fontface = 'bold'),
#   row_names_rot = 90, 
#   row_names_centered = T,
#   row_names_gp = gpar(fontsize=13,fontface='bold'),
#   column_names_gp = gpar(fontsize=13,fontface='bold'),
#   cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#     if (fig_4_bool_final_notmp[j,i]){
#       # grid.text('s', x, y)
#       grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
#     }      
#   },
#   heatmap_legend_param = list(
#     # at = c(-8, -4, 0, 4, 8) / 100,
#     at = c(-4, -2, 0, 2, 4) / 10,
#     title = "\u0394s", 
#     title_gp = gpar(fontsize = 16),
#     legend_height = unit(3.5, "cm"),
#     legend_width = unit(2, "cm"),
#     border='black',
#     lwd=1.7,
#     labels_gp = gpar(fontsize = 14),
#     title_position = "leftcenter-rot"
#   )
# )
# p_figs22_ara0.025_bigscale

### 0.05 arabinose
# Separate the data for ara 0.2
data_part_1_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.05
data_part_2_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.05)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2_notmp) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final_notmp <- left_join(x = data_part_1_notmp, y = data_part_2_notmp, 
                              by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Need to convert to wide formatted data
data_fig_4_final_df_notmp <- data_fig_4_final_notmp %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# Need to convert the dataframe to a matrix
data_fig_4_final_notmp <- as.matrix(data_fig_4_final_df_notmp %>% select(-Position))

rownames(data_fig_4_final_notmp) <- data_fig_4_final_df_notmp$Position

# Get a matrix of true/false values for the synonymous codons (I will reuse the)
fig_4_bool_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.05) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_4_bool_final_notmp <- as.matrix(fig_4_bool_notmp %>% select(-Position))

rownames(fig_4_bool_final_notmp) <- fig_4_bool_notmp$Position

# Need to reorder the columns in the matrices
data_fig_4_final_notmp <- data_fig_4_final_notmp[1:nrow(data_fig_4_final_notmp), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_4_bool_final_notmp <- fig_4_bool_final_notmp[1:nrow(fig_4_bool_final_notmp), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

p_figs23_ara0.05 <- Heatmap(
  t(data_fig_4_final_notmp), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-8, 8, length.out = 7) / 100,
    # breaks = seq(-4, 4, length.out = 7) / 100,
                   colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  show_heatmap_legend = F,
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=24, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  # bottom_annotation = ha2,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_4_bool_final_notmp[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
    }      
  },
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 4, 8) / 100,
    # at = c(-4, -2, 0, 2, 4) / 100,
    title = "\u0394s", 
    title_gp = gpar(fontsize = 16),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 14),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)
p_figs23_ara0.05

# Same heatmap, big scale
# p_figs22_ara0.05_bigscale <- Heatmap(
#   t(data_fig_4_final_notmp), cluster_columns = F, cluster_rows = F, 
#   col = colorRamp2(# breaks = seq(-8, 8, length.out = 7) / 100,
#     breaks = seq(-4, 4, length.out = 7) / 10,
#     colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
#   show_column_names = T, row_names_side = 'left',
#   show_heatmap_legend = F,
#   width=unit(31, 'cm'), height = unit(11.5, 'cm'),
#   border = T,
#   row_title = "Residue",
#   row_title_gp = gpar(fontsize=18, fontface = 'bold'),
#   row_names_rot = 90, 
#   row_names_centered = T,
#   row_names_gp = gpar(fontsize=13,fontface='bold'),
#   column_names_gp = gpar(fontsize=13,fontface='bold'),
#   # bottom_annotation = ha2,
#   cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#     if (fig_4_bool_final_notmp[j,i]){
#       # grid.text('s', x, y)
#       grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
#     }      
#   },
#   heatmap_legend_param = list(
#     # at = c(-8, -4, 0, 4, 8) / 100,
#     at = c(-4, -2, 0, 2, 4) / 10,
#     title = "\u0394s", 
#     title_gp = gpar(fontsize = 16),
#     legend_height = unit(3.5, "cm"),
#     legend_width = unit(2, "cm"),
#     border='black',
#     lwd=1.7,
#     labels_gp = gpar(fontsize = 14),
#     title_position = "leftcenter-rot"
#   )
# )
# p_figs22_ara0.05_bigscale

### 0.4% arabinose
# Separate the data for ara 0.2
data_part_1_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.2)

# Separate the data for ara 0.4
data_part_2_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.4)

## Subtract the scores (just check that the order is the same)
# Change column names
colnames(data_part_2_notmp) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")

# Join
data_fig_4_final_notmp <- left_join(x = data_part_1_notmp, y = data_part_2_notmp, 
                                    by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
) %>% 
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)

# Need to convert to wide formatted data
data_fig_4_final_df_notmp <- data_fig_4_final_notmp %>%
  select(-WT_Residue, -Arabinose, -Arabinose_2, -mean_sel_coeff, -mean_sel_coeff_2) %>%
  pivot_wider(names_from = Residue, values_from = diffNormScore)

# View(data_fig_4_final_notmp %>% filter(Position == 2))

# Need to convert the dataframe to a matrix
data_fig_4_final_notmp <- as.matrix(data_fig_4_final_df_notmp %>% select(-Position))

rownames(data_fig_4_final_notmp) <- data_fig_4_final_df_notmp$Position

# Get a matrix of true/false values for the synonymous codons (I will reuse the)
fig_4_bool_notmp <- data_fig_4_notmp %>%
  filter(Arabinose == 0.4) %>%
  mutate(WT_check = (WT_Residue == Residue)) %>%
  select(-WT_Residue, -Arabinose, -mean_sel_coeff) %>%
  pivot_wider(names_from = Residue, values_from = WT_check)

fig_4_bool_final_notmp <- as.matrix(fig_4_bool_notmp %>% select(-Position))

rownames(fig_4_bool_final_notmp) <- fig_4_bool_notmp$Position

# Need to reorder the columns in the matrices
data_fig_4_final_notmp <- data_fig_4_final_notmp[1:nrow(data_fig_4_final_notmp), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]
fig_4_bool_final_notmp <- fig_4_bool_final_notmp[1:nrow(fig_4_bool_final_notmp), c(1,7,2,19,11,9,12,3,14,20,6,21,17,18,13,15,8,10,16,4,5)]

# Prepare the bottom annotation
data_interfaces_final_notmp <- read_delim('Data/Complete_datasets_TMP0_TMP10/data_annotation_2.txt', delim = '\t')

# col_fun = colorRamp2(c(0, 1), c("white", "#595959"))

# # Prepare the annotation on interfaces
# # Interface A-C will be interface 1
# # Interface A-D will be interface 2
# ha2 <- HeatmapAnnotation(`Interface 1` = data_interfaces_final$`A,C`,
#                          `Interface 2` = data_interfaces_final$`A,D`,
#                          DHF = data_interfaces_final$DHF,
#                          NADPH = data_interfaces_final$NADPH,
#                          `Catalytic residues` = data_interfaces_final$Cat_residues,
#                          `Disordered region` = data_interfaces_final$Disordered_region,
#                          show_annotation_name = T,
#                          annotation_name_gp = gpar(fontface = 'bold', fontsize = 14),
#                          annotation_name_side = 'left',
#                          show_legend = FALSE,
#                          col = list(`Interface 1` = colorRamp2(c(0, 1), c("white", "#009E73")), 
#                                     `Interface 2` = colorRamp2(c(0, 1), c("white", "#0072B2")),
#                                     DHF = colorRamp2(c(0, 1), c("white", "#E69F00")),
#                                     NADPH = colorRamp2(c(0, 1), c("white", "#D55E00")),
#                                     `Catalytic residues` = colorRamp2(c(0, 1), c("white", "#56B4E9")),
#                                     `Disordered region` = colorRamp2(c(0, 1), c("white", "#CC79A7"))
#                          ),
#                          gp = gpar(col = "black")
# )

# Color order
# c('Disordered region', 'Catalytic residues', 'DHF', 'Interface 1', 'Interface 2', 'NADPH', 'Unannotated')
# c('#CC79A7', '#56B4E9', '#E69F00', '#009E73', '#0072B2', '#D55E00', '#000000')

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



p_figs23_ara0.4 <- Heatmap(
  t(data_fig_4_final_notmp), cluster_columns = F, cluster_rows = F, 
  col = colorRamp2(breaks = seq(-8, 8, length.out = 7) / 100,
    # breaks = seq(-4, 4, length.out = 7) / 100,
    colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
  show_column_names = T, row_names_side = 'left',
  show_heatmap_legend = F,
  width=unit(31, 'cm'), height = unit(11.5, 'cm'),
  border = T,
  row_title = "Residue",
  row_title_gp = gpar(fontsize=24, fontface = 'bold'),
  row_names_rot = 90, 
  row_names_centered = T,
  row_names_gp = gpar(fontsize=18,fontface='bold'),
  column_names_gp = gpar(fontsize=18,fontface='bold'),
  bottom_annotation = ha2,
  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    if (fig_4_bool_final_notmp[j,i]){
      # grid.text('s', x, y)
      grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
      
    }      
  },
  heatmap_legend_param = list(
    at = c(-8, -4, 0, 4, 8) / 100,
    # at = c(-4, -2, 0, 2, 4) / 100,
    title = "\u0394s", 
    title_gp = gpar(fontsize = 16),
    legend_height = unit(3.5, "cm"),
    legend_width = unit(2, "cm"),
    border='black',
    lwd=1.7,
    labels_gp = gpar(fontsize = 14),
    title_position = "leftcenter-rot"
  ), 
  column_labels = column_labels
)
p_figs23_ara0.4

# Same figure with bigger scale
# p_figs22_ara0.4_bigscale <- Heatmap(
#   t(data_fig_4_final_notmp), cluster_columns = F, cluster_rows = F, 
#   col = colorRamp2(# breaks = seq(-8, 8, length.out = 7) / 100,
#     breaks = seq(-4, 4, length.out = 7) / 10,
#     colors = rev(brewer.pal(n = 7, name = 'BrBG'))), 
#   show_column_names = T, row_names_side = 'left',
#   show_heatmap_legend = F,
#   width=unit(31, 'cm'), height = unit(11.5, 'cm'),
#   border = T,
#   row_title = "Residue",
#   row_title_gp = gpar(fontsize=18, fontface = 'bold'),
#   row_names_rot = 90, 
#   row_names_centered = T,
#   row_names_gp = gpar(fontsize=13,fontface='bold'),
#   column_names_gp = gpar(fontsize=13,fontface='bold'),
#   bottom_annotation = ha2,
#   cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
#     if (fig_4_bool_final_notmp[j,i]){
#       # grid.text('s', x, y)
#       grid.points(x, y, pch = 19, size = unit(0.75, 'char'))
#       
#     }      
#   },
#   heatmap_legend_param = list(
#     # at = c(-8, -4, 0, 4, 8) / 100,
#     at = c(-4, -2, 0, 2, 4) / 10,
#     title = "\u0394s", 
#     title_gp = gpar(fontsize = 16),
#     legend_height = unit(3.5, "cm"),
#     legend_width = unit(2, "cm"),
#     border='black',
#     lwd=1.7,
#     labels_gp = gpar(fontsize = 14),
#     title_position = "leftcenter-rot"
#   )
# )
# p_figs22_ara0.4_bigscale

## Put the three figures together
ht_list = p_figs23_ara0.01 %v% p_figs23_ara0.025 %v% p_figs23_ara0.05 %v% p_figs23_ara0.4
p_figs23_heatmaps <- grid.grabExpr(
  draw(ht_list,
       row_title_gp = gpar(fontsize=20, fontface = 'bold'),
       ht_gap = unit(1, "cm"))
)

# Panels with bigger scale
# ht_list_bigscale = p_figs22_ara0.01_bigscale %v% p_figs22_ara0.025_bigscale %v% 
#   p_figs22_ara0.05_bigscale %v% p_figs22_ara0.4_bigscale
# p_figs22_heatmaps_bigscale <- grid.grabExpr(
#   draw(ht_list_bigscale,
#        row_title_gp = gpar(fontsize=20, fontface = 'bold'),
#        ht_gap = unit(1, "cm"))
# )

### Add text labels
# ## Arabinose 0.01
# text_fig_ara0.01 <- ggplot() + draw_label(expression(paste('s'[weak], ' - s'[opt], sep = '')),
#                                           x = 0.7, y = 0.45, 
#                                           fontface = 'bold', size = 35, angle = 90, colour = '#fed976') +
#   theme(axis.line = element_blank())
# text_fig_ara0.01
# 
# # Arabinose 0.025
# text_fig_ara0.025 <- ggplot() + draw_label(expression(paste('s'[subopt], ' - s'[opt], sep = '')),
#                                            x = 0.7, y = 0.55, 
#                                            fontface = 'bold', size = 35, angle = 90, colour = '#fd8d3c') +
#   theme(axis.line = element_blank())
# text_fig_ara0.025
# 
# # Arabinose 0.05
# text_fig_ara0.05 <- ggplot() + draw_label(expression(paste('s'[near-opt], ' - s'[opt], sep = '')),
#                                           x = 0.7, y = 0.7, 
#                                           fontface = 'bold', size = 35, angle = 90, colour = '#bd0026') +
#   theme(axis.line = element_blank())
# text_fig_ara0.05
# 
# # Arabinose 0.4
# text_fig_ara0.4 <- ggplot() + draw_label(expression(paste('s'[over], ' - s'[opt], sep = '')),
#                                           x = 0.7, y = 0.8, 
#                                           fontface = 'bold', size = 35, angle = 90, colour = 'black') +
#   theme(axis.line = element_blank())
# text_fig_ara0.4

## Arabinose 0.01
text_fig_ara0.01 <- ggplot() + 
  draw_label(# expression(paste('s'[weak], ' - s'[opt], sep = '')),
    expression(atop(paste(bold('\u0394'), bolditalic(s[weak]), sep = ''),
                    paste(bold('('), 
                          bolditalic(s[weak] - s[opt]), bold(')'), sep = ''))),
    # x = 0.7, y = 0.45, 
    x = 0.7, y = 0.3,
    fontface = 'bold', size = 35, angle = 90, colour = '#fed976') +
  theme(axis.line = element_blank())
text_fig_ara0.01

# Arabinose 0.025
text_fig_ara0.025 <- ggplot() + draw_label(# expression(paste('s'[subopt], ' - s'[opt], sep = '')),
  expression(atop(paste(bold('\u0394'), bolditalic(s[subopt]), sep = ''),
                  paste(bold('('), 
                        bolditalic(s[subopt] - s[opt]), bold(')'), sep = ''))),
  # x = 0.7, y = 0.55, 
  x = 0.7, y = 0.4,
  fontface = 'bold', size = 35, angle = 90, colour = '#fd8d3c') +
  theme(axis.line = element_blank())
text_fig_ara0.025

# Arabinose 0.05
text_fig_ara0.05 <- ggplot() + draw_label(# expression(paste('s'[near-opt], ' - s'[opt], sep = '')),
  expression(atop(paste(bold('\u0394'), bolditalic(s[near-opt]), sep = ''),
                  paste(bold('('), 
                        bolditalic(s[near-opt] - s[opt]), bold(')'), sep = ''))),
  # x = 0.7, y = 0.7, 
  x = 0.7, y = 0.4,
  fontface = 'bold', size = 35, angle = 90, colour = '#bd0026') +
  theme(axis.line = element_blank())
text_fig_ara0.05

# Arabinose 0.4
text_fig_ara0.4 <- ggplot() + draw_label(# expression(paste('s'[over], ' - s'[opt], sep = '')),
  expression(atop(paste(bold('\u0394'), bolditalic(s[over]), sep = ''),
                  paste(bold('('), 
                        bolditalic(s[over] - s[opt]), bold(')'), sep = ''))),
  # x = 0.7, y = 0.8, 
  x = 0.7, y = 0.35,
  fontface = 'bold', size = 35, angle = 90, colour = 'black') +
  theme(axis.line = element_blank())
text_fig_ara0.4




## Add a panel with the distributions of deltaS, separated by region as in fig3B

#### Start figure S23B ####

#### Check the distribution of delta(DMS scores) for critical sites ####

data_fig_3 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 0, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, Arabinose, Secondary_structure, rSASA)

# Separate the data for ara 0.2
data_part_1 <- data_fig_3 %>%
  filter(Arabinose == 0.2)

## Separate the data for ara 0.01
# data_part_2 <- data_fig_3 %>%
#   filter(Arabinose == 0.01)

## Separate the data for ara 0.4
data_part_2 <- data_fig_3 %>%
 filter(Arabinose == 0.4)

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
  pivot_longer(cols = c('A,C', 'A,D', 'DHF', 'NADPH', 'Cat_residues', 'Buried',
                        'Only NADPH', 'Disordered_region','Unannotated'), 
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

# Statistical comparisons
comps <- compare_means(diffNormScore~Site, data = critical_sites_deltaDMS,
                       paired = F) %>%
  mutate(p.format = ifelse(p < 2.2e-16, 'p < 2.2e-16',
                           ifelse(p > 0.01, str_c('p = ', round(as.numeric(p), 2), sep = ''),
                                  sprintf("p = %2.1e", as.numeric(p)))
  )
  # y_pos = rep(c(0.6, 0.7, 0.8), 8)
  )

p_figs23b <- critical_sites_deltaDMS %>% rowwise() %>%
  mutate(mut_id = str_c(WT_Residue, Position, Residue, sep = '')) %>%
  filter(Site != 'Only NADPH') %>%
  mutate(Site = factor(Site, levels = c('Disordered region', 'Catalytic residues',
                                        'DHF binding', 'NADPH binding',
                                        'Buried residues',
                                         'Tetramerization interface', 'Dimerization interface','Unannotated')
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
        plot.margin = margin(t = 1, r = 5, b = 0, l = 9, 'cm')) +
  xlab('') + # ylab('\u0394s') +
  # labs(y = expression(s[weak] - s[opt])) +
  # labs(y = expression(paste(bold('\u0394'), bolditalic(s[weak]),
  #                           bold(' ('), bolditalic(s[weak] - s[opt]), bold(')'), 
  #                           sep = ''))) +
  labs(y = expression(paste(bold('\u0394'), bolditalic(s[over]),
                            bold(' ('), bolditalic(s[over] - s[opt]), bold(')'), 
                            sep = ''))) +
  ylim(-0.08, 0.08) +  
  # ylim(-0.4, 0.4) +  
  annotate('text', 
           x = 3, y = 0.04,
           # x = 3, y = 0.2,
           # label = 's[weak] > s[opt]',
           # label = expression(italic(s[weak] > s[opt])),
           label = expression(italic(s[over] > s[opt])),
           parse = T, size = 8) +
  annotate('text', 
           x = 3, y = -0.04,
           # x = 3, y = -0.2,
           # label = 's[weak] < s[opt]',
           # label = expression(italic(s[weak] < s[opt])),
           label = expression(italic(s[over] < s[opt])),
           parse = T, size = 8)
p_figs23b

#### End figure s23B ####

### Put everything together
# p_figs23_text <- plot_grid(text_fig_ara0.01, text_fig_ara0.025, text_fig_ara0.05, text_fig_ara0.4,
#                            nrow = 4, rel_heights = c(1.1, 1, 1.1, 0.9))
# 
# p_figs23a <- plot_grid(p_figs23_text + theme(plot.margin = margin(t = 0, r = -3, b = 0, l = 8, 'cm')),
#                        p_figs23_heatmaps, ncol = 2, rel_widths = c(0.2, 1), 
#                        labels = c('', 'A'), label_size = 20, label_fontface = 'bold')
# 
# p_figs23b_null <- plot_grid(NULL, p_figs23b, rel_widths = c(0.2, 1), ncol = 2, 
#                             labels = c('', 'B'), label_size = 20, label_fontface = 'bold')
# 
# p_figs23 <- plot_grid(p_figs23a, p_figs23b_null, nrow = 2, rel_heights = c(4, 1.3))
# 
# p_figs23
# 
# ggsave(p_figs23, width = 23, height = 30, dpi = 300, device = cairo_pdf, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_supp_allDiff_noSecStruc_buried.pdf')

### Save the text labels and the panels separately
p_figs23_text <- plot_grid(text_fig_ara0.01, text_fig_ara0.025, text_fig_ara0.05, text_fig_ara0.4,
                           nrow = 4, rel_heights = c(1.1, 1, 1.1, 1))

# p_figs23a <- plot_grid(p_figs23_text + theme(plot.margin = margin(t = 0, r = -3, b = 0, l = 8, 'cm')),
#                        p_figs23_heatmaps, ncol = 2, rel_widths = c(0.2, 1), 
#                        labels = c('', 'A'), label_size = 20, label_fontface = 'bold')

p_figs23 <- plot_grid(p_figs23_heatmaps, 
                      p_figs23b + theme(plot.margin = margin(t = 1, r = 10, b = 0, l = 14, 'cm')),
                      nrow = 2, rel_heights = c(4, 1.3), 
                      labels = c('A', 'B'), label_size = 20, label_fontface = 'bold')

p_figs23
# ggsave(p_figs23, width = 23, height = 30, dpi = 300, device = cairo_pdf, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_supp_allDiff_noSecStruc_buried_panels.pdf')
# ggsave(p_figs23_text, width = 23, height = 31, dpi = 300, device = cairo_pdf, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_supp_allDiff_noSecStruc_buried_text.pdf')

ggsave(p_figs23, width = 23, height = 30, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_supp_allDiff_noSecStruc_buried_panels_over.pdf')
ggsave(p_figs23_text, width = 23, height = 31, dpi = 300, device = cairo_pdf, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_supp_allDiff_noSecStruc_buried_text.pdf')


## Draw the figure with the bigger scale
# p_figs22a_bigscale <- plot_grid(p_figs22_text + theme(plot.margin = margin(t = 0, r = -3, b = 0, l = 8, 'cm')),
#                        p_figs22_heatmaps_bigscale, ncol = 2, rel_widths = c(0.2, 1), 
#                        labels = c('', 'A'), label_size = 20, label_fontface = 'bold')
# 
# p_figs22b_null_bigscale <- plot_grid(NULL, p_figs22b, rel_widths = c(0.2, 1), ncol = 2, 
#                             labels = c('', 'B'), label_size = 20, label_fontface = 'bold')
# 
# p_figs22_bigscale <- plot_grid(p_figs22a_bigscale, p_figs22b_null_bigscale, nrow = 2, rel_heights = c(1, 0.3))
# 
# ggsave(p_figs22_bigscale, width = 23, height = 30, dpi = 300, device = cairo_pdf, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS22_supp_allDiff_noSecStruc_bigscale.pdf')

#### Figure S24: Growth curves for WT and E2R mutants ####

# Load the growth curves with TMP
plate.ind <- 'Data/Growth_curves_May2022/Index_growth_curves_11_05_2022.xlsx'
file.od <- 'Data/Growth_curves_May2022/Growth_curves_11_05_2022.xlsx'

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
plate.ind <- 'Data/Growth_curves_May2022/Index_growth_curves_20_05_2022.xlsx'
file.od <- 'Data/Growth_curves_May2022/Growth_curves_20_05_2022.xlsx'

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
data_figs24 <- bind_rows(data.od1, data.od2)

# Summarise to have only one data point for AUC for each well
data_figs24_summary <- data_figs24 %>% ungroup() %>% 
  group_by(Well, Arabinose, TMP, Mutant) %>%
  summarise(auc = mean(auc_e))

# p_figs23 <- data_figs23_summary %>% rowwise() %>%
#   filter(!(is.na(Mutant))) %>%
#   mutate(TMP = str_c('TMP = ', TMP, ' µg / mL', sep = '')) %>%
#   ggplot(aes(x = as.factor(Arabinose), y = auc, fill = Mutant, colour = Mutant)) +
#   # geom_boxplot(outlier.shape = NA, alpha = 0.4, width = 0.8) +
#   geom_point(position = position_jitterdodge()) +
#   facet_wrap(~TMP, scales = 'free') +
#   theme(axis.title = element_text(size = 20, face = 'bold'), 
#         axis.text = element_text(size = 18), 
#         panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         legend.position = 'top', 
#         legend.justification = 'center', 
#         legend.title = element_text(size = 20), 
#         legend.text = element_text(size = 18), 
#         strip.text = element_text(size = 20, face = 'bold'), 
#         strip.background = element_rect(fill = 'white')) +
#   xlab('Arabinose') + ylab('Growth in liquid culture (AUC)') +
#   ylim(0, 47)
# p_figs23  
# ggsave(plot = p_figs23, device = cairo_pdf, width = 17, height = 7, dpi = 300,
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_AUC_WT_E2R.pdf')
# ggsave(plot = p_figs23, device = 'png', width = 17, height = 7, dpi = 300,
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_AUC_WT_E2R.png')

# Summarize to show the average of the three replicates
data_figs24_sum_curves <- data_figs24 %>% ungroup() %>%
  group_by(Arabinose, TMP, Mutant, time) %>%
  summarise(OD = mean(OD))

# Show the growth curves
# p_figs23_curves <- data_figs23_sum_curves %>% ungroup() %>% rowwise() %>%
#   filter(!(is.na(Mutant)), !(is.na(Arabinose)), time <= 13.5) %>%
#   # filter(time <= 13.5) %>%
#   mutate(
#     TMP = as.factor(TMP),
#     Arabinose = str_c(Arabinose, '% arabinose')
#   ) %>%
#   mutate(Arabinose = factor(Arabinose, levels = c('0% arabinose',
#                                                   '0.001% arabinose', 
#                                                   '0.0025% arabinose', 
#                                                   '0.005% arabinose',
#                                                   '0.01% arabinose',
#                                                   '0.025% arabinose',
#                                                   '0.05% arabinose',
#                                                   '0.2% arabinose',
#                                                   '0.4% arabinose'))) %>%
#   ggplot(aes(x = time, y = OD, colour = Arabinose,
#              group = interaction(Mutant, Arabinose, TMP))) +
#   facet_grid(TMP~Mutant, scales = 'free') +
#   geom_line() +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(face = 'bold', size = 20),
#         axis.text = element_text(size = 18),
#         axis.line = element_line(),
#         legend.position = 'top',
#         legend.title = element_text(size = 20),
#         strip.text = element_text(size = 20, face = 'bold'),
#         strip.background = element_rect(fill = 'white'),
#         legend.text = element_text(size = 18),
#         legend.justification = 0.5) +
#   guides(size = 'none', linetype = 'none', alpha = 'none') +
#   xlab('Time (h)') + ylab('OD / mL') + ylim(0, 6) +
#   geom_vline(xintercept = 13.5, linetype = 'dashed', colour = 'red')
# p_figs23_curves
# ggsave(plot = p_figs23_curves, device = cairo_pdf, width = 21, height = 14, dpi = 300,
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_growthCurves_DfrB1_WT_E2R_t13.5.pdf')
#        # filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_growthCurves_DfrB1_WT_E2R_t20.pdf')
# ggsave(plot = p_figs23_curves, device = 'png', width = 21, height = 14, dpi = 300,
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_growthCurves_DfrB1_WT_E2R_t13.5.png')
#        # filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_growthCurves_DfrB1_WT_E2R_t20.png')

## Try a log scale for the curves
# p_figs23_curves <- data_figs23_sum_curves %>% ungroup() %>% rowwise() %>%
#   filter(!(is.na(Mutant)), !(is.na(Arabinose)), time <= 13.5) %>%
#   # filter(time <= 13.5) %>%
#   mutate(
#     TMP = as.factor(TMP),
#     Arabinose = str_c(Arabinose, '% arabinose')
#   ) %>%
#   mutate(Arabinose = factor(Arabinose, levels = c('0% arabinose',
#                                                   '0.001% arabinose', 
#                                                   '0.0025% arabinose', 
#                                                   '0.005% arabinose',
#                                                   '0.01% arabinose',
#                                                   '0.025% arabinose',
#                                                   '0.05% arabinose',
#                                                   '0.2% arabinose',
#                                                   '0.4% arabinose'))) %>%
#   ggplot(aes(x = time, y = log2(OD), colour = Arabinose,
#              group = interaction(Mutant, Arabinose, TMP))) +
#   facet_grid(TMP~Mutant, scales = 'free') +
#   geom_line() +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(face = 'bold', size = 20),
#         axis.text = element_text(size = 18),
#         axis.line = element_line(),
#         legend.position = 'top',
#         legend.title = element_text(size = 20),
#         strip.text = element_text(size = 20, face = 'bold'),
#         strip.background = element_rect(fill = 'white'),
#         legend.text = element_text(size = 18),
#         legend.justification = 0.5) +
#   guides(size = 'none', linetype = 'none', alpha = 'none') +
#   xlab('Time (h)') + ylab('OD') + # ylim(0, 3) +
#   geom_vline(xintercept = 13.5, linetype = 'dashed', colour = 'red')
# p_figs23_curves
# ggsave(plot = p_figs23_curves, device = cairo_pdf, width = 21, height = 14, dpi = 300,
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_growthCurves_DfrB1_WT_E2R_t13.5_logscale.pdf')
# ggsave(plot = p_figs23_curves, device = 'png', width = 21, height = 14, dpi = 300,
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_growthCurves_DfrB1_WT_E2R_t13.5_logscale.png')

#### Show the percentage of growth recovery ####

# Pivot the area under the curve to calculate the differences
data_figs24_wide <- data_figs24_summary %>% ungroup() %>%
  group_by(Well, Arabinose, Mutant) %>%
  filter(!(is.na(Mutant))) %>%
  pivot_wider(names_from = TMP, values_from = auc, names_prefix = 'AUC_TMP_')

# Calculate the difference between the data with and without TMP
# and then the percentage of growth recovery
data_figs24_wide %<>% mutate(diff_auc = AUC_TMP_10 - AUC_TMP_0)#  %>%
# mutate(pct_recovery = 100 * ( 1 - diff_auc / AUC_TMP_0))

# Add a column for the maximum difference (0% arabinose)
max_diff <- data_figs24_wide %>% filter(Arabinose == 0) %>% ungroup() %>%
  group_by(Arabinose, Mutant) %>%
  summarise(max_diff = median(diff_auc))

# Use a join to add the maximum difference
data_figs24_final <- left_join(x = data_figs24_wide, 
                               y = max_diff %>% ungroup() %>% select(-Arabinose), 
                               by = c('Mutant' = 'Mutant'))

data_figs24_final %<>% mutate(pct_recovery = 100 * (1 - (diff_auc / max_diff)))

# Draw the figure
p_figs24a <- data_figs24_final %>% 
  ggplot(aes(x = as.factor(Arabinose), y = pct_recovery, colour = Mutant)) +
  geom_point(position = position_jitterdodge(), size = 3) +
  # facet_wrap(~TMP, scales = 'free') +
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
  xlab('Arabinose (% m/v)') + ylab('Recovered growth (%)')
p_figs24a
# ggsave(p_figs23a, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_WT_E2R_recovery.pdf')
# ggsave(p_figs23a, device = 'png', width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-05-09_Supp_figures_paper/FigS23_WT_E2R_recovery.png')

## Show figure S24A as means with standard errors of the mean
data_figs24_final_summ <- data_figs24_final %>% ungroup() %>%
  group_by(Arabinose, Mutant) %>%
  summarise(mean_recovery = mean(pct_recovery), 
            sem_recovery = sd(pct_recovery) / sqrt(n()), 
            num_samples = n())

p_figs24a <- data_figs24_final_summ %>% 
  ggplot(aes(x = as.factor(Arabinose), y = mean_recovery, colour = Mutant, 
             ymax = mean_recovery + sem_recovery, 
             ymin = mean_recovery - sem_recovery)) +
  geom_point(size = 1.5) +
  geom_errorbar(width = 0.2) +
  # geom_point(position = position_jitterdodge(), size = 3) +
  # facet_wrap(~TMP, scales = 'free') +
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
  xlab('Arabinose (% m/v)') + ylab('Recovered growth (%)')
p_figs24a

#### Figure S24B: Look at the cost of overexpressing E2R ####

# Need to calculate the median growth recovery at 0.001 arabinose
med_opt_E2R <- data_figs24_final %>% ungroup() %>%
  filter(Mutant == 'E2R', Arabinose == 0.001) %>%
  group_by(Mutant, Arabinose) %>%
  summarise(med_recovery = median(pct_recovery))

cost_e2r <-data_figs24_final %>% mutate(opt_recovery = med_opt_E2R$med_recovery[1])

cost_e2r %<>% mutate(cost = opt_recovery - pct_recovery)

p_figs24b <- cost_e2r %>% 
  filter(Mutant == 'E2R', Arabinose > 0) %>%
  ggplot(aes(x = as.factor(Arabinose), y = cost)) +
  # geom_point(position = position_jitterdodge(), size = 3) +
  geom_jitter(size = 3, width = 0.2) +
  stat_summary(fun = mean, colour = 'red') +
  stat_summary(fun = mean, geom = 'path',
               mapping = aes(group = -1), colour = 'red') +
  # facet_wrap(~TMP, scales = 'free') +
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
  xlab('Arabinose (% m/v)') + ylab('E2R expression cost (% recovered growth)')
p_figs24b

p_figs24 <- plot_grid(p_figs24a, p_figs24b, nrow = 2, labels = c('A', 'B'), 
                      label_size = 20, label_fontface = 'bold')

ggsave(p_figs24, device = cairo_pdf, width = 10, height = 15, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS24_WT_E2R_recovery.pdf')
ggsave(p_figs24, device = 'png', width = 10, height = 15, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS24_WT_E2R_recovery.png')

#### Repeat the ANOVAs with the data without TMP ####

data_fig_2 <- all_data_complete %>% ungroup() %>%
  filter(TMP == 0, Timepoint == 10) %>%
  select(Position, WT_Residue, Residue, mean_sel_coeff, sem_sel_coeff, Arabinose)

data_fig2c <- data_fig_2 %>% ungroup() %>%
  group_by(Position, Arabinose) %>%
  summarise(meanSelCoeff = mean(mean_sel_coeff))

data_fig2c_exp <- data_fig2c %>% 
  mutate(Expression_level = ifelse(Arabinose == 0.01, 'Weak', 
                                   ifelse(Arabinose == 0.025, 'Suboptimal', 
                                          ifelse(Arabinose == 0.05, 'Near-optimal', 
                                                 ifelse(Arabinose == 0.2, 'Optimal', 
                                                        ifelse(Arabinose == 0.4, 'Overexpressed', NA)))))) %>%
  mutate(Expression_level = factor(Expression_level,
                                   levels = c('Weak', 'Suboptimal',
                                              'Near-optimal', 'Optimal', 'Overexpressed')))

## Draw the figure
p_figs10a <- data_fig2c_exp %>% 
  ggplot(aes(x = Expression_level, y = meanSelCoeff, fill = Expression_level)) +
  geom_violin(alpha = 0.8) +
  stat_summary(fun="median", geom="point", size = 5) +
  geom_point() +
  geom_line(aes(group = Position)) +
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
  ylab('Mean s per position')
p_figs10a

#### Figure S10B: Volcano plots ####
data_fig2b <- data_fig_2_final %>% ungroup() %>%
  mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff,
         Arabinose_2 = str_c(Arabinose_2, '% arabinose', sep = ''))


# Relevel the data set
data_fig2b %<>% mutate(Arabinose_2 = factor(Arabinose_2, 
                                            # levels = c('0.025% arabinose', '0.05% arabinose', '0.2% arabinose'))
                                            levels = c('0.01% arabinose', '0.025% arabinose', '0.05% arabinose', '0.4% arabinose'))
)

data_fig2b_exp <- data_fig2b %>% separate(col = Arabinose_2, into = c('Arabinose_num', 'tmp'), 
                                          sep = '% arabinose') %>%
  mutate(Arabinose_num = as.numeric(Arabinose_num), 
         exp_level = ifelse(Arabinose_num == 0.01, 'Weak', 
                            ifelse(Arabinose_num == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose_num == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose_num == 0.2, 'Optimal',
                                                 ifelse(Arabinose_num == 0.4, 'Overexpressed', NA))))))

# Load data for all replicates from both sequencers
all_data_all_reps <- read_delim('../R67_DMS_February2022/Complete_datasets_MiSeq_NovaSeq/all_data_all_reps_bothSequencers.txt', 
                                delim = '\t')

all_data_all_reps_TMP0 <- all_data_all_reps %>% rowwise() %>%
  filter(TMP == 0, Timepoint == 10, WT_Residue != Residue, Residue != '*') %>%
  mutate(Genotype = str_c(WT_Residue, Residue, Position)) %>%
  select(ID, Genotype, Timepoint, TMP, Arabinose, sel_coeff, Sequencer)

genotypes <- unique(all_data_all_reps_TMP0$Genotype)
all_anovas <- c()

## Do an ANOVA of each mutant
for(genotype in genotypes){
  ## Compile a table with each mutant and its p-value from the ANOVA
  # Select data from one genotype to do the ANOVA
  data_genotype <- all_data_all_reps_TMP0 %>% filter(Genotype == genotype) %>% select(-Timepoint, -TMP) %>%
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
  
  # Save values from the ANOVA
  df_arabinose <- anova_test$Df[1]
  df_residuals <- anova_test$Df[2]
  
  sum_sq_arabinose <- anova_test$`Sum Sq`[1]
  sum_sq_residuals <- anova_test$`Sum Sq`[2]
  
  mean_sq_arabinose <- anova_test$`Mean Sq`[1]
  mean_sq_residuals <- anova_test$`Mean Sq`[2]
  
  f_value <- anova_test$`F value`[1]
  p_value <- anova_test$`Pr(>F)`[1]
  
  new_row <- data.frame(Genotype = c(genotype), df_arabinose = c(df_arabinose), df_residuals = c(df_residuals),
                        sum_sq_arabinose = c(sum_sq_arabinose), sum_sq_residuals = c(sum_sq_residuals),
                        f_value = c(f_value), p_val = c(p_value))
  all_anovas <- bind_rows(all_anovas, new_row)
}

## Apply the Benjamini-Hochberg correction for multiple hypotheses
fdr_test <-  p.fdr(pvalues = all_anovas$p_val, adjust.method = 'BH', threshold = 0.05)
all_anovas$p.adj <- fdr_test$fdrs

summary(all_anovas$p.adj)
sum(all_anovas$p.adj < 0.05)

# Prepare the data
temp_data <- data_fig_2 %>% rowwise() %>%
  mutate(ID = str_c(Position, Residue, sep = ''))

# Let's add a column with the difference in scores at low and high expression
score_diff <- data_fig2b %>% 
  filter(Arabinose_2 == '0.01% arabinose') %>%
  mutate(mut_check_diff = ifelse(abs(diffNormScore) > 0.3, TRUE, FALSE))

score_significant <- left_join(x = data_fig2b_exp %>% filter(Arabinose_num == 0.01), 
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
data_fig2d <- left_join(x = temp_data, 
                        y= score_significant %>% select(Position, Residue, p.adj, mut_check),
                        by = c('Position' = 'Position', 'Residue' = 'Residue'))


# Add the expression level
data_fig2d_exp <- data_fig2d %>%
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal',
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA))))))

# Add the data about the minimum change in expression
data_fig2d_exp <- left_join(x = data_fig2d_exp, 
                            y = score_diff %>% select(Position, Residue, diffNormScore, mut_check_diff), 
                            by = c('Position' = 'Position', 'Residue' = 'Residue')) %>%
  rowwise() %>%
  mutate(mut_check_final = and(mut_check, mut_check_diff))

## Try a volcano plot 
p_figs10b <- data_fig2d_exp %>%
  filter(Arabinose == 0.01) %>%
  ggplot(aes(x = diffNormScore, y = -log10(p.adj), colour = mut_check_final)) +
  geom_point() +
  geom_vline(xintercept = -0.3, linetype = 'dashed') +
  geom_vline(xintercept = 0.3, linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
  # xlab('\u0394s') + 
  ylab('-log10(p.adj)') +
  labs(x = expression(s[weak] - s[opt])) +
  theme(axis.title.y = element_text(size = 20, face = 'bold'), 
        axis.title.x = element_text(size = 28, face = 'bold'),
        axis.text = element_text(size = 18), 
        legend.position = 'none')
p_figs10b  

# Arrange the results by position
all_anovas_final <- all_anovas %>% 
  separate(col = Genotype, into = c('mut', 'Position'), sep = 2, remove = F) %>%
  mutate(Position = as.numeric(Position)) %>%
  arrange(Position, Genotype) %>%
  select(-mut, -Position)
write.table(x = all_anovas_final, append = F, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Figures/2022-05-09_Supp_figures_paper/TableSX.ANOVA_individual_mutants_noTMP.csv')

p_figs10 <- plot_grid(p_figs10a, p_figs10b, ncol = 2, labels = c('A', 'B'), label_size = 20, 
                      label_fontface = 'bold')

ggsave(p_figs10, device = cairo_pdf, width = 20, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS10.volcano_plot_signif_mutants_noTMP.pdf')
ggsave(p_figs10, device = 'png', width = 20, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/FigS10.volcano_plot_signif_mutants_noTMP.png')

#### Try figure 1F with the data without TMP ####

# Join the data to add these marks
data_fig1f <- left_join(x = temp_data, 
                        y= score_significant %>% select(Position, Residue, p.adj, mut_check),
                        by = c('Position' = 'Position', 'Residue' = 'Residue'))


# Add the expression level
data_fig1f_exp <- data_fig1f %>% 
  mutate(exp_level = ifelse(Arabinose == 0.01, 'Weak', 
                            ifelse(Arabinose == 0.025, 'Suboptimal', 
                                   ifelse(Arabinose == 0.05, 'Near-optimal', 
                                          ifelse(Arabinose == 0.2, 'Optimal', 
                                                 ifelse(Arabinose == 0.4, 'Overexpressed', 
                                                        NA)))))) %>%
  mutate(exp_level = factor(exp_level, levels = c('Weak', 'Suboptimal', 'Near-optimal', 
                                                  'Optimal', 'Overexpressed')))

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
  geom_violin(alpha = 0.8) +
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
ggsave(p_fig1f, device = cairo_pdf, width = 14, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/Fig1F_noTMP.pdf')
ggsave(p_fig1f, device = 'png', width = 14, height = 10, dpi = 300, 
       filename = 'Figures/2022-05-09_Supp_figures_paper/Fig1F_noTMP.png')

### Do the ANOVA on ranks for all the mutants and all replicates
data_anova <- all_data_all_reps_TMP0 %>% 
  mutate(ID = as.factor(ID),
         Arabinose = as.factor(Arabinose), 
         Genotype = as.factor(Genotype))

## Takes a long time to run
m <- art(sel_coeff ~ Arabinose*Genotype + Error(ID), data=data_anova)
summary(m)
anova_general <- anova(m)

anova_general_table <- data.frame(anova_general) %>%
  mutate(Pr..F. = ifelse(Pr..F. == 0, '<2.2e-16', Pr..F.))
write.table(x = anova_general_table, append = F, quote = F, sep = '\t', row.names = F, col.names = T,
            file = 'Figures/2022-05-09_Supp_figures_paper/TableSX.ANOVA_general_noTMP.csv')



# #### These figures will be used for Figure 5 ####
# #### Fig. S13: deltaS vs S #### (this is actually probably going to the main text)
# 
# # Mutants to annotate
# list_mut <- joined_sets %>% select(Position, WT_Residue, Residue) %>% ungroup() %>%
#   unique() %>% mutate(Genotype = str_c(WT_Residue, Position, Residue, sep = ''))
# 
# # Separate the data for ara 0.2
# data_part_1 <- data_fig_4 %>%
#   filter(Arabinose == 0.2)
# 
# # Separate the data for ara 0.01
# data_part_2 <- data_fig_4 %>%
#   filter(Arabinose == 0.01)
# 
# ## Subtract the scores (just check that the order is the same)
# # Change column names
# colnames(data_part_2) <- c("Position", "WT_Residue", "Residue", "mean_sel_coeff_2", "Arabinose_2")
# 
# 
# # Join
# data_fig_4_final <- left_join(x = data_part_1, y = data_part_2, 
#                               by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 'Residue' = 'Residue')
# ) %>% 
#   mutate(diffNormScore = mean_sel_coeff_2 - mean_sel_coeff)
# 
# mutants_highlight <- data_fig_4_final %>% rowwise() %>%
#   mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
#                             ifelse(Residue == '*', 'Stop', 'Missense')), 
#          Genotype = str_c(WT_Residue, Position, Residue, sep = '')) %>%
#   filter(Genotype %in% list_mut$Genotype)
# 
# # Rename WT
# list_mut %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))
# mutants_highlight %<>% mutate(Genotype = ifelse(Genotype == 'E2E', 'WT', Genotype))
# 
# # Get data to represent the stop and WT codons
# summary_wt_stop <- data_fig_4_final %>% rowwise() %>% ungroup() %>%
#   mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
#                             ifelse(Residue == '*', 'Stop', 'Missense'))) %>% 
#   filter(mut_check != 'Missense') %>%
#   group_by(mut_check) %>%
#   summarise(mean_s = mean(mean_sel_coeff), sem_s = sd(mean_sel_coeff) / sqrt(n()), 
#             mean_ds = mean(diffNormScore), sem_ds = sd(diffNormScore) / sqrt(n()))
# 
# data_fig5 <- data_fig_4_final %>% rowwise() %>%
#  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
#                            ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
#  filter(mut_check == 'Missense') %>% 
#  mutate(sem_s = 0, sem_ds = 0) %>% # These are single points
#  select(mut_check, mean_sel_coeff, sem_s, diffNormScore, sem_ds)
# 
# colnames(data_fig5) <- c('mut_check', 'mean_s','sem_s', 'mean_ds', 'sem_ds')
# 
# data_fig_5 <- rbind(data_fig5, summary_wt_stop)
# 
# ## Let's try a different figure: delta(DMS scores) vs DMS scores
# p <- data_fig_4_final %>% rowwise() %>%
#   mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
#                             ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
#   # filter(mut_check == 'Missense') %>%
#   ggplot(aes(x = mean_sel_coeff, y = diffNormScore)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
#   geom_point(size = 3, aes(colour = mut_check, alpha = mut_check)) +
#   # geom_point(# inherit.aes = F, 
#   #            data = summary_wt_stop, 
#   #            aes(x = mean_s, y = mean_ds), size = 5) +
#   # geom_errorbar(# inherit.aes = F,
#   #               data = summary_wt_stop, size = 1,
#   #               aes(x = mean_s, y = mean_ds,
#   #                 ymax = mean_ds + sem_ds, ymin = mean_ds - sem_ds, colour = mut_check)) +
#   # geom_errorbarh(# inherit.aes = F, 
#   #                data = summary_wt_stop, size = 1,
#   #                aes(x = mean_s, y = mean_ds, 
#   #                    xmax = mean_s + sem_s, xmin = mean_s - sem_s, colour = mut_check)) +
#   geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
#                    box.padding = 0.4) +
#   xlab('Sel. coeff. (optimal expression)') + 
#   scale_colour_manual(values = c('grey', 'red', 'blue')) +
#   scale_alpha_manual(values = c(0.3, 1, 1), guide = 'none') +
#   labs(colour = 'Mutation type', y = expression(s[weak] - s[opt])) +
#   theme(
#         axis.title = element_text(face = 'bold', size = 20), 
#         axis.text = element_text(size = 18), 
#         legend.position = 'top',
#         legend.justification = 'center', 
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 18)
#   ) + 
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
#   annotate('text', x = -0.7, y = 0.2, 
#            # label = 'Higher fitness at low expression s[weak] > s[opt]', parse = T) +
#            label = 's[weak] > s[opt]', parse = T, size = 8) +
#   annotate('text', x = -0.7, y = -0.5, 
#            # label = 'Higher fitness at optimal expression sweak < sopt', parse = T)
#            label = 's[weak] < s[opt]', parse = T, size = 8)
# 
# p
# ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_deltaDMSscore.pdf')
# 
# ## Rename the column for mutants to highlight
# mutants_highlight %<>% mutate(mean_s = mean_sel_coeff, mean_ds = diffNormScore)
# 
# # New version
# p_fig5a <- data_fig_5 %>% rowwise() %>%
#   # mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
#   #                           ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
#   mutate(mut_check = factor(mut_check, levels = c('Missense', 'Stop', 'WT'))) %>%
#   ggplot(aes(x = mean_s, y = mean_ds)) +
#   geom_hline(yintercept = 0, linetype = 'dashed', size = 1) +
#   geom_point(aes(alpha = mut_check, size = mut_check, # shape = mut_check, 
#                  colour = mut_check)) +
#   geom_errorbar(aes(ymax = mean_ds + sem_ds, ymin = mean_ds - sem_ds, 
#                     colour = mut_check)) +
#   geom_errorbarh(aes(xmax = mean_s + sem_s, xmin = mean_s - sem_s, 
#                      colour = mut_check)) +
#   geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
#                  box.padding = 0.4, size = 4, fontface = 'bold') +
#   xlab('Sel. coeff. (optimal expression)') + 
#   scale_size_manual(values = c(3, 5, 5)) +
#   guides(size = 'none', alpha = 'none') +
#   # scale_shape_manual(values = c(21, 24, 25)) +
#   scale_colour_manual(values = c('grey', 'red', 'blue')) +
#   # scale_fill_manual(values = c('grey', 'red', 'blue')) +
#   scale_alpha_manual(values = c(0.5, 1, 1)) +
#   labs(colour = 'Mutation type', y = expression(s[weak] - s[opt])) +
#   theme(
#     axis.title.x = element_text(face = 'bold', size = 20), 
#     axis.title.y = element_text(face = 'bold', size = 28),
#     axis.text = element_text(size = 18), 
#     legend.position = 'top',
#     legend.justification = 'center', 
#     legend.title = element_text(size = 20),
#     legend.text = element_text(size = 18)
#   ) + 
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4) +
#   annotate('text', x = -0.7, y = 0.2, 
#            # label = 'Higher fitness at low expression s[weak] > s[opt]', parse = T) +
#            label = 's[weak] > s[opt]', parse = T, size = 8) +
#   annotate('text', x = -0.7, y = -0.5, 
#            # label = 'Higher fitness at optimal expression sweak < sopt', parse = T)
#            label = 's[weak] < s[opt]', parse = T, size = 8)
# 
# p_fig5a
# ggsave(p_fig5a, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_deltaDMSscore.pdf')
# 
# 
# 
# ## Similar figure but coloring by secondary structure 
# data_secStruc_new <- left_join(x = data_fig_4_final, 
#                                y = all_data_complete %>% select(Position, Secondary_structure), 
#                                by = c('Position' = 'Position')) %>%
#   unique() %>% rowwise() %>%
#   mutate(Secondary_structure = ifelse(Secondary_structure == 'Missing', 'Disordered', Secondary_structure))
# 
# # Add secondary structure to the mutants that I need to highlight
# mutants_highlight_secStruc <- left_join(x = mutants_highlight,
#                                         y = data_secStruc_new %>% select(Position, WT_Residue, Residue, 
#                                                                          Secondary_structure), 
#                                         by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
#                                                'Residue' = 'Residue'))
# 
# p <- data_secStruc_new %>% rowwise() %>%
#   mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
#                             ifelse(Residue == '*', 'Stop', 'Missense')), 
#          Secondary_structure = factor(Secondary_structure, levels = c('none', '3/10 helix', 'Bend', 
#                                                                       'Beta ladder', 'Disordered', 
#                                                                       'H-bonded turn'))) %>%
#   ggplot(aes(x = mean_sel_coeff, y = diffNormScore)) +
#   # geom_point(size = 3, aes(colour = Secondary_structure)) +
#   geom_hex() +
#   facet_wrap(~Secondary_structure, nrow = 2) +
#   geom_label_repel(data = mutants_highlight_secStruc, aes(label = Genotype, group = Secondary_structure), 
#                    box.padding = 0.4) +
#   xlab('s (optimal expression)') + ylab('\u0394s') + 
#   scale_colour_manual(values = c('#000000', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7')) +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         axis.title = element_text(face = 'bold', size = 20), 
#         axis.text = element_text(size = 18), 
#         # legend.position = 'top',
#         # legend.justification = 'center', 
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 18)
#   ) + 
#   # labs(colour = 'Secondary structure') +
#   labs(fill = 'Count') +
#   scale_fill_gradientn(colours = terrain.colors(10)) +
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4)
# p
# ggsave(p, device = cairo_pdf, width = 17, height = 10, dpi = 300, 
#        filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_deltaDMSscore_secStruc.pdf')
# 
# 
# ## Similar figure but colored by RSA class ##
# data_RSA_new <- left_join(x = data_fig_4_final, 
#                                y = all_data_complete %>% select(Position, rSASA), 
#                                by = c('Position' = 'Position')) %>%
#   unique() %>% rowwise() %>%
#   mutate(RSA_class = ifelse(Position <= 21, 'Unavailable',
#                             ifelse(rSASA <= 0.25, 'Buried', 'Exposed')))
# 
# p <- data_RSA_new %>% rowwise() %>%
#   mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
#                             ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
#   ggplot(aes(x = mean_sel_coeff, y = diffNormScore)) +
#   geom_point(size = 3, aes(colour = RSA_class)) +
#   geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
#                    box.padding = 0.4) +
#   xlab('s (optimal expression)') + ylab('\u0394s') + 
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         axis.title = element_text(face = 'bold', size = 20), 
#         axis.text = element_text(size = 18), 
#         legend.position = 'top',
#         legend.justification = 'center', 
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 18)
#   ) + 
#   labs(colour = 'Solvent accessibility') +
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4)
# p
# ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
#        filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_deltaDMSscore_RSA.pdf')
# 
# 
# 
# ## Similar figure but colored by protein regions ##
# data_regions <- read_delim('Data/Complete_datasets_TMP0_TMP10/data_annotation_2.txt', delim = '\t') %>%
#   mutate(Disordered_region = ifelse(Position %in% c(2:21), 1, 0)) %>% rowwise() %>%
#   mutate(Unannotated = ifelse(sum(DHF, NADPH, `A,C`, `A,D`, Disordered_region, Cat_residues) == 0, 1, 0)) %>%
#   pivot_longer(cols = c('DHF', 'NADPH', 'A,C', 'A,D', 'Disordered_region', 'Cat_residues', 'Unannotated'), 
#                names_to = 'Region', values_to = 'Value') %>%
#   filter(Value == 1)
#   # mutate(Region = ifelse(DHF == 1, 'DHF', 
#   #                        ifelse(NADPH == 1, 'NADPH', 
#   #                               ifelse(`A,C` == 1, 'Interface 1', 
#   #                                      ifelse(`A,D` == 1, 'Interface 2', 
#   #                                             ifelse(Disordered_region == 1, 'Disordered region', 
#   #                                                    'Unannotated'))))))
# 
# data_regions_new <- left_join(x = data_fig_4_final, 
#                           y = data_regions %>% select(Position, Region),
#                           by = c('Position' = 'Position')) %>%
#   unique()
# 
# # Add regions to highlighted mutants
# mutants_highlight_regions <- left_join(x = mutants_highlight,
#                                         y = data_regions_new %>% select(Position, WT_Residue, Residue, 
#                                                                          Region), 
#                                         by = c('Position' = 'Position', 'WT_Residue' = 'WT_Residue', 
#                                                'Residue' = 'Residue'))  %>%
#   mutate(Region = ifelse(Region == 'A,C', 'Interface 1', 
#                          ifelse(Region == 'A,D', 'Interface 2', 
#                                 ifelse(Region == 'Disordered_region', 'Disordered region',
#                                        ifelse(Region == 'Cat_residues', 'Catalytic residues', Region)))))
# 
# p_fig5b <- data_regions_new %>% rowwise() %>%
#   mutate(Region = ifelse(Region == 'A,C', 'Interface 1', 
#                          ifelse(Region == 'A,D', 'Interface 2', 
#                                 ifelse(Region == 'Disordered_region', 'Disordered region',
#                                        ifelse(Region == 'Cat_residues', 'Catalytic residues', Region))))) %>% 
#   # If I don't want to include unannotated
#   filter(Region != 'Unannotated') %>%
#   mutate(Region = factor(Region, levels = c('Catalytic residues', 'DHF', 'Interface 1', 
#                                             'Interface 2', 'NADPH', 'Disordered region' 
#   ))) %>%
#   ggplot(aes(x = mean_sel_coeff, y = diffNormScore)) +
#   geom_point(size = 3, aes(colour = Region)) +
#   # geom_hex() +
#   geom_label_repel(data = mutants_highlight_regions %>% 
#                      # When I don't want to include the unannotated
#                      filter(Region != 'Unannotated') %>%
#                      mutate(Region = factor(Region, levels = c('Catalytic residues', 'DHF', 'Interface 1', 
#                                                                'Interface 2', 'NADPH', 'Disordered region' 
#                                                                ))),
#                    aes(label = Genotype, group = Region), 
#                    box.padding = 0.4, size = 4, fontface = 'bold') +
#   facet_wrap(~Region, nrow = 3, scales = 'free') +
#   xlab('s (optimal expression)') +
#   labs(y = expression(s[weak] - s[opt])) +
#   scale_colour_manual(values = c('#56B4E9', '#E69F00', '#009E73', '#0072B2', '#D55E00', '#CC79A7'# ,
#                                  # '#000000'
#                                  )) +
#   theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
#         panel.grid.minor = element_blank(), 
#         axis.title.x = element_text(face = 'bold', size = 20), 
#         axis.title.y = element_text(face = 'bold', size = 28),
#         axis.text = element_text(size = 18), 
#         legend.position = 'none',
#         legend.justification = 'center', 
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 18), 
#         strip.text = element_text(size = 20, face = 'bold'), 
#         strip.background = element_rect(fill = 'white'), 
#         axis.line=element_line()
#   ) + 
#   labs(colour = 'Site') +
#   xlim(-0.9, 0.4) + ylim(-0.9, 0.4)
# p_fig5b
# ggsave(p_fig5b, device = cairo_pdf, width = 14, height = 21, dpi = 300, 
#        filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_deltaDMSscore_regions.pdf')
# 
# #### Put the figures together ####
# 
# p_fig5 <- plot_grid(p_fig5a, p_fig5b, ncol = 2, labels = c('A', 'B'), label_fontface = 'bold')
# ggsave(p_fig5, device = cairo_pdf, width = 20, height = 10, dpi = 300, 
#        filename = 'Figures/2022-04-06_Figures_paper/5.Fig5_no_unannotated.pdf')

#### Stuff I am not using anymore ####

## Alternatives to the figure with deltaS vs S
p <- data_fig_4_final %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
  ggplot(aes(x = mean_sel_coeff, y = mean_sel_coeff_2)) +
  geom_point(size = 3, aes(colour = mut_check)) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4) +
  xlab('s (optimal expression)') + ylab('s (weak expression)') + 
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.position = 'top',
        legend.justification = 'center', 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)
  ) + 
  labs(colour = 'Mutation type') +
  xlim(-1, 0.5) + ylim(-1, 0.5)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_DMSscore.pdf')


p <- data_secStruc_new %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense')), 
         Secondary_structure = factor(Secondary_structure, levels = c('none', '3/10 helix', 'Bend', 
                                                                      'Beta ladder', 'Disordered', 
                                                                      'H-bonded turn'))) %>%
  ggplot(aes(x = mean_sel_coeff, y = mean_sel_coeff_2)) +
  # geom_point(size = 3, aes(colour = Secondary_structure)) +
  geom_hex() +
  facet_wrap(~Secondary_structure, nrow = 2) +
  geom_label_repel(data = mutants_highlight_secStruc, aes(label = Genotype, group = Secondary_structure), 
                   box.padding = 0.4) +
  xlab('s (optimal expression)') + ylab('s (weak expression)') + 
  scale_colour_manual(values = c('#000000', '#56B4E9', '#009E73', '#0072B2', '#D55E00', '#CC79A7')) +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        # legend.position = 'top',
        # legend.justification = 'center', 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)
  ) + 
  # labs(colour = 'Secondary structure') +
  labs(fill = 'Count') +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  xlim(-1, 0.5) + ylim(-1, 0.5)
p
ggsave(p, device = cairo_pdf, width = 17, height = 10, dpi = 300, 
       filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_DMSscore_secStruc.pdf')


p <- data_RSA_new %>% rowwise() %>%
  mutate(mut_check = ifelse(Residue == WT_Residue, 'WT',
                            ifelse(Residue == '*', 'Stop', 'Missense'))) %>%
  ggplot(aes(x = mean_sel_coeff, y = mean_sel_coeff_2)) +
  geom_point(size = 3, aes(colour = RSA_class)) +
  geom_label_repel(data = mutants_highlight, aes(label = Genotype), 
                   box.padding = 0.4) +
  xlab('s (optimal expression)') + ylab('s (weak expression)') + 
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.position = 'top',
        legend.justification = 'center', 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)
  ) + 
  labs(colour = 'Solvent accessibility') +
  xlim(-1, 0.5) + ylim(-1, 0.5)
p
ggsave(p, device = cairo_pdf, width = 10, height = 7, dpi = 300, 
       filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_DMSscore_RSA.pdf')

p <- data_regions_new %>% rowwise() %>%
  mutate(Region = ifelse(Region == 'A,C', 'Interface 1', 
                         ifelse(Region == 'A,D', 'Interface 2', 
                                ifelse(Region == 'Disordered_region', 'Disordered region',
                                       ifelse(Region == 'Cat_residues', 'Catalytic residues', Region))))) %>% 
  ggplot(aes(x = mean_sel_coeff, y = mean_sel_coeff_2)) +
  geom_point(size = 3, aes(colour = Region)) +
  geom_label_repel(data = mutants_highlight_regions, aes(label = Genotype, group = Region), 
                   box.padding = 0.4) +
  facet_wrap(~Region, nrow = 2) +
  xlab('s (optimal expression)') + ylab('s (weak expression)') + 
  scale_colour_manual(values = c('#56B4E9', '#E69F00', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#000000')) +
  theme(panel.grid.major = element_line(colour="#8c8c8c", linetype = 'dashed'),
        panel.grid.minor = element_blank(), 
        axis.title = element_text(face = 'bold', size = 20), 
        axis.text = element_text(size = 18), 
        legend.position = 'top',
        legend.justification = 'center', 
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)
  ) + 
  labs(colour = 'Site') +
  xlim(-1, 0.5) + ylim(-1, 0.5)
p
ggsave(p, device = cairo_pdf, width = 24, height = 14, dpi = 300, 
       filename = 'Figures/2022-04-07_supp_figures/FigS13.DMSscore_vs_DMSscore_regions.pdf')
