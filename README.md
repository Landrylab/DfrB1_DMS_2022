Scripts and data used for "Epistasis between promoter activity and coding mutations shapes gene evolvability" by Angel F. Cisneros*, Isabelle Gagnon-Arsenault*, Alexandre K. Dubé, Philippe C. Després, Pradum Kumar, Christian R. Landry

\* indicates equal contributions

Overview of each folder:

Scripts:

- Sequencing data processing: These scripts demultiplex the sequencing data, merge reads, and align to the reference sequence. They populate the folders at Data/MiSeq_analysis and Data/NovaSeq_analysis with the intermediate files of the preprocessing steps. The final script generates the complete dataset used to produce the figures.

- In_silico_mutagenesis: Preprocess PDB file 2rk1 and perform random mutagenesis using mutatex.

- Random_forest: Script used to train and validate the random forest regressor for expression-dependent differences in fitness effects.

- Evolutionary_analysis: Script used to calculate the Shannon entropy based on the multiple sequence alignment.

- Structural_analysis: Scripts used to organize AlphaFold2 pLDDT values, DSSP data, and prepare PDB files for visualization.

- Scripts_for_figures: Scripts used to generate the main and supplementary figures shown in the paper.

Data: Contains data used for all the experimental assays and computational analyses described in the paper.

Figures: PDF versions of main and supplementary figures shown in the paper.

To reanalyze the data, first download the sequencing data (accession numbers SRR19419448 and SRR19419449) in the corresponding folder (Data/Sequencing_data).


