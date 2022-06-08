This folder contains scripts for structural analyses done on PDB: 2rk1. Scripts are provided both in the form of 
Jupyter notebooks and standalone Python scripts.

- AlphaFold2_pLDDT_values: takes the best AF2 model generated for DfrB1 and organizes the pLDDT values in a table.

- extract_DSSP_data: reads the results of the DSSP analysis and organizes them in a table.

- Replace_b_factors: Reads the tables with the maximum, minimum, and mean deltaS_weak values and uses them to replace
the b-factors in the PDB file for the DfrB1 structure.
