This folder contains the scripts used to run in-silico mutagenesis with MutateX using a SLURM queue manager. 

Scripts must be executed sequentially specifying the arguments mentioned at the top of each script to generate
the results in the Data/In_silico_mutagenesis_results/006_gathered_data folder. call_interfaces_helper.py is
called automatically by 004_call_interfaces to load some required functions.

The only parameters that should be adjusted are the command line arguments used when executing each script, the 
SLURM job parameters in scripts 003_foldx_repair_slurm_multiple and 005_mutatex_slurm.sh, a virtual python
environment with mutatex (if necessary), and the paths to mutatex templates in 005_mutatex_slurm.sh
