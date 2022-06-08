#!/bin/bash

# This code will receive an input PDB file for mutagenesis and the option
# to do it considering multimers or not. Input arguments are:
# $1 = path to the input PDB (without the pdb extension)
# $2 = path to the output directory
# $3 = multimer or not (0 for no multimers, 1 for multimers)
# $4 = number of cores to use
# $5 = amount of memory to request (Mb)

# Organize the output directory
mkdir -p $2
cp $1.pdb $2

cp mutation_list.txt $2
cd $2


prot=$(basename $1)

# Set the flag for multimers or no multimers
if [ $3 -eq 0 ]
then
    mut_arg='--no-multimers'
fi

cat > mut_$3_${prot}.sbatch << EOF
#!/bin/bash

#SBATCH -D $PWD
#SBATCH -J mut_$3_${prot}
#SBATCH -o mut_$3_${prot}.out
#SBATCH -c $4
#SBATCH -p medium
#SBATCH --time=2-00:00
#SBATCH --mem=$5

unset $PYTHONPATH

# Load virtual environment if necessary
source /path/to/virtual/env

mutatex ${prot}.pdb \
        -p $4 \
        -m mutation_list.txt \
        -x $FOLDX_BINARY \
        -f suite5 \
        -R /path/to/repair_runfile_template.txt \
        -M /path/to/mutate_runfile_template.txt \
        -I /path/to/interface_runfile_template.txt \
        -B $mut_arg


EOF

sbatch mut_$3_${prot}.sbatch

