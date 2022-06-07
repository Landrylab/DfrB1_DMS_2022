Script used to visualize the deltaS values on the protein structure with ChimeraX after using the replace_b_factor_normscores script to replace the b factors with deltaS values.


## Load the PDB file and adjust view of the structure

# Set theme
preset pub

# Make ribbons wider
cartoon style thickness 1.2

## Load file with DHF + NADPH (Open Menu)

### Add color
## For files with the maximum deltas
color by bfactor range 0, 0.4  palette #F5F5F5:#F6E8C3:#D8B365:#8C510A key true

## For files with the mean deltas
color by bfactor range -0.4, 0.4  palette #01665E:#5AB4AC:#C7EAE5:#F5F5F5:#F6E8C3:#D8B365:#8C510A key true

# For files with the minimum s
color by bfactor range -0.4, 0  palette #01665E:#5AB4AC:#C7EAE5:#F5F5F5 key true

# Slightly longer legend for the figures with Min. s
key size 0.35, 0.04
key pos 0.4, 0.1
key fontSize 40
key ticks true

## Add legend title
# For files with the maximum deltas
2dlab text 'Max. \u0394  ' size 48 x 0.47 y 0.15 bold true
2dlab text 's  ' size 48 x 0.63 y 0.15 bold true italic true

# For files with the minimum deltas
2dlab text 'Min. \u0394  ' size 48 x 0.47 y 0.15 bold true
2dlab text 's  ' size 48 x 0.62025 y 0.15 bold true italic true

# For files with the mean deltas
2dlab text 'Mean \u0394  ' size 48 x 0.47 y 0.15 bold true
2dlab text 's  ' size 48 x 0.646 y 0.15 bold true italic true

## Color ligand
# Set color to #FF7A00 (click on the model color)

# Color ligand by element (#2 depends on the order in which the files were loaded)
color #2 by het 


