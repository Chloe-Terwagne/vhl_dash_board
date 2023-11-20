# Run it using: pymol -cq isolate_chains.pml

# Load the PDB file
fetch 1lm8

# Remove HETAM from chain V
select HETATM_V, chain V and hetatm
remove HETATM_V

# Remove all atoms except chain V and renumber
select V, chain V
save 1LM8_vhl_isolated.pdb, V
delete all

# Remove all atoms except chains V, C, and H and renumber --------------------------------# Load the PDB file
fetch 1lm8

# Remove HETAM from chain V
select HETATM_V, chain V and hetatm
remove HETATM_V

# Remove all atoms except chains V, C, and H and renumber
set retain_order
select V, chain V
alter V, rank=-1
sort
select VCH, chain V+C+H
save 1LM8_vch_isolated.pdb, VCH
delete all

# Remove all atoms except chains H ------------------------------------------------------------------
fetch 1lm8

# Select HIF
select H, chain H
save 1LM8_h_isolated.pdb, H

# Remove all atoms except chains H ------------------------------------------------------------------
fetch 1lm8

# Select ELOC-C
select C, chain C
save 1LM8_C_isolated.pdb, C