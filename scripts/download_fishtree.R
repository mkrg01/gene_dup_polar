library(tidyverse)
library(fishtree)
library(ape)

# Get a phylogenetic tree for ray-finned fishes by fishtree
tree <- fishtree_phylogeny()

# Save
write.tree(tree, file='data/raw_data/tree/fishtree.nwk')