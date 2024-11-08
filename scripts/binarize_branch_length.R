library(tidyverse)
library(ape)
library(argparse)

# Input
parser <- ArgumentParser(description='')
parser$add_argument('-t', '--tree', help='input a tree')
parser$add_argument('-o', '--output', help='output directory', default='output_tree')
args = parser$parse_args()

# Load a tree
tree <- ape::read.tree(args$tree)

# Set all branch lengths to 1
tree$edge.length <- rep(1, length(tree$edge.length))

# Write the tree in a newick format
ape::write.tree(tree, file=args$output)
