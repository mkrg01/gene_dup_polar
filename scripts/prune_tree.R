library(tidyverse)
library(argparse)
library(ape)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--tree', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
tree <- read.tree(args$tree)

# Prune
tree <- keep.tip(tree, sp_list_df$sci_name)
write.tree(tree, args$output)