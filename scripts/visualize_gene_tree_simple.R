library(tidyverse)
library(argparse)
library(ape)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

# Parse
parser <- ArgumentParser()
parser$add_argument('--seqinfo', type='character')
parser$add_argument('--tree', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
tree <- read.tree(args$tree)
seqinfo_df <- read_tsv(args$seqinfo)

# Visualize
highlight <- tree$tip.label[grep("^Dissostichus_mawsoni", tree$tip.label)]
print(highlight)
polar_fish <- seqinfo_df %>%
    filter(polar == TRUE) %>%
    pull(seq_id)

p <- ggtree(tree) +
	geom_tiplab(aes(label=label, color=ifelse(label %in% highlight, "red", ifelse(label %in% polar_fish, "cyan3", "black"))), size=1.0) +
    geom_treescale(x=0.1, y=length(seqinfo_df$seq_id)*0.55) + 
    hexpand(0.03) +
    scale_color_identity()

ggsave(args$output, plot=p, width=8, height=11)