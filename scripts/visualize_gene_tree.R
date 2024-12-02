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
	geom_tiplab(aes(label=label, color=ifelse(label %in% highlight, "#FF0800", ifelse(label %in% polar_fish, "#00AAFF", "black"))), size=1.0) +
    geom_treescale(x=0, y=length(seqinfo_df$seq_id)*0.85) + 
    hexpand(0.03) +
    scale_color_identity() +
    geom_point2(aes(subset=!is.na(as.numeric(label)), fill=cut(as.numeric(label), c(0, 70, 90, 100))), shape=21, size=0.8) +
    theme_tree(legend.position=c(0.10, 0.95)) + 
    scale_fill_manual(values=c("white", "grey", "black"), guide='legend', name='Bootstrap support (BS)', breaks=c('(90,100]', '(70,90]', '(0,70]'), labels=expression(BS>=90,70 <= BS * " < 90", BS < 70))

ggsave(args$output, plot=p, width=8, height=11)