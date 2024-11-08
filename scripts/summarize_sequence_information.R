library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--DE_ratio', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
DE_ratio_df <- read_tsv(args$DE_ratio)

seqinfo_df <- DE_ratio_df %>%
	left_join(sp_list_df, by = 'sci_name') %>%
    write_tsv(args$output)