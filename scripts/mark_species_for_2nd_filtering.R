library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--input', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$input)

selected_df <- sp_list_df %>%
	filter(select_1st == 'TRUE')

# 4. Downsampling
# If only one trait exists in a family, select species with the highest BUSCO complete score.
# If different traits exist in a family, select species with the highest BUSCO complete score for each trait.
selected_df <- selected_df %>%
	group_by(family) %>%
	filter(if (any(polar)) polar else TRUE) %>%
	slice_max(complete, n=1, with_ties = FALSE) %>%
	ungroup()

sp_list_df %>%
	mutate(select_2nd = organism_id %in% selected_df$organism_id) %>%
    arrange(sci_name) %>%
	write_tsv(args$output)