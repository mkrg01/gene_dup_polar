library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--polar_fish', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
polar_fish_df <- read_csv(args$polar_fish)

# Join
sp_list_df <- sp_list_df %>%
    left_join(polar_fish_df, by='sci_name') %>%
	write_tsv(args$output)