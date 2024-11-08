library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--latitude', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
latitude_df <- read_tsv(args$latitude)

# Join
sp_list_df <- sp_list_df %>%
    left_join(latitude_df, by='sci_name') %>%
	mutate(lat_abs_max = pmax(abs(lat_max_fb), abs(lat_min_fb))) %>%
	write_tsv(args$output)