library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--input', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$input)

# Filter
# 1. Select species in Teleostei
# 2. Select species in FishTree
# 3. Select species with max latitude
selected_df <- sp_list_df %>%
	filter(infraclass == 'Teleostei') %>%
	filter(tree == 'TRUE') %>%
	filter(!is.na(lat_abs_max))

# Add information for filtering
sp_list_df %>%
	mutate(select_1st = organism_id %in% selected_df$organism_id) %>%
	write_tsv(args$output)