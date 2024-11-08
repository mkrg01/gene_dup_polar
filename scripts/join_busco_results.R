library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--busco', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
busco_df <- read_tsv(args$busco)

busco_df <- busco_df %>%
    mutate(Input_file = str_replace(Input_file, '.fa$', '')) %>%
    rename(organism_id = Input_file, busco_dataset = Dataset, complete = Complete, single = Single, duplicated = Duplicated, fragmented = Fragmented, missing = Missing) %>%
    select(organism_id, busco_dataset, complete, single, duplicated, fragmented, missing)

# Join
sp_list_df %>%
    left_join(busco_df, by='organism_id') %>%
	write_tsv(args$output)