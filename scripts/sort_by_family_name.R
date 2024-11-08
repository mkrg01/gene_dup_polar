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
sp_list_df %>%
    arrange(family) %>%
    write_tsv(args$output)