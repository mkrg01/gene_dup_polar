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
    add_row(sci_name = 'Dissostichus_mawsoni', polar = TRUE, lat_min_fb = -78, lat_max_fb = -45, lat_abs_max = 78) %>% # https://fishbase.se/summary/Dissostichus-mawsoni.html
    write_tsv(args$output)