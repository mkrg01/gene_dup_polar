library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--select1', type='character')
parser$add_argument('--select2', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
df1 <- read_tsv(args$select1)
df2 <- read_tsv(args$select2)

# Filter
OG_list_1 <- df1 %>%
	filter(fold > 2) %>%
	filter(n_sp_with_more_than_one_gene_polar > n_sp_polar * 0.5) %>% 
    .$OG_id
OG_list_2 <- df2 %>%
	filter(fold > 2) %>%
	filter(n_sp_with_more_than_one_gene_polar > n_sp_polar * 0.5) %>%
    .$OG_id

# Merge
df1 <- df1 %>%
    dplyr::select(OG_id, OG_name, fold) %>%
    rename(fold_all = fold)
df2 <- df2 %>%
    dplyr::select(OG_id, fold) %>%
    rename(fold_downsampling = fold)
df <- df1 %>% 
    left_join(df2, by='OG_id')

# Output
df %>% 
    filter(OG_id %in% OG_list_1) %>%
    arrange(desc(fold_all)) %>% 
    write_tsv(args$output)