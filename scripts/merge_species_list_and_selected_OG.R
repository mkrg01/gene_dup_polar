library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--gene_list', type='character')
parser$add_argument('--selected_OG', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_df <- read_tsv(args$sp_list)
gene_df <- read_tsv(args$gene_list)
selected_OG_df <- read_tsv(args$selected_OG)

# Summarise the number of genes by OG by species (row: organism_id, column: OG_id)
selected_organism_id_list <- sp_df$organism_id
selected_OG_list <- selected_OG_df$OG_id
n_gene_by_OG_df <- gene_df %>%
	select(OG_id, organism_id) %>%
    filter(OG_id %in% selected_OG_list) %>%
    filter(organism_id %in% selected_organism_id_list) %>%
	group_by(OG_id, organism_id) %>%
	summarise(num=n()) %>%
	ungroup() %>%
	pivot_wider(names_from = OG_id, values_from = num, values_fill = 0) %>%
    select(organism_id, selected_OG_list)

# Merge and output
sp_df %>%
    left_join(n_gene_by_OG_df, by='organism_id') %>%
    write_tsv(args$output)