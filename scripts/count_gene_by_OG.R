library(tidyverse)
library(argparse)
library(ape)

# Parse
parser <- ArgumentParser()
parser$add_argument('--input', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
gene_df <- read_tsv(args$input)

# Summarise the number of genes by OG by species (row: OG, column: species)
n_gene_by_OG <- gene_df %>%
	select(OG_id, organism_id) %>%
	group_by(OG_id, organism_id) %>%
	summarise(num=n()) %>%
	ungroup() %>%
	filter(!is.na(OG_id)) %>%
	pivot_wider(names_from = organism_id, values_from = num, values_fill = 0)

# Add OG_name
gene_df %>%
	select(OG_id, OG_name) %>%
	distinct() %>%
	filter(!is.na(OG_id)) %>%
	left_join(n_gene_by_OG, by='OG_id') %>%
	write_tsv(args$output)