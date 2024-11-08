library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--n_gene_by_OG', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
n_gene_by_OG_df <- read_tsv(args$n_gene_by_OG)

polar_sp <- sp_list_df %>%
	filter(polar == 'TRUE') %>%
	.$organism_id
nonpolar_sp <- sp_list_df %>%
	filter(polar == 'FALSE') %>%
	.$organism_id
n_gene_by_OG_df <- n_gene_by_OG_df %>%
	mutate(n_gene_polar = rowSums(select(., polar_sp))) %>%
	mutate(n_gene_nonpolar = rowSums(select(., nonpolar_sp))) %>%
	mutate(n_sp_polar = length(polar_sp)) %>%
	mutate(n_sp_nonpolar = length(nonpolar_sp))

gene_more_than_one_df <- n_gene_by_OG_df %>%
	mutate(across(polar_sp, ~ .x > 1)) %>%
	mutate(across(nonpolar_sp, ~ .x > 1)) %>%
	mutate(n_sp_with_more_than_one_gene_polar = rowSums(select(., polar_sp))) %>%
	mutate(n_sp_with_more_than_one_gene_nonpolar = rowSums(select(., nonpolar_sp))) %>%
	select(OG_id, n_sp_with_more_than_one_gene_polar, n_sp_with_more_than_one_gene_nonpolar)

n_gene_by_OG_df %>%
	left_join(gene_more_than_one_df, by='OG_id') %>%
	select(OG_id, OG_name, n_gene_polar, n_gene_nonpolar, n_sp_with_more_than_one_gene_polar, n_sp_with_more_than_one_gene_nonpolar, n_sp_polar, n_sp_nonpolar) %>%
	mutate(fold = (n_gene_polar / n_sp_polar) / (n_gene_nonpolar / n_sp_nonpolar)) %>%
	arrange(desc(fold)) %>%
	write_tsv(args$output)