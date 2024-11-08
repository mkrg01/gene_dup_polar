library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--genes', type='character')
parser$add_argument('--OGs', type='character')
parser$add_argument('--OG2genes', type='character')
parser$add_argument('--output', type='character')
# parser$add_argument('--output_ver', type='character')
# parser$add_argument('--output_met', type='character')
# parser$add_argument('--output_euk', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
genes_df <- read_tsv(args$genes, col_names=FALSE)
OGs_df <- read_tsv(args$OGs, col_names=FALSE)
OG2genes_df <- read_tsv(args$OG2genes, col_names=FALSE)

gene_orgid_df <- genes_df %>%
	filter(X2 %in% sp_list_df$organism_id) %>%
	select(X1, X2) %>%
	rename(gene_id=X1, organism_id=X2)
OGs_df <- OGs_df %>%
	rename(OG_id=X1, OG_level=X2, OG_name=X3)
OG2genes_df <- OG2genes_df %>%
	rename(OG_id=X1, gene_id=X2)

# Actinopterygii OG (tax id: 7898)
OGs_act_df <- OGs_df %>%
	filter(OG_level==7898)
OG2genes_act_df <- OG2genes_df %>%
	inner_join(OGs_act_df, by='OG_id')
gene_orgid_df %>%
	left_join(OG2genes_act_df, by='gene_id') %>%
	select(gene_id, organism_id, OG_id, OG_name) %>%
	write_tsv(args$output)

# # Vertebrata OG (tax id: 7742)
# OGs_ver_df <- OGs_df %>%
# 	filter(OG_level==7742)
# OG2genes_ver_df <- OG2genes_df %>%
# 	inner_join(OGs_ver_df, by='OG_id')
# gene_orgid_df %>%
# 	left_join(OG2genes_ver_df, by='gene_id') %>%
#	select(gene_id, organism_id, OG_id, OG_name) %>%
# 	write_tsv(args$output_ver)

# # Metazoa OG (tax id: 33208)
# OGs_met_df <- OGs_df %>%
# 	filter(OG_level==33208)
# OG2genes_met_df <- OG2genes_df %>%
# 	inner_join(OGs_met_df, by='OG_id')
# gene_orgid_df %>%
# 	left_join(OG2genes_met_df, by='gene_id') %>%
#	select(gene_id, organism_id, OG_id, OG_name) %>%
# 	write_tsv(args$output_met)

# # Eukaryota OG (tax id: 2759)
# OGs_euk_df <- OGs_df %>%
# 	filter(OG_level==2759)
# OG2genes_euk_df <- OG2genes_df %>%
# 	inner_join(OGs_euk_df, by='OG_id')
# gene_orgid_df %>%
# 	left_join(OG2genes_euk_df, by='gene_id') %>%
#	select(gene_id, organism_id, OG_id, OG_name) %>%
# 	write_tsv(args$output_euk)