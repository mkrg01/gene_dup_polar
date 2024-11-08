library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--GO', type='character')
parser$add_argument('--OG_list', type='character')
parser$add_argument('--selected_OG', type='character')
parser$add_argument('--GO_category', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
GO_df <- read_tsv(args$GO, col_names=FALSE) %>%
    rename(OG_id=X1, GO_category=X2, GO_id=X3, n_gene=X4)
OG_df <- read_tsv(args$OG_list)
selected_OG_df <- read_tsv(args$selected_OG)

# Select only OG that at least one species after downsampling have
OG_list <- OG_df %>%
    filter(n_gene_polar + n_gene_nonpolar > 0) %>%
    .$OG_id
# Select only OG with GO term
GO_df <- GO_df %>%
    filter(OG_id %in% OG_list)

# Hypergeometric test
n_all_OG <- length(OG_list)
n_selected_OG <- length(selected_OG_df$OG_id)
GO_mf_bp_cc <- GO_df %>%
    filter(GO_category == args$GO_category)
GO_mf_bp_cc_count <- GO_mf_bp_cc %>%
    select(GO_id) %>%
    group_by(GO_id) %>%
    summarise(n_all_OG_per_GO=n()) %>%
    ungroup()
GO_mf_bp_cc_count_selected_OG <- GO_mf_bp_cc %>%
    filter(OG_id %in% selected_OG_df$OG_id) %>%
    select(GO_id, OG_id) %>%
    group_by(GO_id) %>%
    summarise(n_selected_OG_per_GO=n(), OG_id=paste(OG_id, collapse=', ')) %>%
    ungroup()
GO_mf_bp_cc_summary <- GO_mf_bp_cc_count_selected_OG %>%
    left_join(GO_mf_bp_cc_count, by = 'GO_id') %>%
    mutate(p_value = phyper(n_selected_OG_per_GO, n_all_OG_per_GO, n_all_OG - n_all_OG_per_GO, n_selected_OG, lower.tail=FALSE)) %>%
    mutate(adjusted_p_value = p.adjust(p_value, method = 'BH')) %>%
    arrange(desc(n_selected_OG_per_GO), p_value) %>%
    select(GO_id, n_selected_OG_per_GO, n_all_OG_per_GO, p_value, adjusted_p_value, OG_id) %>%
    write_tsv(args$output)