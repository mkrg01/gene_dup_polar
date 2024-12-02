library(tidyverse)
library(argparse)
library(ape)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--gene', type='character')
parser$add_argument('--selected_OG', type='character')
parser$add_argument('--tree', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
gene_df <- read_tsv(args$gene)
selected_OG_df <- read_csv(args$selected_OG)
tree <- read.tree(args$tree)

# Create a heatmap (row: species, column: OG)
sp_OG <- gene_df %>%
	select(OG_id, organism_id) %>%
	group_by(OG_id, organism_id) %>%
	summarise(num=n()) %>%
	ungroup() %>%
	filter(!is.na(OG_id)) %>%
	pivot_wider(names_from = OG_id, values_from = num, values_fill = 0)

sp_OG <- sp_OG %>%
	inner_join(sp_list_df, by = 'organism_id') %>%
	select(sci_name, all_of(selected_OG_df$OG_id)) %>%
	pivot_longer(!sci_name, names_to = 'OG_id', values_to = 'n_gene') %>%
	mutate(OG_id = factor(OG_id, levels=selected_OG_df$OG_id))

sp_OG <- sp_OG %>%
	left_join(selected_OG_df, by='OG_id') %>%
	select(sci_name, OG_name, n_gene) %>%
	mutate(OG_name = factor(OG_name, levels=selected_OG_df$OG_name))

# Convert to category
gene_category_df <- sp_OG %>%
	mutate(n_gene = case_when(n_gene > 4 ~ '>4', n_gene == 4 ~ '4', n_gene == 3 ~ '3', n_gene == 2 ~ '2', n_gene == 1 ~ '1', n_gene == 0 ~ '0')) %>%
	mutate(n_gene = factor(n_gene, levels=c('>4', '4', '3', '2', '1', '0')))

# Visualize
tree_df <- as_tibble(tree)
tree_trait_df <- full_join(tree_df, sp_list_df, by=c('label'='sci_name'))
tree_trait <- as.treedata(tree_trait_df)
polar_fish <- sp_list_df %>%
    filter(polar == TRUE) %>%
    pull(sci_name)
p <- ggtree(tree_trait) +
	geom_tiplab(aes(label=paste0('italic(', label, ')~italic()'), color=ifelse(label %in% polar_fish, "#00AAFF", "black")), parse=TRUE, size=1.0) +
	scale_color_identity()
p <- p +
	new_scale_fill() +
	geom_fruit(data=gene_category_df, geom=geom_tile, mapping=aes(y=sci_name, x=OG_name, fill=n_gene), offset=0.4, pwidth=2.5, width=8, axis.params=c(axis='x', text.angle=90, text.size=1, hjust=1, line.color='white')) +
	scale_fill_manual(name='Copy number', values=c('>4'='#000000', '4'='#334E5C', '3'='#668FA3', '2'='#99C2D6', '1'='#CCE7F5', '0'='#FFFFFF')) +
	ylim(-length(tree$tip.label)/4, length(tree$tip.label)+1)

ggsave(args$output, plot=p, width=6, height=6.5)