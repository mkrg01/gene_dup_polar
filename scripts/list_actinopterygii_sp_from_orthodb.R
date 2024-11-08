library(tidyverse)
library(argparse)
library(ape)
library(taxize)

# Parse
parser <- ArgumentParser()
parser$add_argument('--levels', type='character')
parser$add_argument('--species', type='character')
parser$add_argument('--level2species', type='character')
parser$add_argument('--tree', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
levels_df <- read_tsv(args$levels, col_names=FALSE)
species_df <- read_tsv(args$species, col_names=FALSE)
level2species_df <- read_tsv(args$level2species, col_names=FALSE)
tree <- read.tree(args$tree)

# Select Actinopterygii
act_id <- levels_df %>%
	filter(X2 == 'Actinopterygii') %>%
	pull(X1)
organism_id_list <- level2species_df %>%
	filter(grepl(paste0(',', act_id, ','), X4)) %>%
	pull(X2)
sp_list_df <- species_df %>%
	filter(X2 %in% organism_id_list) %>%
	select(X1, X2, X3, X4) %>%
	rename(tax_id=X1, organism_id=X2, sp=X3, assembly_id=X4)

# Format species names
sp_list_df$sp <- gsub(' ', '_', sp_list_df$sp)
spl <- stringr::str_split(sp_list_df$sp, "_")
spl_list <- c()
for (i in 1:length(sp_list_df$sp)) spl_list <- append(spl_list, paste0(spl[i][[1]][1], '_', spl[i][[1]][2]))
sp_list_df$sci_name <- spl_list

# Retrieve taxonomic ranks
rank_df <- tibble(tax_id=double(0), superclass=character(0), class=character(0), subclass=character(0), infraclass=character(0), order=character(0), family=character(0))
for (tax_id in sp_list_df$tax_id) {
	taxize_df <- classification(tax_id, db = 'ncbi')[[1]]
	superclass <- filter(taxize_df, rank == 'superclass')$name; if(length(superclass)==0) {superclass <- NA}
	class <- filter(taxize_df, rank == 'class')$name; if(length(class)==0) {class <- NA}
	subclass <- filter(taxize_df, rank == 'subclass')$name; if(length(subclass)==0) {subclass <- NA}
	infraclass <- filter(taxize_df, rank == 'infraclass')$name; if(length(infraclass)==0) {infraclass <- NA}
	order <- filter(taxize_df, rank == 'order')$name; if(length(order)==0) {order <- NA}
	family <- filter(taxize_df, rank == 'family')$name; if(length(family)==0) {family <- NA}
	df <- tibble(tax_id=tax_id, superclass=superclass, class=class, subclass=subclass, infraclass=infraclass, order=order, family=family)
	rank_df <- bind_rows(rank_df, df)
	Sys.sleep(5)
}
sp_list_df <- sp_list_df %>%
	left_join(rank_df, by='tax_id') %>%
	select(sci_name, assembly_id, organism_id, tax_id, superclass, class, subclass,infraclass, order, family)

# Match the fishtree
sp_list_df <- sp_list_df %>%
    mutate(tree = sci_name %in% tree$tip.label)

# Save
sp_list_df %>%
	write_tsv(args$output)