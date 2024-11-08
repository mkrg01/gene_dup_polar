library(tidyverse)
library(argparse)
library(ape)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--gene', type='character')
parser$add_argument('--selected_OG', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
gene_df <- read_tsv(args$gene)
selected_OG_df <- read_tsv(args$selected_OG)

if (!dir.exists(args$output)) {
    dir.create(args$output)
}

# Create a heatmap (row: species, column: OG)
sp_OG <- gene_df %>%
	dplyr::select(OG_id, organism_id) %>%
	group_by(OG_id, organism_id) %>%
	summarise(num=n()) %>%
	ungroup() %>%
	filter(!is.na(OG_id)) %>%
	pivot_wider(names_from = OG_id, values_from = num, values_fill = 0) %>%
    dplyr::select(organism_id, all_of(selected_OG_df$OG_id))

sp_list_df <- sp_list_df %>%
    mutate(polar_fish = ifelse(polar == TRUE, 'True', 'False'))

trait_og_df <- sp_list_df %>%
    dplyr::select(sci_name, organism_id, lat_abs_max, polar_fish) %>%
	inner_join(sp_OG, by = 'organism_id')

og_id_list <- selected_OG_df$OG_id
for (og_id in og_id_list) {
    og_name <- selected_OG_df %>%
        filter(OG_id == og_id) %>%
        .$OG_name
    p <- ggplot(trait_og_df, aes(x=lat_abs_max, y=!!sym(og_id), color = polar_fish)) +
        geom_point(size=0.5) +
        scale_color_manual(values = c(True = "cyan3", False = "black")) +
        xlab('Absolute latitude (max)') +
        ylab('Copy number') +
        ggtitle(og_id) +
        theme(legend.position = "none")
    
    ggsave(paste0(args$output, '/', og_id, '.pdf'), plot=p, width=2.5, height=2.5)
}
