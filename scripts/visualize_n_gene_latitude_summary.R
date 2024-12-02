library(tidyverse)
library(argparse)

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

# Create a heatmap (row: species, column: OG)
sp_OG <- gene_df %>%
	dplyr::select(OG_id, organism_id) %>%
    filter(OG_id %in% selected_OG_df$OG_id) %>%
	group_by(OG_id, organism_id) %>%
	summarise(num=n()) %>%
	ungroup() %>%
    pivot_wider(names_from = OG_id, values_from = num, values_fill = 0) %>%
    pivot_longer(!organism_id, names_to='OG_id', values_to='num')

sp_list_df <- sp_list_df %>%
    mutate(Polar_fish = ifelse(polar == TRUE, 'True', 'False'))

trait_og_df <- sp_OG %>%
    inner_join(sp_list_df, by = 'organism_id') %>%
    dplyr::select(organism_id, OG_id, num, lat_abs_max, Polar_fish)

trait_og_df$OG_id <- factor(trait_og_df$OG_id, levels=selected_OG_df$OG_id)
trait_og_df$Polar_fish <- factor(trait_og_df$Polar_fish, levels=c('True', 'False'))

p <- ggplot(trait_og_df) +
    geom_point(aes(x=lat_abs_max, y=num, color = Polar_fish), size=0.5) +
    facet_wrap(~OG_id, ncol=4, scales='free') +
    scale_color_manual(values = c(True = '#00AAFF', False = 'black')) +
    xlab('Max absolute latitude') +
    ylab('Copy number')

ggsave(args$output, plot=p, width=9, height=10)