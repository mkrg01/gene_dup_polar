library(tidyverse)
library(argparse)

# Parse
parser <- ArgumentParser()
parser$add_argument('--input', type='character')
parser$add_argument('--output', type='character')
args <- parser$parse_args()

# Input
pmcmc_df <- read_tsv(args$input)

pmcmc_df %>%
    arrange(pMCMC_polar) %>%
    filter(pMCMC_polar_corrected < 0.1) %>%
	write_tsv(args$output)