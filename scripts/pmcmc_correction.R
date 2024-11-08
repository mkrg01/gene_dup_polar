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
    mutate(pMCMC_polar_corrected = p.adjust(pMCMC_polar, method = 'BH')) %>%
    mutate(pMCMC_latitude_corrected = p.adjust(pMCMC_latitude, method = 'BH')) %>%
	write_tsv(args$output)