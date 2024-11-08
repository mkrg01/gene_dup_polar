library(tidyverse)
library(argparse)
library(ape)
library(phytools)
library(MCMCglmm)
library(foreach)
library(doParallel)

# Parse
parser <- ArgumentParser()
parser$add_argument('--sp_list', type='character')
parser$add_argument('--gene_list', type='character')
parser$add_argument('--target_OG', type='character')
parser$add_argument('--tree', type='character')
parser$add_argument('--output', type='character')
parser$add_argument('--threads', type='integer')
args <- parser$parse_args()

# Input
sp_list_df <- read_tsv(args$sp_list)
gene_df <- read_tsv(args$gene_list)
tree <- read.tree(args$tree)
target_OG_df <- read_tsv(args$target_OG)
threads <- args$threads
target_OG_list <- target_OG_df$OG_id

og_df <- gene_df %>%
        dplyr::select(OG_id, organism_id) %>%
        filter(OG_id %in% target_OG_list) %>%
        group_by(OG_id, organism_id) %>%
        summarise(num=n()) %>%
        ungroup() %>%
        filter(!is.na(OG_id)) %>%
        pivot_wider(names_from = OG_id, values_from = num, values_fill = 0)

trait_og_df <- sp_list_df %>%
    dplyr::select(sci_name, organism_id, polar, lat_abs_max) %>%
    left_join(og_df, by='organism_id') %>%
    dplyr::select(-organism_id)

if (!is.ultrametric(tree)) {
    tree <- force.ultrametric(tree)
}
Ainv <- inverseA(tree)$Ainv
prior <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002)))
data <- data.frame(trait_og_df)


cl <- makeCluster(threads)
registerDoParallel(cl)

res_df <- foreach(i = 1:length(target_OG_list), .combine = dplyr::bind_rows, .packages = c('MCMCglmm', 'tibble', 'ggplot2')) %dopar% {
    target_OG <- target_OG_list[i]

    # Polar
    formula <- as.formula(paste0('X', target_OG, '~polar'))

    # Test run to identify appropriate "nitt", "burnin", and "thin" parameters
    m1 <- MCMCglmm(
        fixed = formula,
        random = ~sci_name,
        family = "poisson",
        prior = prior,
        data = data,
        ginverse = list(sci_name = Ainv),
        nitt = 10000,
        burnin = 0,
        thin = 1
    )
    unit_df <- tibble(step = 1:length(m1$VCV[,2]), unit = as.numeric(m1$VCV[,2]))
    p <- ggplot(data=unit_df) +
        geom_line(aes(x=step, y=unit))
    ggsave(paste0('data/genome/mcmcglmm/mcmcglmm_testrun_polar_', target_OG, '.pdf'), width=5, height=5, plot=p)

    # Run
    m2 <- MCMCglmm(
        fixed = formula,
        random = ~sci_name,
        family = "poisson",
        prior = prior,
        data = data,
        ginverse = list(sci_name = Ainv),
        nitt = 10005000,
        burnin = 5000,
        thin = 1000
    )

    # Latitude
    formula_latitude <- as.formula(paste0('X', target_OG, '~lat_abs_max'))

    # Test run
    m3 <- MCMCglmm(
        fixed = formula_latitude,
        random = ~sci_name,
        family = "poisson",
        prior = prior,
        data = data,
        ginverse = list(sci_name = Ainv),
        nitt = 10000,
        burnin = 0,
        thin = 1
    )
    unit_df <- tibble(step = 1:length(m3$VCV[,2]), unit = as.numeric(m3$VCV[,2]))
    p <- ggplot(data=unit_df) +
        geom_line(aes(x=step, y=unit))
    ggsave(paste0('data/genome/mcmcglmm/mcmcglmm_testrun_latitude_', target_OG, '.pdf'), width=5, height=5, plot=p)
    
    # Run
    m4 <- MCMCglmm(
        fixed = formula_latitude,
        random = ~sci_name,
        family = "poisson",
        prior = prior,
        data = data,
        ginverse = list(sci_name = Ainv),
        nitt = 10005000,
        burnin = 5000,
        thin = 1000
    )

    # Save
    save(m2, file = paste0('data/genome/mcmcglmm/mcmcglmm_result_polar_', target_OG, '.RData'))
    save(m4, file = paste0('data/genome/mcmcglmm/mcmcglmm_result_latitude_', target_OG, '.RData'))
    res_polar <- summary(m2)
    res_latitude <- summary(m4)
    tibble(OG_id = target_OG, DIC_polar = res_polar$DIC, DIC_latitude = res_latitude$DIC, pMCMC_polar = res_polar$solutions[2,5], pMCMC_latitude = res_latitude$solutions[2,5])
}

target_OG_df %>%
    left_join(res_df, by = 'OG_id') %>%
    write_tsv(args$output)

stopCluster(cl)