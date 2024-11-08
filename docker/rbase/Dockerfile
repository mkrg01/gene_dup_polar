FROM rocker/tidyverse:4.3.3

# install by apt-get
RUN apt-get update && apt-get install -y \
        locales \
 && rm -rf /var/lib/apt/lists/*

# Generate C.UTF-8 locale
RUN locale-gen C.UTF-8

# Set environment variables
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV LANGUAGE C.UTF-8

# Install from GitHub
RUN installGithub.r \
        ropensci/taxize

# Install from CRAN
RUN install2.r --deps TRUE --repos 'https://cloud.r-project.org' \
        BiocManager \
        argparse \
        phytools \
        ape \
        fishtree

# Install by BioManager
RUN Rscript -e "BiocManager::install(c('treeio', 'ggtree', 'ggtreeExtra', 'Biostrings', 'ggmsa'))"
