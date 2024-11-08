FROM rocker/tidyverse:4.4.1

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

# Install from CRAN
RUN install2.r --deps TRUE --repos 'https://cloud.r-project.org' \
        argparse \
        ape \
        phytools \
        MCMCglmm \
	foreach \
	doParallel
