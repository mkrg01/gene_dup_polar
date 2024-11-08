#!/bin/bash

# odb11v0_all_fasta.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_all_fasta.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_all_fasta.tab.gz
gunzip data/raw_data/orthodb/odb11v0_all_fasta.tab.gz

# odb11v0_all_og_fasta.tab.gz
# wget -nv https://data.orthodb.org/download/odb11v0_all_og_fasta.tab.gz -P data/raw_data/orthodb/
# md5sum data/raw_data/orthodb/odb11v0_all_og_fasta.tab.gz
# gunzip data/raw_data/orthodb/odb11v0_all_og_fasta.tab.gz

# odb11v0_levels.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_levels.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_levels.tab.gz
gunzip data/raw_data/orthodb/odb11v0_levels.tab.gz

# odb11v0_species.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_species.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_species.tab.gz
gunzip data/raw_data/orthodb/odb11v0_species.tab.gz

# odb11v0_level2species.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_level2species.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_level2species.tab.gz
gunzip data/raw_data/orthodb/odb11v0_level2species.tab.gz

# odb11v0_genes.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_genes.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_genes.tab.gz
gunzip data/raw_data/orthodb/odb11v0_genes.tab.gz

# odb11v0_gene_xrefs.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_gene_xrefs.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_gene_xrefs.tab.gz
gunzip data/raw_data/orthodb/odb11v0_gene_xrefs.tab.gz

# odb11v0_OGs.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_OGs.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_OGs.tab.gz
gunzip data/raw_data/orthodb/odb11v0_OGs.tab.gz

# odb11v0_OG2genes.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_OG2genes.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_OG2genes.tab.gz
gunzip data/raw_data/orthodb/odb11v0_OG2genes.tab.gz

# odb11v0_OG_xrefs.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_OG_xrefs.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_OG_xrefs.tab.gz
gunzip data/raw_data/orthodb/odb11v0_OG_xrefs.tab.gz

# odb11v0_OG_pairs.tab.gz
wget -nv https://data.orthodb.org/download/odb11v0_OG_pairs.tab.gz -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/odb11v0_OG_pairs.tab.gz
gunzip data/raw_data/orthodb/odb11v0_OG_pairs.tab.gz

# README.txt
wget -nv https://data.orthodb.org/download/README.txt -P data/raw_data/orthodb/
md5sum data/raw_data/orthodb/README.txt