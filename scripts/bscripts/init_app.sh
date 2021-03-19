#!/bin/bash

# First step -- Enter in users work directory
cd /app/user_work

# Next, check if directories exist
## DATA
if [ ! -d "data" ]; then
  # Create dir
  mkdir data data/example_PBMC

  # Download data
  wget https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
  tar zxvf pbmc3k_filtered_gene_bc_matrices.tar.gz
  mv filtered_gene_bc_matrices/hg19/* data/example_PBMC
  rm -rf filtered_gene_bc_matrices pbmc3k_filtered_gene_bc_matrices.tar.gz
fi
## RDS_files
if [ ! -d "RDS_files" ]; then
  # Create dir
  mkdir RDS_files
fi

if [ ! -d "www" ]; then
  # Create dir
  cp -R /app/www .
fi

# Open server
R -e "shiny::runApp('/app', host = '0.0.0.0', port = 3838, launch.browser = F)"
