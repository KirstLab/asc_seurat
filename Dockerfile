# Ubuntu latest
FROM kirstlab/asc_seurat:dynverse

# Owner
MAINTAINER Felipe Marques de Almeida <marques.felipe@aluno.unb.br>
SHELL ["/bin/bash", "-c"]

# Set workdir
WORKDIR /app
COPY www /app/www

# Install CRAN packages
# biocmanager
RUN R -e 'install.packages("BiocManager", dep = T)'
RUN R -e 'library(BiocManager)'
# tidyverse
RUN apt-get install -y r-base-dev xml2 libxml2-dev libssl-dev libcurl4-openssl-dev unixodbc-dev
RUN R -e 'install.packages("tidyverse", dep = T)'
RUN R -e 'library(tidyverse)'
# seurat
RUN R -e 'install.packages("Seurat", dep = T)'
RUN R -e 'library(Seurat)'
RUN R -e 'install.packages("SeuratObject", dep = T)'
RUN R -e 'library(SeuratObject)'
# patchwork
RUN R -e 'install.packages("patchwork", dep = T)'
RUN R -e 'library(patchwork)'
# ggplot2
RUN R -e 'install.packages("ggplot2", dep = T)'
RUN R -e 'library(ggplot2)'
# svglite
RUN apt-get install -y r-cran-svglite
# circlize
RUN R -e 'install.packages("circlize", dep = T)'
RUN R -e 'library(circlize)'
# reactable
RUN R -e 'install.packages("reactable", dep = T)'
RUN R -e 'library(reactable)'
# sctransform
RUN R -e 'install.packages("sctransform", dep = T)'
RUN R -e 'library(sctransform)'
# shiny
RUN R -e 'install.packages("shiny", dep = T)'
RUN R -e 'library(shiny)'
RUN R -e 'install.packages("shinyWidgets", dep = T)'
RUN R -e 'library(shinyWidgets)'
RUN R -e 'install.packages("shinyFeedback", dep = T)'
RUN R -e 'library(shinyFeedback)'
RUN R -e 'install.packages("shinycssloaders", dep = T)'
RUN R -e 'library(shinycssloaders)'
# rclipboard
RUN R -e 'install.packages("rclipboard", dep = T)'
RUN R -e 'library(rclipboard)'
# future
RUN R -e 'install.packages("future", dep = T)'
RUN R -e 'library(future)'
# ggthemes
RUN R -e 'install.packages("ggthemes", dep = T)'
RUN R -e 'library(ggthemes)'
# metap (must install multtest first)
RUN R -e 'BiocManager::install("multtest")'
RUN R -e 'library(multtest)'
RUN R -e 'install.packages("metap", dep = T)'
RUN R -e 'library(metap)'
# DT
RUN R -e 'install.packages("DT", dep = T)'
RUN R -e 'library(DT)'
# dplyr
RUN R -e 'install.packages("dplyr", dep = T)'
RUN R -e 'library(dplyr)'
# hdf5r
RUN apt-get install -y libhdf5-dev
RUN R -e 'install.packages("hdf5r", dep = T)'
RUN R -e 'library(hdf5r)'

# Install Bioconductor packages
# complex heatmap
RUN apt-get install -y  r-cran-cluster r-bioc-complexheatmap
RUN R -e 'BiocManager::install("ComplexHeatmap")'
RUN R -e 'library(ComplexHeatmap)'
# tradeseq
RUN R -e 'BiocManager::install("tradeSeq")'
RUN R -e 'library(tradeSeq)'
# single cell experiment
RUN R -e 'BiocManager::install("SingleCellExperiment")'
RUN R -e 'library(SingleCellExperiment)'
# slingshot
RUN R -e 'BiocManager::install("slingshot")'
RUN R -e 'library(slingshot)'
# biomart
RUN R -e 'BiocManager::install("biomaRt")'
RUN R -e 'library(biomaRt)'
# topgo
RUN R -e 'BiocManager::install("topGO")'
RUN R -e 'library(topGO)'

# Install Docker
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add - && \
	add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu focal stable" && \
	apt update && \
	apt install -y docker-ce

# Configuring Docker
RUN usermod -aG docker root

# Get server files
COPY global.R /app/global.R
COPY server.R /app/server.R
COPY ui.R /app/ui.R

# expose port
EXPOSE 3838

# Fix permissions
RUN chmod a+rwx -R /app/*

# Init image
COPY /scripts/bscripts/init_app.sh /app/init_app.sh
CMD ./init_app.sh
