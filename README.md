[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![readthedocs](https://readthedocs.org/projects/asc-seurat/badge/?version=latest)](https://asc-seurat.readthedocs.io/en/latest/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4623183.svg)](https://doi.org/10.5281/zenodo.4623183)

<p align="center">
  <!-- <a href="https://github.com/othneildrew/Best-README-Template">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a> -->

  <h2 align="center">Asc-Seurat</h2>

  <p align="center">
    <h3 align="center"> Analytical single-cell Seurat-based web application</h3>
    <br />
    <a href="https://asc-seurat.readthedocs.io/en/latest/index.html"><strong>See the documentation »</strong></a>
    <br />
    <br />
    <a href="https://github.com/KirstLab/asc_seurat/issues">Report Bug</a>
    ·
    <a href="https://github.com/KirstLab/asc_seurat/issues">Request Feature</a>
  </p>
</p>




<!-- ABOUT THE PROJECT -->
## About Asc-Seurat


Asc-Seurat (Analytical single-cell Seurat-based web application) is a web application based on [Shiny](https://shiny.rstudio.com/). Pronounced as “ask Seurat”, it provides a click-based, easy-to-install, and easy-to-use interface that allows the execution of all steps necessary for scRNA-seq analysis. It integrates many of the capabilities of the [Seurat](https://satijalab.org/seurat/) and [Dynverse](https://dynverse.org/) and also allows an instantaneous functional annotation of genes of interest using [BioMart](http://www.biomart.org/).

<p align="center">
<img src="https://github.com/KirstLab/asc_seurat/raw/main/docs/images/asc_seurat_workflow.png" width="700">
</p>

<p align="center">
<strong>Asc-Seurat workflow overview.</strong>
</p>

<!-- GETTING STARTED -->

## Prerequisites

To install Asc-Seurat, you need to have Docker installed on your machine. Docker needs to be properly installed and configured in the user's machine. Check the installation instructions provided by Docker at https://docs.docker.com/engine/install.

## Release notes

| :warning: WARNING          |
|:---------------------------|
| Please be aware that Asc-Seurat uses multiple R packages and that many of those are in continuous development. While the docker version of Asc-Seurat is stable, it may become outdated as the packages on wich it relies on are updated. [Here](https://asc-seurat.readthedocs.io/en/latest/packages_version.html) you can find a list of the packages used by Asc-Seurat and their versions.|


* v2.2 - **Released on February 8th, 2022**

    - Add the capacity to load a clustered dataset in the tab for the individual sample analysis.
  	- Add the capacity to load a clustered dataset in the tab for the integrated sample analysis.
  	- Genes identified as mitochondrial genes via the regex expression are now shown to the users.
  	- Changes the color scheme of the dynverse's plots to match the color scheme used by Seurat's plots.
  	- Small changes in the interface to improve usability.

  	- Fix a bug in the download of markers identified for multiple clusters in the integrated dataset. If a gene was identified as a marker in multiple clusters, a number was appended in the gene's name.
  	- Fix a bug that caused the app to crash when searching for conserved markers in an integrated dataset, and the gene was not expressed in one or more of the samples.
  	- Fix a bug where plots were exported with a dark background.
  	- Fix a bug in the advanced plots that caused expressed genes not to be identified. When using integrated datasets, the function now looks for the RNA assay instead of the integrated assay.
  	- Fix a bug where the app would crash when downloading the plots generated in the trajectory inference tab.

* v2.1 - **Released on May 26th, 2021**.

    - Changes the assay used for differential expression analysis and visualization to "RNA" when using SCTransform normalization. Therefore, "SCT" assay is used for the steps until clustering the data.
    - Changes the output of the differential expression analysis to the format required for the visualization tools.

* v2.0 - **Released on May 19th, 2021**.

    - Inclusion of SCTransform normalization
    - Addition of stacked violin plots
    - Addition of multiple-genes dot plot
    - Improvements on the user interface
    - Improvements in the app stability
    - Fix of minor bugs.

* v1.0 - **Released on March 19th, 2021**.

    - Release of Asc-Seurat.

## Installation

Download Asc-Seurat's docker image.
   ```sh
    # Download the docker image:
    docker pull kirstlab/asc_seurat
   ```

<!-- USAGE EXAMPLES -->

## Quickstart

```sh
# Starts Asc-Seurat on MacOS or Linux
docker run -v $(pwd):/app/user_work -v /var/run/docker.sock:/var/run/docker.sock -d --name Asc_Seurat --rm -p 3838:3838 kirstlab/asc_seurat

# Starts Asc-Seurat using Windows CMD
docker run -v %cd%:/app/user_work -v /var/run/docker.sock:/var/run/docker.sock -d --name Asc_Seurat --rm -p 3838:3838 kirstlab/asc_seurat

# Starts Asc-Seurat using Windows Powershell
docker run -v ${PWD}:/app/user_work -v /var/run/docker.sock:/var/run/docker.sock -d --name Asc_Seurat --rm -p 3838:3838 kirstlab/asc_seurat
```

:heavy_check_mark: Then, open your preferred web browser and paste the address http://localhost:3838/

## Usage

Please, visit https://asc-seurat.readthedocs.io/en/latest/ for a complete documentation.

## Reference
Pereira WJ, Almeida FM, Balmant KM, Rodriguez DC, Triozzi PM, Schmidt HW, Dervinis C, Pappas Jr. GJ, Kirst M. [Asc‑Seurat: analytical single‑cell Seurat‑based web application](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04472-2). BMC Bioinformatics 22, 556 (2021).

Pereira WJ, Almeida FM, Balmant KM, Rodriguez DC, Triozzi PM, Schmidt HW, Dervinis C, Pappas Jr. GJ, Kirst M. [Asc-Seurat – Analytical single-cell Seurat-based web application](https://www.biorxiv.org/content/10.1101/2021.03.19.436196v1). biorxiv, 2021.

<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.
