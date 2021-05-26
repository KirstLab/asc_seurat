[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![readthedocs](https://readthedocs.org/projects/asc-seurat/badge/?version=latest)](https://asc-seurat.readthedocs.io/en/latest/) [![DOI](https://zenodo.org/badge/340176419.svg)](https://zenodo.org/badge/latestdoi/340176419)

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

## Realese notes

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

:heavy_check_mark: Then, open your preferred web browser and paste the address https://localhost:3838/

## Usage

Please, visit https://asc-seurat.readthedocs.io/en/latest/ for a complete documentation.

## Reference
[1] Pereira WJ, Almeida FM, Balmant KM, Rodriguez DC, Triozzi PM, Schmidt HW, Dervinis C, Pappas Jr. GJ, Kirst M. [Asc-Seurat – Analytical single-cell Seurat-based web application](https://www.biorxiv.org/content/10.1101/2021.03.19.436196v1). biorxiv, 2021.

<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.
