.. _loading_data:

****************************************
Loading the data of an individual sample
****************************************

Location of your dataset
========================

For Asc-Seurat to read your sample dataset(s), they need to be located in a subdirectory inside the :code:`data/` directory. The :code:`data/` directory will be created during the installation and contains a subdirectory with an example dataset called :code:`example_PBMC/`. This dataset is from the publicly available `10×'s Peripheral Blood Mononuclear Cells (PBMC) <https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>`_ and contains 2700 cells.

.. figure:: images/data_tree.png
   :width: 50%
   :align: center

   Organization of the :code:`data/` directory.

Therefore, to add your dataset, create a subdirectory inside :code:`data/` containing the counts' matrix (*matrix.mtx.gz*), cell barcodes (*barcodes.tsv.gz*), and gene names (*features.tsv.gz*).

Asc-Seurat provides separated environments (tabs) for the analysis of a single sample and the integrated analysis of multiple samples. Below are the instructions to load your data depending on if you are using one or multiple samples.

Loading the data
================

To analyze an individual sample, select the second tab in the web application, named :any:`One sample`. Then, select the sample to analyze and set the initial criteria to exclude cells that should not be load, as shown below.

 After inserting your datasets in the :code:`data/` directory, the samples will be available to load in Asc-Seurat, as shown below.

.. figure:: images/loading_one_sample.png
   :width: 90%
   :align: center

   Example of how to load an individual sample for analysis and of the requested initial parameters.

In the first box to the left, it is possible to select the sample to use. However, there are a few parameters that you need to provide before loading your data. This step is based on Seurat's functions `CreateSeuratObject <https://www.rdocumentation.org/packages/Seurat/versions/3.1.4/topics/CreateSeuratObject>`_ and `PercentageFeatureSet <https://satijalab.org/seurat/reference/PercentageFeatureSet.html>`_. Between parenthesis, we list the name of the parameter in the CreateSeuratObject function.

Below is a description of these parameters:

 * **Project name**: Sets the name for the project. This will appear in some of the plots but it is not required (project).
 * **Min. number of cells expressing a gene**: Include genes only if they are detected in at least this many cells (min.cells).
 * **Min. number of genes a cell must express to be included**: Include cells only if they expressed at least this number of genes (min.features).
 * **Regex to identify mitochondrial genes**: Here, the regular expression (`Regex <https://en.wikipedia.org/wiki/Regular_expression>`_) is a sequence of characters that is used to identify the genes belonging to the mitochondrial genome (pattern). For example, when using the human genome, this sequence should be "^MT-".

After setting the parameters described above, click on the button :guilabel:`Load data of the selected sample` to start the analysis. A violin plot showing the distribution of cells will appear. This plot can then be use to set more restrictive parameters for :ref:`quality control <quality_control>`.
