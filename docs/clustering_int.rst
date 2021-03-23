.. _clustering_int:

**********
Clustering
**********

After filtering your integrated data to remove low quality cells, Asc-Seurat allows the clustering of the remaining cells according to their expression profiles. However, before the clustering, some steps are necessary.

A series of steps are executed before the clustering steps, including normalization, scaling and dimensional reduction via PCA. Moreover, users need to decide how many dimensions are to be used during the clustering. Asc-Seurat provides an elbow plot to inform this decision. A description of these steps, and of the required parameters, is shown below.

Normalization and detection of the most variable genes
======================================================

.. note::

    So far, Asc-Seurat only allows the `LogNormalize <https://satijalab.org/seurat/reference/LogNormalize.html>`_ normalization. `SCTransform <https://satijalab.org/seurat/reference/SCTransform.html>`_ will be added soon.


Asc-Seurat allows the normalization using Seurat's `LogNormalize <https://satijalab.org/seurat/reference/LogNormalize.html>`_ function. Users have the option to change the scaling factor if necessary but it is typically not needed. In the same window (see the image below), users can select what method should be used to identify the most variable genes, and how many of the most variable genes should be used during the dimension reduction (PCA).

The most variable genes are genes that exhibit high cell-to-cell variation in the dataset and therefore are more informative. We use Seurat's function `FindVariableFeatures <https://satijalab.org/seurat/reference/FindVariableFeatures.html>`_ The default setting should work well for the majority of cases.

.. figure:: images/normalization_settings.png
   :alt: Quality control.
   :width: 100%
   :align: center


Dimensional reduction (PCA)
===========================

As mentioned above, only the most variable genes are use in the PCA (Principal Components Analysis). The PCA will be executed using Seurat's function `RunPCA <https://satijalab.org/seurat/reference/RunPCA.html>`_ and, after its conclusion, an `elbow plot <https://satijalab.org/seurat/reference/ElbowPlot.html>`_ is generated automatically, to help users to decide how many PCs should be include to inform the clustering step.

As shown below, users can easily download the elbow plot. Also, users should set the number of PCs to include during clustering in the windows at the right side of the plot. For the PBMC integrated dataset, 20 PCs will be used during the clustering.

.. figure:: images/PCA_int.png
   :alt: Quality control.
   :width: 100%
   :align: center

Clustering of cells
====================

The next step is the clustering of the cells. For that, Asc-Seurat used both `FindNeighbors <https://satijalab.org/seurat/reference/FindNeighbors.html>`_ and `FindClusters <https://satijalab.org/seurat/reference/FindClusters.html>`_ functions of Seurat package.

Before the execution, however, users need to set a value for the resolution parameter. The resolution is an important parameter to evaluate because it determines the profile and number of clusters identified for a dataset. Selecting larger values will favor splitting cells into more clusters while selecting a smaller value has the opposite effect. Quoting from `Seurat's tutorial: <https://satijalab.org/seurat/archive/v1.4/pbmc3k_tutorial.html>`_ "We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets".

.. tip::

	There is no easy way to define an optimal value for the resolution parameter. Users need to try different values and evaluate the resulting clusters according with the expectation for their cells population. Visualizing the expression profile of cell-type specific markers can provide a hint if the chosen value is too small or too large.

After the execution of the clustering step, three plots are generated for cluster visualization, all of them using the Uniform Manifold Approximation and Projection (UMAP) technique. The first plot shows the clustering of the whole dataset colored by cluster. The second plot shows the same plot but cells are colored by sample. The third plot shows the clustering of the cells of each sample, with one subplot per sample.

.. figure:: images/clustering_int.png
   :alt: Quality control.
   :width: 100%
   :align: center

   Clustering of the PBMC integrated dataset using 20 PCs and a resolution value of 0.5.

Selecting clusters of interest
------------------------------

In some cases, it is interesting to select or exclude some clusters of cells from your dataset before executing the subsequent steps. This is useful, for example, when users desire explore a developmental trajectory of a specific group of cell types.

Asc-Seurat makes this step simple. Users only need to select the cluster(s) to keep or to exclude and start the reanalysis of the remaining cells by clicking on :guilabel:`Reanalyze after selection/exclusion of clusters`, see below.

.. figure:: images/excluding_cells_p1.png
   :alt: Quality control.
   :width: 100%
   :align: center

   Asc-Seurat make it easy to select or exclude a cluster (or clusters) of cells. In this example, we exclude all cells belonging to the cluster 0.

Asc-Seurat will then execute the steps with the new set of cells up to the PCA. Then, users need to evaluate the elbow plot and decide the number of PCs to use in the clustering of the new set of cells. Also, users can either keep the same value for the resolution parameter or modify it before clicking on :guilabel:`Rn the clustering analysis` to start the clustering once more.

.. figure:: images/clustering_int_2.png
  :alt: Quality control.
  :width: 100%
  :align: center

  Clustering of the PBMC integrated dataset after excluding cells belonging to the cluster 0 from the original dataset.


.. warning::

	The cluster's numbering will change every time that you select or exclude cluster(s).
