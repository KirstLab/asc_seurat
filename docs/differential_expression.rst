.. _differental_expression:

***********************************************************
Markers identification and differential expression analysis
***********************************************************

After clustering the cells, users may be interested in identifying genes specifically expressed in one cluster (markers) or in genes that are differentially expressed among clusters of interest. Asc-Seurat can apply multiple algorithms to identify gene markers for individual clusters or to identify differentially expressed genes (DEGs) among clusters, using Seurat's functions `FindMarkers <https://satijalab.org/seurat/reference/FindMarkers.html>`_ and `FindAllMarkers <https://satijalab.org/seurat/reference/FindConservedMarkers.html>`_.

Asc-Seurat allows users to filter gene markers and DEGs by the fold change and minimal percentage of cells expressing a gene in the cluster(s). Moreover, users can define a significance level to exclude genes based on the adjusted p-value (see below).

.. figure:: images/DE_one_sample_1.png
   :width: 65%
   :align: center

   Example of Asc-Seurat's interface showing the settings to the search for gene markers for each of the clusters using the Wilcox test.

.. figure:: images/DE_one_sample_2.png
   :width: 65%
   :align: center

   Example of Asc-Seurat's interface showing the settings to the search for markers for a specific cluster (cluster 0).

.. figure:: images/DE_one_sample_3.png
   :width: 65%
   :align: center

   Example of Asc-Seurat's interface showing the settings to search for DEGs genes among clusters 0, 2, and 3.

An iterative table will be available after executing the search for marker or DEGs, showing the significant genes. Moreover, users can download the list of significant markers or DEGs as a csv file.

.. figure:: images/DEG_table.png
   :width: 100%
   :align: center

   List of the ten most significant markers identified for cluster 3 of the PBMC dataset (as defined in :ref:`clustering`).

The list of genes in the csv can then be used to visualize their gene expression in a series of plots, as shown in the section :ref:`expression_visualization`.
