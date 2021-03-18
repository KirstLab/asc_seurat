.. _differental_expression:

***********************************************************
Markers identification and differential expression analysis
***********************************************************

After clustering the cells, users might be interest on identifying genes that are specifically expressed on one cluster (markers) or on genes that are differentially expressed among clusters of interest. Asc-Seurat allows the application of multiple algorithms to identify gene markers for individual clusters or to identify differentially expressed genes (DEGs) among clusters. For that, it deploys Seurat's functions `FindMarkers <https://satijalab.org/seurat/reference/FindMarkers.html>`_ and `FindAllMarkers <https://satijalab.org/seurat/reference/FindConservedMarkers.html>`_.

On Asc-Seurat, allows users to filter gene markers and DEGs by the Fold change and minimal percentage of cells expressing a gene in the cluster(s). Moreover, users can define a significance level to exclude genes based on the adjusted p-value (see below).

.. figure:: images/DE_one_sample_1.png
   :width: 80%
   :align: center

   Example of Asc-Seurat's interface showing the settings to the search for gene markers for each of the clusters using the Wilcox test.

.. figure:: images/DE_one_sample_2.png
   :width: 80%
   :align: center

   Example of Asc-Seurat's interface showing the settings to the search for markers for a specific clusters (cluster 0).

.. figure:: images/DE_one_sample_3.png
   :width: 80%
   :align: center

   Example of Asc-Seurat's interface showing the settings to the search for DEGs genes among clusters 0, 2 and 3.

After executing the search for marker or DEGs, an iterative table will be available showing the significant genes. Moreover, user can download the list of significant markers or DEGs as a csv file.

.. figure:: images/DEG_table.png
   :width: 80%
   :align: center

The list of genes contained in the csv can then be used to visualize their gene expression in a series of plots, as shown in the section :ref:`expression_visualization`.
