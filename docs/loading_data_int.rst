.. _loading_data_int:

****************************************************
Loading the data and integration of multiple samples
****************************************************

To analyze multiple samples, select the third tab in the web application, named :guilabel:`Integration of multiple samples`.

.. note::

    The integration is based on Seurat's functions `FindIntegrationAnchors <https://www.rdocumentation.org/packages/Seurat/versions/4.0.0/topics/FindIntegrationAnchors>`_ and `IntegrateData <https://www.rdocumentation.org/packages/Seurat/versions/4.0.0/topics/IntegrateData>`_. For more information, see `Seurat's integration tutorial <https://satijalab.org/seurat/articles/integration_introduction.html>`_ and `Stuart,T. et al. (2019) <https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8>`_.

For the integration of multiple samples, the process is a little different. You still need to add your datasets in the :code:`data/` directory, creating a subdirectory for each sample. However, you also need to provide a configuration file containing the parameter values for each sample. During the installation, an example file named *configuration_file_for_integration_analysis.csv* will be created in your directory and can be used as a model to create your own file.

.. tip::

	The integration of samples can be biased if the parameters are not chosen appropriately. Therefore, it is recommended to explore each sample separately in the tab :guilabel:`One sample`, defining adequate parameters to remove poor quality cells before the integration.

Your configuration file must have 6 columns and a header (the column names are not restricted). They are used to specify what cells should be kept for each sample while loading the data before the integration.

Also, the columns need to be in a specific order as listed below.

 #. **Subdirectory name**: The name of the subdirectories containing your datasets. Each sample must have a unique name for its subdirectory, even if they are replicates.
 #. **Sample name (any name you prefer)**: Your choice of name for each sample. If you have replicates and want them to be considered as one in the plots and analysis, use the same name for all replicates.
 #. **Min. number of cells expressing a gene**: Include genes only if they are detected in at least this many cells.
 #. **Min. number of genes a cell must express to be included**: Include cells only if they expressed at least this number of genes.
 #. **Maximum number of genes a cell can express and still be included**: Remove cells that express more than this number of genes. This is useful to remove cells that you suspect are doublets.
 #. **Maximum percentage of genes belonging to the mitochondrial genome**: Here, the regular expression (`Regex <https://en.wikipedia.org/wiki/Regular_expression>`_) is a sequence of characters that is used to identify the genes belonging to the mitochondrial genome. For example, when using the human genome, this sequence should be "^MT-".
