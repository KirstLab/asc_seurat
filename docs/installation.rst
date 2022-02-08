.. _installation:

************
Installation
************

Dependencies
============

Asc-Seurat relies on multiple R packages and their dependencies (See :ref:`references`). However, we provide a Docker image that contains all necessary software and packages.

To install Asc-Seurat, it is necessary to have Docker installed on the machine. Docker needs to be correctly installed and configured in the user's machine. Check the installation instructions provided by Docker at https://docs.docker.com/engine/install.

.. warning::

   Single-cell RNA-seq data analysis can be resource-consuming. By default, Docker will use (allocate) only a fraction of your RAM memory. A minimum requirement of 8 Gb of RAM memory was necessary to analyze a dataset containing around eight thousand cells during our tests. Therefore, users need to adjust the amount of allocated memory according to their dataset. Please visit: https://docs.docker.com/docker-for-mac/space/ (MAC) or https://docs.docker.com/docker-for-windows/ (Windows) to learn how to make this adjustment.

Image download
--------------

After installing Docker, users can download the Docker image containing Asc-Seurat by executing the command below. The installation is quick and straightforward. After that, everything is set.

.. code-block:: bash

    # Download the docker image:
    docker pull kirstlab/asc_seurat

Starting Asc-Seurat
===================

After downloading the image, users can start the app on their working directory. See below for the instructions on how to start the app in the different operational systems.

.. note::

  During the first execution, some folders will be created in the working directory. They include the folders :code:`data/` and :code:`RDS_files/` that users will use to store their datasets, allowing Asc-Seurat to read them.

  Always start the run inside the working directory to be able to use the data inside these folders.

For macOS and Linux
-------------------
.. _test_tip_as_target:
.. tip::

    The code below will automatically update Asc-Seurat to the latest version. You can download and execute a specific version of Asc-Seurat by adding the version's tag to the image's name, i.e., replace :code:`kirstlab/asc_seurat` by :code:`kirstlab/asc_seurat:v2.1` to use v2.1.

.. code-block:: bash

   # Create the working directory
   mkdir my_project
   cd my_project

   # Starts Asc-Seurat
   docker pull kirstlab/asc_seurat && docker run -v "$(pwd):/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock -v /tmp:/tmp -d --name Asc_Seurat --rm -p 3838:3838 kirstlab/asc_seurat

.. note::

    After executing the "docker run" command, open your preferred web browser and paste the address http://localhost:3838/. Asc-Seurat should be ready.

If users want to kill the Docker container, run the command below.

.. code-block:: bash

   docker kill Asc_seurat

For Windows
-----------

To run Asc-Seurat on Windows via Docker, it is necessary to use Windows 10. Moreover, Windows Subsystem for Linux (WSL) needs to be installed. Before running Asc-Seurat, users must guarantee that Docker and its WSL 2 components are correctly installed and running. For that, check the two (sequential) tutorials below:

1. `Docker installation info <https://docs.docker.com/docker-for-windows/install/>`_
2. `Define windows WSL 2 as default <https://docs.microsoft.com/en-us/windows/wsl/install-win10#step-5---set-wsl-2-as-your-default-version>`_ (If you followed the link above correctly, you only need to execute step 5 of this tutorial).

The tutorials above contain all the necessary information to install Docker on Windows. However, it is also possible to find video tutorials on YouTube. Check the following link for an example: https://youtu.be/5nX8U8Fz5S0 .

After certifying that everything is working, Asc-Seurat can be started using the commands below:

.. tip::

    The code below will automatically update Asc-Seurat to the latest version. You can download and execute a specific version of Asc-Seurat by adding the version tag to the image's name, i.e., replace :code:`kirstlab/asc_seurat` by :code:`kirstlab/asc_seurat:v.1.0` to use v1.0.

.. code-block:: bash

    # Create the working directory
    mkdir my_project
    cd my_project

    # If using Windows CMD
    docker pull kirstlab/asc_seurat && docker run -v "%cd%:/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock -v /tmp:/tmp -d --rm -p 3838:3838 kirstlab/asc_seurat

    # If using Windows Powershell
    docker pull kirstlab/asc_seurat && docker run -v "${PWD}:/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock -v /tmp:/tmp -d --rm -p 3838:3838 kirstlab/asc_seurat

.. note::

    After executing the "docker run" command, open your preferred web browser and paste the address http://localhost:3838/. Asc-Seurat should be ready.

If users want to kill the Docker container, run the command below.

.. code-block:: bash

   docker kill Asc_seurat
