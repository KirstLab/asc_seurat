# A Guide for Running asc_seurat Container on a Cluster

Please follow these steps:

1. Copy the asc_seurat.job file to your home directory by running this command:
   ```
   cp /storage/hpc/group/ircf/software/singularity_asc_seurat/asc_seurat.job ~/asc_seurat.job
   ```

2. Go to your home directory by executing:
   ```
   cd
   ```

3. Check for an available node using this command:
   ```
   sinfo --state=idle
   ```

4. Adjust the "asc_seurat.job" file according to your requirements for Partition, Node, Memory, etc.:
    ```
    #!/bin/bash
    ##SBATCH -p r630-hpc3
    ##SBATCH -w lewis4-r630-hpc3-node548
    #SBATCH -p Gpu
    #SBATCH -t 0-02:00  # time (days-hours:minutes)
    #SBATCH --ntasks-per-node=10
    #SBATCH --mem=100G
    #SBATCH --output=/home/%u/log_asc_seurat.job.%j
    ##SBATCH --mail-user=youremail@missouri.edu  # email address for notifications
    ##SBATCH --mail-type=END,FAIL  # which type of notifications to send
    #SBATCH -J asc_seurat
    ```

5. Submit the job to the SLURM scheduler by running:
   ```
   sbatch asc_seurat.job
   ```

6. Check the job log and follow the instructions provided in it:
   ```
   cat log_asc_seurat.job.*
   ```