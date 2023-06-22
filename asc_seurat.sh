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

module load singularity

readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
echo ${PORT}

cat 1>&2 <<END

#*** Credit: Sarayut (Nine) Winuthayanon, https://www.linkedin.com/in/winuthayanons/ ***#

How to Run asc_seurat on the Slurm Cluster: 

1. Execute either of the following commands in your local computer terminal:

    ssh -N -L ${PORT}:${HOSTNAME}:${PORT} ${USER}@lewis.rnet.missouri.edu

    ssh -N -L ${PORT}:${HOSTNAME}:${PORT} ${USER}@lewis42.rnet.missouri.edu

    ssh -N -L ${PORT}:${HOSTNAME}:${PORT} ${USER}@lewis4-dtn.rnet.missouri.edu

    ssh -N -L ${PORT}:${HOSTNAME}:${PORT} ${USER}@lewis4-dtn1.rnet.missouri.edu

2. In your browser, visit 

    http://localhost:${PORT}

3. Check job status on login node:

        sacct

    # To find jobs were issued from May 2023:

        sacct -S 2023-05-01

To terminate the job:
1. Close the browser
2. On the login node:

    scancel -f ${SLURM_JOB_ID}

END

SIF_FILE=asc_seurat.sif
MYDATA=/storage/hpc/data/${USER}/mydata
mkdir -p ${MYDATA}
mkdir -p var/run

singularity exec \
  -B ${MYDATA}:/app/user_work \
  -B var/run:/var/run \
  --writable-tmpfs \
  --cleanenv ${SIF_FILE} \
  /app/init_app.sh ${PORT}

printf 'asc_seurat exited' 1>&2