docker login
# docker exec -it gtc_dashboard "sh"
# docker system prune -a
docker build -t winuthayanon/php_db_ldap:alpine . -f Dockerfile_Alpine
# docker build -t winuthayanon/php_db_ldap:ubuntu . -f Dockerfile_Ubuntu
docker push winuthayanon/php_db_ldap:alpine
# docker push winuthayanon/php_db_ldap:ubuntu
singularity pull docker://winuthayanon/php_db_ldap:alpine
singularity run -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html php_db_ldap_alpine.sif
singularity inspect php_db_ldap_alpine.sif
# singularity instance start -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html php_db_ldap_alpine.sif php_db_ldap_alpine
# singularity instance list
# singularity instance stop php_db_ldap_alpine

sudo singularity build php_db_ldap.sif php_db_ldap.def
singularity inspect php_db_ldap.sif
singularity run -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html php_db_ldap.sif
# nohup singularity run -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html php_db_ldap.sif > /dev/null 2>&1 &
singularity instance start -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html php_db_ldap.sif php_db_ldap
singularity instance list
singularity instance stop php_db_ldap
singularity instance start php_db_ldap
singularity shell instance://php_db_ldap


singularity build --fakeroot php_db_ldap.sif php_db_ldap.def
singularity build php_db_ldap.sif php_db_ldap.def
singularity run --bind $PWD/data/html:/var/www/html php_db_ldap.sif
sudo singularity build php_db_ldap.sif php_db_ldap.def



singularity run --bind $PWD/data/html:/var/www/html php_db_ldap.sif
sudo singularity build php_db_ldap.sif php_db_ldap.def
singularity run php_db_ldap.sif
singularity run --writable-tmpfs -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs php_db_ldap.sif
singularity run --writable-tmpfs -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html php_db_ldap.sif
singularity run --writable-tmpfs -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html -B tmp:/tmp php_db_ldap.sif
singularity run -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html -B tmp/sessions:/var/www/sessions php_db_ldap.sif

docker exec -it gtc_dashboard "sh"

docker build -t winuthayanon/gtc_asp2php:alpine . -f Dockerfile_Alpine
docker push winuthayanon/gtc_asp2php:alpine
singularity pull docker://winuthayanon/gtc_asp2php:alpine
sudo singularity build gtc_asp2php_alpine.sif gtc_asp2php_alpine.def
singularity inspect gtc_asp2php_alpine.sif
singularity inspect -d gtc_asp2php_alpine.sif
singularity inspect -r gtc_asp2php_alpine.sif
mkdir -p run var/lib/nginx/logs
singularity run -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html gtc_asp2php_alpine.sif
# nohup singularity run -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html gtc_asp2php_alpine.sif > /dev/null 2>&1 &
singularity instance start -B run:/run -B var/lib/nginx/logs:/var/lib/nginx/logs -B data/html:/var/www/html gtc_asp2php_alpine.sif gtc_asp2php
singularity instance list
singularity instance stop gtc_asp2php
singularity instance start gtc_asp2php
singularity shell instance://gtc_asp2php

# https://developer.nvidia.com/cuda-12-0-0-download-archive?target_os=Linux&target_arch=x86_64&Distribution=Ubuntu&target_version=22.04&target_type=deb_network
docker system prune -a
docker build -t winuthayanon/gpu:12.0.0-devel-ubuntu22.04 . -f Dockerfile_12.0.0-devel-ubuntu22.04
docker push winuthayanon/gpu:12.0.0-devel-ubuntu22.04

docker build -t winuthayanon/gpu:12.1.1-runtime-ubuntu22.04 . -f Dockerfile_12.1.1-runtime-ubuntu22.04
docker push winuthayanon/gpu:12.1.1-runtime-ubuntu22.04

docker build -t winuthayanon/gpu:12.1.1-devel-ubuntu22.04 . -f Dockerfile_12.1.1-devel-ubuntu22.04
docker push winuthayanon/gpu:12.1.1-devel-ubuntu22.04

docker pull winuthayanon/gpu:ubuntu
singularity pull docker://winuthayanon/gpu:ubuntu
singularity run --nv gpu_ubuntu.sif
singularity pull docker://winuthayanon/gpu:12.1.1-runtime-ubuntu22.04
singularity run --nv gpu_ubuntu.sif
singularity pull docker://winuthayanon/gpu:12.1.1-devel-ubuntu22.04
singularity run --nv gpu_ubuntu.sif

module load cuda
module list
module unload cuda
module list
module avail
module load cuda
module load singularity

docker pull nvcr.io/nvidia/tensorflow:23.05-tf2-py3
singularity pull docker://nvcr.io/nvidia/tensorflow:23.05-tf2-py3

# https://github.com/dizcza/docker-hashcat
singularity cache clean
sinfo -p Gpu -o %n,%G
srun --pty -p Gpu --time=0-02:00 --mem=100G --ntasks-per-node 10 /bin/bash
srun -p Gpu -N1 -n20 -t 0-02:00 --mem=100G --gres gpu:1 --pty /bin/bash
srun -p Gpu -N1 -n20 -t 0-04:00 --mem=100G --gres gpu:1 --pty /bin/bash
srun -p Gpu --account ircf -N1 -n20 -t 0-2:00 --mem=100G --gres gpu:1 --pty /bin/bash

singularity pull docker://dizcza/docker-hashcat:latest
singularity pull docker://dizcza/docker-hashcat:cuda

export VERSION=1.20.5 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc

sudo apt install build-essential libssl-dev uuid-dev libgpgme-dev squashfs-tools libseccomp-dev wget pkg-config git cryptsetup-bin libglib2.0-dev

./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

# https://guiesbibtic.upf.edu/recerca/hpc/running-singularity-containers-with-gpu
singularity run --nv /soft/singularity/tensorflow_20.08-tf2-py3.sif python -c "import tensorflow as tf; print('Num GPUs Available: ',len(tf.config.experimental.list_physical_devices('GPU'))); print('Tensorflow version: ',tf.__version__)"

# https://docs.sylabs.io/guides/latest/user-guide/gpu.html
singularity pull docker://tensorflow/tensorflow:latest-gpu
rsync -avPe "ssh -p 50022" ssc@beenplus.com:./git/gpu/*.sif .
singularity run --nv tensorflow_latest-gpu.sif
srun -p Gpu --account ircf -N1 -n20 -t 0-2:00 --mem=100G --gres gpu:1 --pty /bin/bash
srun -p Gpu --account ircf -N1 -n20 -t 0-2:00 --mem=100G --gres gpu:3 --pty /bin/bash
srun -p Gpu --account ircf -N2 -n20 -t 0-2:00 --mem=100G --gres gpu --pty /bin/bash
SINGULARITYENV_CUDA_VISIBLE_DEVICES=0 singularity run --nv tensorflow_latest-gpu.sif
SINGULARITYENV_CUDA_VISIBLE_DEVICES=0,1,2 singularity run --nv tensorflow_latest-gpu.sif
SINGULARITYENV_CUDA_VISIBLE_DEVICES=0 singularity run --nv gpu_ubuntu.sif
SINGULARITYENV_CUDA_VISIBLE_DEVICES=0,1,2 singularity run --nv gpu_ubuntu.sif
SINGULARITYENV_CUDA_VISIBLE_DEVICES=0,1,2 singularity run --nv gpu_12.1.1-runtime-ubuntu22.04.sif
SINGULARITYENV_CUDA_VISIBLE_DEVICES=0,1,2 singularity run --nv gpu_12.1.1-devel-ubuntu22.04.sif
SINGULARITYENV_CUDA_VISIBLE_DEVICES=all singularity run --nv gpu_12.1.1-devel-ubuntu22.04.sif
export CUDA_VISIBLE_DEVICES=0,1,2
singularity run --nv gpu_12.1.1-devel-ubuntu22.04.sif

Singularity> 
cd bin
ln -sf /usr/bin/hashcat python
export PATH=$HOME/bin:$PATH
python -I
python -b -D2 -m 1800

CUDA API (CUDA 12.0)
====================
* Device #1: Tesla V100-PCIE-32GB, 32190/32500 MB, 80MCU
* Device #2: Tesla V100-PCIE-32GB, 32190/32500 MB, 80MCU
* Device #3: Tesla V100-PCIE-32GB, 32190/32500 MB, 80MCU

OpenCL API (OpenCL 2.0 pocl 1.8  Linux, None+Asserts, RELOC, LLVM 11.1.0, SLEEF, DISTRO, POCL_DEBUG) - Platform #1 [The pocl project]
=====================================================================================================================================
* Device #4: pthread-Intel(R) Xeon(R) Gold 6238 CPU @ 2.10GHz, skipped

Benchmark relevant options:
===========================
* --opencl-device-types=2
* --optimized-kernel-enable

--------------------------------------------------------------------
* Hash-Mode 1800 (sha512crypt $6$, SHA512 (Unix)) [Iterations: 5000]
--------------------------------------------------------------------

Speed.#1.........:   345.9 kH/s (74.25ms) @ Accel:16384 Loops:256 Thr:32 Vec:1
Speed.#2.........:   345.9 kH/s (74.24ms) @ Accel:16384 Loops:256 Thr:32 Vec:1
Speed.#3.........:   346.8 kH/s (74.04ms) @ Accel:16384 Loops:256 Thr:32 Vec:1
Speed.#*.........:  1038.6 kH/s

# 3090
Speed.#1.........:   476.5 kH/s (69.97ms) @ Accel:8 Loops:256 Thr:1024 Vec:1

Singularity> python
Python 2.7.15+ (default, Jul  9 2019, 16:51:35)
[GCC 7.4.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> from tensorflow.python.client import device_lib
>>> print(device_lib.list_local_devices())

singularity pull docker://kirstlab/asc_seurat
singularity pull docker://winuthayanon/asc_seurat
sudo singularity build asc_seurat.sif asc_seurat.def

docker build -t winuthayanon/asc_seurat:20230530 . -f Dockerfile_20230530
docker push winuthayanon/asc_seurat:20230530
docker build -t winuthayanon/asc_seurat:20230608 . -f Dockerfile_20230608
docker push winuthayanon/asc_seurat:20230608
# sudo singularity build asc_seurat.sif asc_seurat.def
sudo singularity build --force asc_seurat.sif asc_seurat.def
singularity inspect asc_seurat.sif
singularity inspect -d asc_seurat.sif
singularity inspect -r asc_seurat.sif
singularity instance list
singularity instance stop asc_seurat

docker run -v "$(pwd):/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock \
-v /tmp:/tmp -d --name Asc_Seurat --rm -p 3838:3838 winuthayanon/asc_seurat

mkdir -p user_work var/run tmp
singularity instance start -B user_work:/app/user_work -B var/run:/var/run \
-B tmp:/tmp asc_seurat.sif asc_seurat
singularity run -B user_work:/app/user_work -B var/run:/var/run \
-B tmp:/tmp asc_seurat.sif


docker build -t winuthayanon/asc_seurat:20230615 . -f Dockerfile_20230615
docker push winuthayanon/asc_seurat:20230615
docker run -v "$(pwd):/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock \
-v /tmp:/tmp -d --name Asc_Seurat --rm -p 3838:3838 winuthayanon/asc_seurat:20230615
docker run -v "$(pwd):/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock \
-v /tmp:/tmp --name Asc_Seurat --rm -p 3838:3838 winuthayanon/asc_seurat:20230530
docker run -v "$(pwd):/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock \
-v /tmp:/tmp --name Asc_Seurat --rm -p 3838:3838 winuthayanon/asc_seurat:20230615

docker run -v "$(pwd):/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock \
-w /app/user_work \
-it winuthayanon/asc_seurat:20230615 \
R -e "shiny::runApp('/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"

docker run -v "$(pwd):/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock \
-p 4949:4949 \
-w /app/user_work \
-it winuthayanon/asc_seurat:20230615 \
R -e "shiny::runApp('/app', host = '0.0.0.0', port = 4949, launch.browser = FALSE)"

docker run -v "$(pwd)/nine:/app/user_work" -v /var/run/docker.sock:/var/run/docker.sock \
-p 4949:4949 \
-it winuthayanon/asc_seurat:20230615 \
R -e "shiny::runApp('/app', host = '0.0.0.0', port = 4949, launch.browser = FALSE)"

# sudo singularity build asc_seurat.sif asc_seurat.def
sudo singularity build --force asc_seurat.sif asc_seurat.def
singularity inspect asc_seurat.sif
singularity inspect -d asc_seurat.sif
singularity inspect -r asc_seurat.sif
singularity instance list
singularity instance stop asc_seurat
singularity run -B user_work:/app/user_work -B var/run:/var/run \
-B tmp:/tmp asc_seurat.sif

# PORT=4848
readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
echo ${PORT}
SIF_FILE=asc_seurat.sif
mkdir -p mydata
singularity exec \
  -B mydata:/app/user_work \
  -B /var/run:/var/run \
  --writable-tmpfs \
  --cleanenv ${SIF_FILE} \
  /app/init_app.sh && \
  R -e "shiny::runApp('/app', host = '0.0.0.0', port = ${PORT}, launch.browser = FALSE)"

readonly PORT=$(python -c 'import socket; s=socket.socket(); s.bind(("", 0)); print(s.getsockname()[1]); s.close()')
echo ${PORT}
SIF_FILE=asc_seurat.sif
mkdir -p mydata
singularity exec \
  -B mydata:/app/user_work \
  -B /var/run:/var/run \
  --writable-tmpfs \
  --cleanenv ${SIF_FILE} \
  /app/init_app.sh ${PORT}