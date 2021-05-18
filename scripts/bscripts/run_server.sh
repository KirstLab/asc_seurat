#!/bin/bash

# Simple script to start the Asc_Seurat shiny server using docker
#
# Author: Felipe M. Almeida (almeidafmarques@outlook.com)

# Help
Help()
{
	# Display Help
cat << EOF

Simple script for starting the ASC-seurat shiny app.
It must be executed in the working dir (where the data dir for analysis are stored)

 Syntax: run_server.sh [-h|s|d] [-p <number>] [-t <tag>]

 Options:

 -h/--help			Print this help
 -p/--port			Set out a custom PORT for server listening. [ Default: 3838 ]
 -s/--start			Start asc_seurat server
 -r/--remove			Removes the started asc_seurat server
 -d/--download			Only download the required docker image
 -t/--tag           Tag name pointing to Asc-Seurat version.

EOF
}

# Defaults
TAG="v2.0"
PORT="3838"
DOCKER_IMAGE=kirstlab/asc_seurat:${TAG}
DOCKER_CONTAINER_NAME=Asc_Seurat

# Download image
Download()
{
	# Run docker pull
	docker pull "$DOCKER_IMAGE"
}

# Start server
Start()
{

	# Tells user
	echo "The server has started in: http://localhost:${PORT}/"

	# Start server in current directory
	echo "When finished, run the command:"
	echo -e    "\tUsing the container name: docker rm -f $DOCKER_CONTAINER_NAME"
	echo -e -n "\tUsing the container id:   docker rm -f "
	docker run \
		-v $(pwd):/app/user_work \
		-v /var/run/docker.sock:/var/run/docker.sock \
		--detach --rm \
		--name "$DOCKER_CONTAINER_NAME" \
		--publish 3838:"$PORT" \
		"$DOCKER_IMAGE"

}

# Remove server
Remove()
{
	echo -e "Removing server $DOCKER_CONTAINER_NAME!"
	docker rm -f $DOCKER_CONTAINER_NAME
}

# No arguments given
if [ $# -eq 0 ] ; then
	Help
	exit
fi

# Get positional arguments
POSITIONAL=()
while [[ $# -gt 0 ]]
do
ARGS="$1"
case $ARGS in
    -h|--help)
      Help
      shift
      ;;
	-d|--download)
		Download
		shift
		;;
	-s|--start)
		Start
		shift
		;;
	-r|--remove)
		Remove
		shift
		;;
	-p|--port)
		if [ "$2" ]; then
	      if [ "$2" -eq "$2"  ] 2>/dev/null ; then
        PORT=$2
        shift 2
      else
        echo -e '\nERROR: "-p/--port" requires a numeric argument'
        echo -e "argument parsed: $2 \n"
        exit 1
      fi
    else
	      echo -e '\nERROR: "-p/--port" requires a numeric argument\n'
	fi
		;;
    -t|--tag)
		#if [ "$2" ]; then
	    TAG=$2
        DOCKER_IMAGE=kirstlab/asc_seurat:${TAG}
        shift 2
    #else
    #	fi
		;;
    *)
      printf "******************************\n"
      printf "Error: Invalid argument $1\n"
      printf "******************************\n"
      exit 1
esac
done
