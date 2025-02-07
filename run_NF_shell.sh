#!/bin/bash

# Rstudio server needs write access to /run and /home/rstudio
# Thus, create two writeable directories
mkdir -p $HOME/bind/run
mkdir -p $HOME/bind/rstudio
mkdir -p $HOME/bind/tmp

FILE=$HOME/vroni/rocker4csbg.sif

# execute rserver in the container
# bind the directories created above
# start a server on localhost, port 8888
singularity shell \
  --bind $HOME/bind/run:/run \
  --bind $HOME/bind/rstudio:/home/rstudio \
  --bind $HOME/bind/tmp:/tmp \
  --bind /vscratch/nfortelny/:/media/vroni \
  $FILE

#  --bind /usr/local/AGFORTELNY:/media/AGFORTELNY \
