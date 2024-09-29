#!/bin/bash

# number of threads for OpenMP
num_threads=6

# let the shell use the maximum amount of stack memory
ulimit -s unlimited
# set mem for thread stack
export OMP_STACKSIZE=256m

#thread pinning settings
export OMP_PLACES="{0}:${num_threads}:1"
export OMP_PROC_BIND=true
export OMP_DISPLAY_ENV=true

# LD_LIBRARY_PATH for so python interpreter can find libs
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$PWD/lib

#export OMP_DISPLAY_AFFINITY=true
#export OMP_NUM_THREADS=${num_threads}

# eventually want to specify if on 32-bit or 64 bit
# on ubuntu, we can use
# uname -a | grep 64
# bit64=$(echo $?)
# if 1, we are on 64. else 32
# this is needed to define memory alignment masks for 
# aligned allocators.
