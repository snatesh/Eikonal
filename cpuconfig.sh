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
export OMP_NUM_THREADS=${num_threads}
