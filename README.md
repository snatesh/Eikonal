# Eikonal
A modal spectral method for the Eikonal equation on simplicial tesselations. 

## About
The implementation offers a backend library for computing near-optimal Gaussian-like quadratures on the triangle (see `ngjquad_opt.cpp` for driver code), built on top of `NLOPT`, `BLAS` and `LAPACK`, my own work, theoretical results due to Koornwinder regarding a family of orthogonal polynomials on the triangle, and numerical results due to Vioreanu and Rokhlin on constructing near-optimal Gaussian-like quadratures over convex regions (VR-quadrature). 

I term the VR-quadrature approach "Discretize-Orthogonalize-Optimize" (DOO), while my approach is more "Orthogonalize-Discretize-Optimize" (ODO). This nomenclature draws analogy from similar jargon in variational/infinite-dimensional optimzation - Discretize-then-Optimize (DtO) or Optimize-then-Discretize (OtD) - as both concern a choice of an analytical or numerical first step.

The generatable quadrature rules are used to represent to high accuracy and with high efficiency functions appearing in Eikonal models under the Koornwinder polynomial basis family. 

Recent results by Townsend, Oliver and Vasil enable the construction of sparse differential operators on the standard triangle, interoperating within the Koornwinder polynomial family parameter space, and leading to low-storage, banded discrete operators. 

At the time of writing, the following capabilities exist:
- representing and manipulating functions in a modal sense under Koornwinder polynomial expansions (global approximation)
- differentiation up to any order via promotion and ladder operators acting only on modes of a function
- evaluation at a point in real space of a function represented in the modal basis (coefficients). This could be accelerated with Clenshaw's algorithm, generalized to 2D.
- Implementation of analytical entries to the Jacobi matrices appearing in the 2D-recurrence relation for the orthogonal Koornwinder polynomials, generalized for any choice of parameters.
- Quadrature discovery with optimization parameters from `nlopt` exposed for such use. In particular, I developed an interlaced scheme where first a constrained non-linear optimization problem is solved, then, if needed, a Newton relaxation, followed by another non-linear problem with a different objective, and a last on-demand Newton relaxation.

   - The gist is:
      With z=(x,y,w) in R^(3N), minimize over z norm(P(x,y).w - e1)^2, where e1 is the first basis vector in R^M with N < M, P(x,y) is the NxM Vandermonde matrix, (x,y) is constrained to lie within the triangle, and sum(w) is constrained to 1, followed by a Newton relaxation to find roots of P(z).w - e1 = 0. Then the argmin is passed to another constrained non-linear problem, wherein we try to minimize cond(P(x,y)) (in the 2-norm) with w fixed, and constrain norm(P(x,y).w-e1) to be at most the objective minimum from the first optimization problem. This seems to find Gaussian-like quadrature rules at the limit proposed by Vioreanu and Rokhlin (in terms of how much N < M), but with the added benefit of O(1) condition number for the interpolation matrix. This latter feature is important for modal discretizations of PDEs, as when we impose boundary conditions *n*odally, the interpolation operator on those nodes must be well conditioned if we hope for the assembled systems to be invertible.

For example, the image below depicts the output when attempting to generate a 45 (n = 9) node quadrature rule which can exactly integerate 105 polynomials (m = 14). In this case, we do not employ any Newton relaxation.
 
![alt text](https://github.com/snatesh/Eikonal/blob/main/testing/testdata/output.png?raw=true)

The red nodes are our new found abscissa, while the black ones are the initial points. We see the new nodes are well conditioned, and we correctly integrate the test function
sin(x^2+y^2) over the standard right triangle with vertices (0,0), (0,1), (1,0) (plug it into Wolfram alpha for reference). 

There remains much to be done for this project, on both paper and metal! 

## Docker Containerization ##
It is likely most simple to install and use the libarary from within a `Docker` container. 
For user convenience, a `Dockerfile` is provided which can generate a `Docker` image with
all required dependencies and the `Eikonal` library installed. Note, the `Docker Engine` must
be installed on your system (See https://docs.docker.com/engine/), along with `git`.

Administrative privileges may be required to run `Docker` commands. On `Linux` systems, you must add your user to the  `Docker` group by executing:

```shell
sudo groupadd docker
sudo usermod -aG docker <user> # replace <user> with your username
newgrp docker # or log out and log back in
```

Since you need `root` privilege to execute the above commands, you could instead forgo their execution, and prefix `docker` commands with `sudo`.

Once these dependencies are met, executing (on `Linux` or `Unix` systems) the following commands will 
build the Docker base image for the project and test the installation within a running 
instance of the image (a `Docker` container):

```shell
git clone git@github.com:snatesh/Eikonal.git
cd Eikonal
docker build -t ngj_tri_opt:latest .
docker run -it ngj_tri_opt:latest bash
cd Eikonal
make test
```
Above, the `docker run` command with `bash` post-fixed will instantiate a 
`bash` shell running in the container, within which you can read/write/execute
files as you please. The last two commands above change directory into
the `/Eikonal` folder within the container, compile the 
`gtest` code (located in `/Eikonal/testing/gtest.cpp`, and run the tests therein.

As another example, to generate near optimal Gaussian-like quadrature with unit
weight sum and nodes interior to the standard triangle from order `n=4` Koornwinder 
polynomials that can integrate order `m=6` polynomials exactly, 
try:

```shell
cd /Eikonal/bin && ./ngjquad_opt
```

The initial implemenation for such a task can be found in `/Eikonal/src/ngjquad_opt.cpp`.

## Dependencies ##
- `g++` compiler (tested on V13.2.0) 
- `build-essentials` (particularly GNU `make` utility)
- google testing framework (available via dpkg and apt as `libgtest-dev`)
- `dh-autoreconf`,`autoreconf`,`autotools` for `SNOPT` installation (currently not used)
-  `NLOPT` (open source non linear optimization library, `https://github.com/stevengj/nlopt.git`)
- `cblas` and `lapack` - Installation is easiest via package manager as:
```shell
sudo apt install libopenblas-openmp-dev liblapacke-dev
```
   which ensures the `C` wrapper to lapack (in headers `lapacke.h`) is 
   installed in a sane location.

-  Thje c++ compoiler should have support for the `OpenMP` shared memory parallelization library, and
   the library must exist on your system. That is, `omp.h` and `libomp.so` must exist somewhere in the filesystem,
   the compiler must understand openMP directives, and the linker should be able to find and link to `libomp`
   given the `-fopenmp` flag. Most modern compilers will ship with the header and library files,
   as well as support for `OpenMP` directives. In case the files don't make it, you can use (on `Linux` with `dpkg`)
```shell
sudo apt install libomp-dev
```
### Thread settings for OpenMP ###
An example threading config file `cpuconfig.sh` is included to show some of the environment variables that `OpenMP` exposes 
for users to set from the shell. In general, using the number of physical cores on the system improves the performance of
most algorithsm that `Eikonal` uses (over hyperthreading), while the binding of threads spawned by OpenMP to those physical core IDs is
something with which you should experiment on your system by setting the `OMP_PROC_BIND` and `OMP_PLACES` environement variables. 
Proc-binding essentially disables hyperthreading when you set the number of threads to less than or equal to the number of cores per socket.

I find that enabling threads significantly reduces the convergence time for quadrature search when `n,m` are large enough (>~10), so using and playing around with `OpenMP` settings is well worth the effort if you need high order quadrature.

```shell
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
```

Note, the `OMP_NUM_THREADS` variable is set here, though this will change. It is not advisable to set such an environmental variable if you link to other programs which also use `OpenMP`. They may have their own tested/working heuristics for setting the number of threads. and fixing it in the shell context can mess up the performance of their threaded functions when called from within the same context. The `num_threads` variable set above will eventually be passed to the program at runtime (set by calling `omp_set_num_threads`), while `OMP_NUM_THREADS` will be unset/empty.

## NLOPT Installation ##
The open-source nonlinear optimization libary which we use is `nlopt`, by our favorite `FFTW` co-creator Steven Johnson! 
For our purposes, it is used to generate near-optimal Gaussian-like quadrature on the triangle, and is to be investigated for 
use in solving non-linear PDEs/PDE-constrained optimization problems. 

To install NLOPT, simply execute

```shell
git clone https://github.com/stevengj/nlopt.git
cd nlopt
mkdir build && cd build
cmake ..
make
sudo make isntall 
```

which (on `Linux` systems) will by default install header and shared library files in
`/usr/local/include` and `/usr/local/lib`.

### Threading in `nlopt` ###
Following the discussion of threading with `OpenMP` above, it seems that `NLOPT` is threaded using the lower level `pthreads` library, but responds to the `OMP_NUM_THREADS` variable setting (i.e. if `num_threads=6`, CPU utilization will not exceed 600%). I have yet to test whether this response is in favor of performance, in terms of time or memory use. The faster we can go, the farther we can push the order of generated quadratures, so this is worth looking into at some point.


