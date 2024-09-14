# Eikonal
A modal spectral method for the Eikonal equation applied to continuous shortest path problems on simplicial tesselations. 

The primarily goal in developing this library is for use in solving continuous shortest path problems, with application in navigation, and more interestingly (to me) infrastructure design. For example, consider finding the optimal path a tunnel boring machine (for underground transportation infrastructure) should take in going from point A to point B, given obstructing inclusions provided by geological/civil survey data (feasibility/difficulty/cost varies depending on soil/rock composition, existing civil works, etc.), or undocumented inclusions encountered during excavation. The library can admit fast computations for such purposes, possibly in real-time/on-line, even if we assume the machine has freedom to maneuver 3-dimensionally. These are the sort of problems on which I'd like to apply the routines contained herein, which are modeled by appropriate configuration of the right-hand-side in Eikonal equations. The twist is that we aim to minimize the solution of the equation (time) by changing the boundary with pre-defined/newly-found constraints. I will think of some toy problems in 2D, and see if I can get things to work with some amusing movies to illustrate!

## About
The implementation offers a backend library for computing near-optimal Gaussian-like quadratures on the triangle (see `ngjquad_opt.cpp` for driver code), built on top of `NLOPT`, `BLAS` and `LAPACK`, my own work, theoretical results due to Koornwinder regarding a family of orthogonal polynomials on the triangle, and numerical results due to Vioreanu and Rokhlin on constructing near-optimal Gaussian-like quadratures over convex regions (VR-quadrature). 

I term the VR-quadrature approach "Discretize-Orthogonalize-Optimize" (DOO), while my approach is more "Orthogonalize-Discretize-Optimize" (ODO). This nomenclature draws analogy from similar jargon in variational/infinite-dimensional optimzation - Discretize-then-Optimize (DtO) or Optimize-then-Discretize (OtD) - as both concern a choice of an analytical or numerical first step.

The generatable quadrature rules are used to represent to high accuracy and with high efficiency functions appearing in Eikonal models under the Koornwinder polynomial basis family. 

Recent results by Townsend, Oliver and Vasil enable the construction of sparse differential operators on the standard triangle, interoperating within the Koornwinder polynomial family parameter space, and leading to low-storage, banded discrete differential operators. The library provides an implementation of these discrete PDOs. They can be composed together to form solvers for a variety of linear PDEs of arbitrary order, even with variable coefficients (coming soon (i hope)). Non-linearity is handled by passing an appropriately defined problem to `nlopt` in terms of these operators acting on minimization variables (solution coefficients in a Koornwinder expansion).

At the time of writing, the following capabilities exist:
- representing and manipulating functions in a modal sense under Koornwinder polynomial expansions (global approximation)
- differentiation up to any order via promotion and ladder operators acting only on modes of a function
- evaluation at a point in real space of a function represented in the modal basis (coefficients). This could be accelerated with Clenshaw's algorithm, generalized to 2D.
- Implementation of analytical entries to the Jacobi matrices appearing in the 2D-recurrence relation for the orthogonal Koornwinder polynomials, generalized for any choice of parameters.
- Quadrature discovery with optimization parameters from `nlopt` exposed for such use. In particular, I developed an interlaced scheme where first a constrained non-linear optimization problem is solved, then, if needed, a Newton relaxation, followed by another non-linear problem with a different objective, and a last on-demand Newton relaxation.

   - The gist is:
      With `z=(x,y,w)` in `R^(3N)`, minimize over `z` `norm(P(x,y).w - e1)^2`, where `e1` is the first basis vector in `R^M` with `N < M`, `P(x,y)` is the `NxM` Vandermonde matrix, `(x,y)` is constrained to lie within the triangle, and `sum(w)` is constrained to 1, followed by a Newton relaxation to find roots of `P(z).w - e1 = 0`. Then the argmin is passed to another constrained non-linear problem, wherein we try to minimize `cond(P(x,y))` (in the `2-norm`) with `w` fixed, and constrain `norm(P(x,y).w-e1)` to be at most the objective minimum from the first optimization problem. This seems to find Gaussian-like quadrature rules at the limit proposed by Vioreanu and Rokhlin (in terms of how much `N < M`), but with the added benefit of `O(1)` condition number for the interpolation matrix. This latter feature is important for modal discretizations of PDEs, as when we impose boundary conditions *n*odally, the interpolation operator on those nodes must be well conditioned if we hope for the assembled systems to be invertible.

For example, the image below depicts the output when attempting to generate a 45 (n = 9) node quadrature rule which can exactly integerate 105 polynomials (m = 14). In this case, we do not employ any Newton relaxation.
 
![alt text](https://github.com/snatesh/Eikonal/blob/main/testing/testdata/output.png?raw=true)

The red nodes are our new found abscissa, while the black ones are the initial points. We see the new nodes are well conditioned, and we correctly integrate the test function
sin(x^2+y^2) over the standard right triangle with vertices (0,0), (0,1), (1,0) (plug it into Wolfram alpha for reference). 

There remains much to be done for this project, on both paper and metal! 

## TODO Ideas for Quadrature Generation
We seeks roots to the function `F(x,y,w) = P(x,y).w - I = 0`, where `P` is the transpose of the vandermonde matrix `P(x,y)^T`, and `(x,y)` are a set of points in the standard right triangle. There is some theoretical limit on the solvability of the equation in terms of the number of points (#columns of `P`) and the number of integrals for which we want exact quadrature (#rows of `P` or length of `I`). We know it is not solvable at the optimal Gaussian quadrature limit, as the Jacobi matrices for Koornwinder polynomials on the simplex are not a commuting family. The limit is explored numerically for more general convex regions by Rokhlin and Vioreanu, and an empirical cutoff formula is provided in their work. Assuming we are under this limit in terms of the system size, and further assuming a root exists (ignoring conditioning of `P`), then `||F||^2` in the `2-norm` should also admit a zero at some point `(x,y,w)`. Rather than trying to minimize `||F||^2` using non-linear optimization techniques like SQP, or any method where quadratic approximations to sub-problems are used, we might be as well successful building a quadratic approximation of the objective `G=||F||^2` and using the Newton direction `-(\nabla^2 G(x,y,w))^{-1}\nabla G(x,y,w)` to find the root. The only "approximation" to be made here in the model is for the gradient and hessian entries concerning variable `w`. We can hopefully use simple centered finite difference schemes there. Variables `x,y` enter only in `P`, and we should know by now how to take derivatives of polynomials. Ill-conditioning in the argmin can be handled with our current approach (`nlopt` call on `min cond(P)`). 

   - the first step is seeing if incorporating analytical gradient information even helps convergence in the `nlopt` routines currently used. After all, if it isn't broken, don't fix it, and we are already able to generate high order quadrature with a large increase in efficiency relative to the past (days->hours on the same OTS equipment). I will need to check if the ones currently in use even need derivatives (many are derivative-free algorithms), and whether there are others that require derivative info. 
   - if `nlopt` responds positively to analytical gradient info, then we can just plug and play the Hessian into my current Newton method implementation, discarding the backtracking line search part. Let's see what happens!

Newton's method has fast quadratic convergence, but there is difficulty of step length control for gradient descent which leads to constraint violations for the quadrature rule. If we can get Newton descent right, it may open up the possibility for on-line `p`-adaptivity in PDE discretizations using the Koornwinder bases (which also reminds me of the possiblity of `(a,b,c)` adaptivity - interoperating between basis families parameterized by `(a,b,c)`).

## TODO Ideas for generalizing to d>2
The only challenge with all of this has been in deriving matrix entries for the 3-term matrix recurrence relations defining spaces of polynomials orthogonal to all spaces of lower order (for `d > 1`). We need these entries to assemble Jacobi matrices, and we need those to initialize any descent method for finding Gaussian-like quadrature/cubature. By definition, `Jn_i.P = x_i.P + [0\\A_{n-1}_i \bb{P}_n]`. In `2D`, if `{x_1,x_2}_j` is a set of nodes admitting a well conditioned interpolation matrix on the triangle, then we should be able to solve to mach_eps for the entries in `Jn_1,Jn_2`. All we would need to do so (and, in fact, already have) are the entries to `A_{n-1}_{1,2}`. In higher dimensions, this reduces the courageous pencil/paper work by quite a bit, even if we use symbolic algebra systems. Then, we can consider JEVD techniques used in SP problems like BSS to generalize the initialization method for `d>2` (complexification doesn't work for 3D, which is used in 2D for the eigenvalue matching problem. Maybe we can use quaternions and projection to dim under for 3D? Is there even an eigenvalue solver that supports quaternions??).

### Computing approximate joint eigenvalue decompositions
Consider orthogonal polynomials in dimension $d$ defined in a convex region of $\mathbb{R}^d$. Suppose that the Jacobi matrices $J_k \in \mathbb{R}^{\text{dim}\Pi_{n-1}^d$$, $k=0,\cdots d-1$ corresponding to a 3-term recurrence relation on subspaces of homogoneous polynomials with total degree $m=0,\cdots n-1$.  

## Docker Containerization
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
the `/Eikonal` folder in the container context, compile the 
`gtest` code (located in `/Eikonal/testing/gtest.cpp`), and run the tests therein.

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

-  The c++ compoiler should have support for the `OpenMP` shared memory parallelization library, and
   the library must exist on your system. That is, `omp.h` and `libomp.so` must exist somewhere in the filesystem,
   the compiler must understand `OpenMP` directives, and the linker should be able to find and link to `libomp`
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


