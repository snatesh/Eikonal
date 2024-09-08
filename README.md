# Eikonal
A modal spectral method for the Eikonal equation on simplicial tesselations

## Docker Containerization ##
It is likely most simple to install and use the libarary from within a `Docker` container. 
For user convenience, a `Dockerfile` is provided which can generate a `Docker` image with
all required dependencies and the `Eikonal` library installed. Note, the `Docker Engine` must
be installed on your system (See https://docs.docker.com/engine/), along with `git`.

Once those dependencies are met, executing (on `Linux` or `Unix` systems) the following commands will 
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
bash shell within the container, within which you can read/write/execute
files as you please. The last two commands above change directory into
the `Eikonal` folder within the container (`/Eikonal`), compile the 
`gtest` code, and run the tests therein.

As another example, to generate Gaussian-like quadrature from order `n=4` Koornwinder 
polynomials that can integrate order `m=6` polynomials exactly (`n` < `m`), 
try (from the `Eikonal` directory)

```shell
cd bin
./ngj_quad_opt
```

## Dependencies ##
- `g++` compiler (tested on V13.2.0) 
- `build-essentials` (particularly GNU `make` utility)
- google testing framework (available via dpkg and apt as `libgtest-dev`)
- `dh-autoreconf`,`autoreconf`,`autotools` for `SNOPT` installation (currently not used)
-  `NLOPT` (open source non linear optimization librarh, `https://github.com/stevengj/nlopt.git`)
- `cblas` and `lapack` - Installation is easiest via package manager as:
```shell
sudo apt install libopenblas-openmp-dev liblapacke-dev
```
   which ensures the `C` wrapper to lapack (in headers `lapacke.h`) is 
   installed in a sane location.

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

which (on `Linux` systems) will by default install header and include files in
`/usr/local/include` and `/usr/local/lib`.

## Docker Containerization ##
It is likely most simple to install and use the libarary from within a `Docker` container. 
For user convenience, a `Dockerfile` is provided which can generate a `Docker` image with
all required dependencies and the `Eikonal` library installed. Note, the `Docker Engine` must
be installed on your system (See https://docs.docker.com/engine/). 

Once that dependency is met, executing (on `Linux` or `Unix` systems) the following commands will 
build the Docker base image for the project and test the installation within a running 
instance of the image (Docker container).
```shell


