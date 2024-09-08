# Eikonal
A modal spectral method for the Eikonal equation on simplicial tesselations

## Dependencies ##
- g++ compiler (tested on V13.2.0) 
- build-essentials (particularly GNU `make` utility)
- google testing framework (available via dpkg and apt as `libgtest-dev`)
- `dh-autoreconf`,`autoreconf`,`autotools` for `SNOPT` installation
-  NLOPT (open source non linear optimization lib, `https://github.com/stevengj/nlopt.git`)
- cblas and lapack installation is easiest via package manager as:
```shell
sudo apt install libopenblas-openmp-dev liblapacke-dev
```

## NLOPT Installation ##
The open-source nonlinear optimization libary which we use is `nlopt`, by our favorite `FFTW` co-creator Steven Johnson! For our purposes, it is used to generate near-optimal Gaussian-like quadrature on the triangle, and is to be investigated for use in solving non-linear PDEs/PDE-constrained optimization problems. 

To install NLOPT, simply execute

```shell
git clone https://github.com/stevengj/nlopt.git
cd nlopt
mkdir build && cd build
cmake ..
make
sudo make isntall 
```

which (on Linux systems) will by default install header and include files in
`/usr/local/include` and `/usr/local/lib`.

## Docker Containerization ##
It is likely most simple to install and use the libarary from within a docker container. 
For user convenience, a Dockerfile is provided which can generate a Docker image with
all required dependencies and the Eikonal library installed. Note, the `Docker Engine` must
be installed on your system (See https://docs.docker.com/engine/). 

Once that dependency is met, executing (on Linux or Unix systems) the following commands will 
build the Docker base image for the project and test the installation within a running 
instance of the image (Docker container).
```shell


