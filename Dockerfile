# syntax=docker/dockerfile:1
FROM ubuntu:24.04

# update aptitude and install dependencies available therein
RUN apt update && apt install -y python3 python3-numpy python3-matplotlib gcc g++ build-essential make cmake libgtest-dev libopenblas-dev liblapacke-dev git openssh-client python3-pip libopenblas-openmp-dev
# copy local version of Eikonal repo
RUN mkdir Eikonal
COPY LICENSE README.md Makefile Eikonal/
COPY testing Eikonal/testing/
COPY include Eikonal/include/
COPY src Eikonal/src/
# clone and build nlopt 
RUN git clone https://github.com/stevengj/nlopt.git
RUN cd nlopt && mkdir build && cd build && cmake .. && make && make install
# sym link for shared lib to be picked up by g++ linker
RUN ln -s /usr/local/lib/libnlopt.so.0 /usr/lib/libnlopt.so.0
# try building eikonal libs and exec
RUN cd /Eikonal && make 
