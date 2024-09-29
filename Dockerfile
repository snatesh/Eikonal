# syntax=docker/dockerfile:1
FROM ubuntu:24.04

# update aptitude and install dependencies available therein
RUN apt update && apt install -y python3 python3-numpy python3-matplotlib gcc g++ build-essential make cmake libgtest-dev libopenblas-dev liblapacke-dev git openssh-client python3-pip libomp-dev libopenblas-openmp-dev swig
# copy local version of Eikonal repo
RUN mkdir Eikonal
COPY LICENSE README.md CMakeLists.txt Eikonal/
COPY testing Eikonal/testing/
COPY include Eikonal/include/
COPY src Eikonal/src/
COPY wrapper Eikonal/wrapper/
COPY python Eikonal/python/
# clone and build nlopt 
RUN git clone https://github.com/stevengj/nlopt.git
RUN cd nlopt && mkdir build && cd build && cmake -DPython_EXECUTABLE=/usr/bin/python3.12 .. && make && make install
# export env paths so linkers pick up relevant libraries
RUN export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib
RUN export PYTHONPATH=${PYTHONPATH}:/usr/local/lib/python3.12/site-packages
# sym link for shared lib to be picked up by g++ linker
RUN ln -s /usr/local/lib/libnlopt.so.0 /usr/lib/libnlopt.so.0
# try building eikonal libs and exec
RUN cd /Eikonal && mkdir build && cd build && cmake .. && make && make install
RUN export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/Eikonal/lib 
RUN ctest --verbose
