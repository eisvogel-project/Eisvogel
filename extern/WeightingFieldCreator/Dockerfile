# docker build -t meep .
FROM ubuntu:mantic

VOLUME /home/eisvogel

RUN apt-get update
RUN apt-get install -y pkg-config cmake
RUN apt-get install -y libhdf5-mpi-dev libharminv-dev libfftw3-dev libgsl-dev libgslcblas0
RUN apt-get install -y python3-matplotlib python3-mpi4py python3-meep-mpi-default libmeep-mpi-default-dev