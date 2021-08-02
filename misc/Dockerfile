FROM ubuntu:20.04
COPY ./EM /EM 
RUN apt-get update && apt-get -y install curl make gcc gfortran gnupg2 software-properties-common
RUN curl http://xmi-apt.tomschoonjans.eu/xmi.packages.key | apt-key add -
#RUN add-apt-repository "deb [arch=amd64] http://xmi-apt.tomschoonjans.eu/ubuntu focal stable"
RUN add-apt-repository -s "deb [arch=amd64] http://xmi-apt.tomschoonjans.eu/ubuntu focal stable"
RUN apt-get update && apt-get -y install libxrl11 libxrl11-dev libxrlf03-11 libxrl-perl python3-libxrl11 libnetcdff-dev
RUN apt-get update && apt-get -y install cmake apt-utils fftw3-dev gdb
RUN cd /EM/src/ && mkdir build && cd build && cmake .. && make
RUN echo 'alias xrfpm="/EM/src/build/bin/xrfpm"' >> ~/.bashrc
RUN echo 'alias medgen="/EM/src/build/bin/medgen"' >> ~/.bashrc
