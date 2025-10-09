# creates an image containing an installed instance of seqwish

# load base image
FROM ubuntu:22.04

# update system
RUN apt-get update > /dev/null

# install dependencies
RUN apt-get -qqy install zlib1g zlib1g-dev libomp-dev

# install build tools
RUN apt-get -qqy install build-essential cmake git

# copy over current directory to container
ADD . /seqwish

# set working directory
WORKDIR /seqwish

# build
RUN cmake -H. -Bbuild && cmake --build build -- -j3

# cleanup
RUN apt-get -qy autoremove

# add seqwish to the PATH
ENV PATH /seqwish/bin:$PATH
