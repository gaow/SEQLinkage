FROM ubuntu:trusty

MAINTAINER Gao Wang (wang.gao@columbia.edu)

# Both RV-NPL and RarePedSim run python 2.7
# Install miniconda2 and dependencies: https://hub.docker.com/r/conda/miniconda2/dockerfile
RUN apt-get -qq update && apt-get -qq -y install software-properties-common \
    && add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get -qq update && apt-get -qq -y install wget build-essential gfortran bzip2 swig gcc-5 g++-5 libboost1.55-all-dev \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 1 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 1

RUN wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /opt/miniconda2 \
    && rm -rf /tmp/miniconda.sh

ENV PATH=/opt/miniconda2/bin:${PATH}
RUN conda update conda \
    && conda install -y matplotlib scipy \
    && conda clean --all --yes

RUN conda install -c bioconda tabix \
    && conda clean --all --yes

# Install cstatgen
RUN cd /tmp && wget https://github.com/statgenetics/cstatgen/archive/master.tar.gz && tar xzvf master.tar.gz && cd cstatgen-master && python setup.py install && rm -rf /tmp/*
# Install SEQLinkage
ARG DUMMY=unknown
RUN cd /tmp && wget https://github.com/gaow/SEQLinkage/archive/master.tar.gz && tar xzvf master.tar.gz && cd SEQLinkage-master && python setup.py install && rm -rf /tmp/*
RUN mkdir /.SEQLinkage/ && cd /.SEQLinkage/ && wget http://bioinformatics.org/spower/download/.private/genemap.txt && mkdir ./bin && cd ./bin && wget http://bioinformatics.org/spower/download/.private/linux/mlink http://bioinformatics.org/spower/download/.private/linux/unknown http://bioinformatics.org/spower/download/.private/linux/makeped http://bioinformatics.org/spower/download/.private/linux/pedcheck
ENV PYTHON_EGG_CACHE=/tmp/Python-Eggs

# Default command
CMD ["bash"]

# To build
# docker build --build-arg DUMMY=`date +%s` -t gaow/seqlink .
# docker push gaow/seqlink