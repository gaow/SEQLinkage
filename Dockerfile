FROM debian:jessie

RUN apt-get update && \
	apt-get -y install bzip2 wget && \
  apt-get clean

RUN wget --quiet http://bioinformatics.org/spower/download/.private/SEQDocker.tar.bz2
RUN mkdir -p /opt/miniconda3/envs/snake
RUN tar jxf SEQDocker.tar.bz2 -C /opt/miniconda3/envs/snake
RUN rm -f SEQDocker.tar.bz2

RUN useradd -ms /bin/bash seqx
USER seqx
RUN mkdir -p /home/seqx/data
WORKDIR /home/seqx/data

RUN wget --quiet http://bioinformatics.org/spower/download/.private/SEQLinkageResources.tar.bz2
RUN tar jxf SEQLinkageResources.tar.bz2 -C /home/seqx
RUN rm -f SEQLinkageResources.tar.bz2

ENV PATH /home/seqx/.SEQLinkage/bin:/opt/miniconda3/envs/snake/bin:$PATH
ENV PYTHONPATH /opt/miniconda3/envs/snake/lib/python2.7/site-packages:$PYTHONPATH

RUN spower -h > /dev/null
RUN seqlink -h > /dev/null