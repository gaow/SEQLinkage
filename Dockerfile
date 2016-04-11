FROM debian:jessie

RUN apt-get update && \
	apt-get -y install bzip2 wget

RUN echo "Downloading SEQDocker ..."
RUN wget --quiet http://bioinformatics.org/spower/download/.private/SEQDocker.tar.bz2
RUN wget --quiet http://bioinformatics.org/spower/download/.private/SEQLinkageResources.tar.bz2

RUN echo "Install SEQDocker ..."
RUN mkdir -p /opt/miniconda3/envs/snake
RUN tar jxf SEQDocker.tar.bz2 -C /opt/miniconda3/envs/snake
RUN tar jxf SEQLinkageResources.tar.bz2 -C $HOME
ENV PATH $HOME/.SEQLinkage/bin:/opt/miniconda3/envs/snake/bin:$PATH
ENV PYTHONPATH /opt/miniconda3/envs/snake/lib/python2.7/site-packages:$PYTHONPATH
RUN spower -h > /dev/null
RUN seqlink -h > /dev/null
RUN echo "Installation complete!"