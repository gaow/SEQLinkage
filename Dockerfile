FROM debian:jessie

RUN apt-get update && \
	apt-get -y install bzip2 wget

RUN echo "Downloading SEQDocker ..."
RUN wget --quiet http://bioinformatics.org/spower/download/.private/SEQDocker.tar.bz2
RUN wget --quiet http://bioinformatics.org/spower/download/.private/SEQLinkageResources.tar.bz2

RUN echo "Install SEQDocker ..."
RUN mkdir -p /opt/gaow
RUN tar jxf SEQDocker.tar.bz2 -C /opt/gaow
RUN tar jxf SEQLinkageResources.tar.bz2 -C $HOME
ENV PATH $HOME/.SEQLinkage/bin:/opt/gaow/bin:$PATH
ENV PYTHONPATH /opt/gaow/lib/python2.7/site-packages:$PYTHONPATH
RUN echo "Installation complete!"