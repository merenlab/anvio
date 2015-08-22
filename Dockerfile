FROM ubuntu:trusty
 
ENV ANVIO_VERSION 1.0.2
 
RUN apt-get update
RUN apt-get install -y wget \
    cython \
    gsl-bin libgsl0-dbg libgsl0-dev libgsl0ldbl \
    hmmer \
    python-numpy \
    python-pip \
    python-scipy \
    sqlite3 \
    zlib1g-dev
 
RUN mkdir -p /tmp/build && \
    cd /tmp/build && \
    wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.tar.gz && \
    tar -zxvf v2.6.2.tar.gz && \
    cd Prodigal-2.6.2/ && \
    make && \
    cp prodigal /usr/bin/ && \
    cd - && \
    rm -rf /tmp/build
 
RUN pip install bottle==0.12.8 \
    hcluster==0.2.0 \
    ete2==2.3.6 \
    scikit-learn==0.16.1 \
    django==1.8.4 \
    pysam==0.8.3 \
    cython==0.23
 
RUN pip install anvio==$ANVIO_VERSION
 
COPY tests /anvio-tests
 
RUN echo "export PS1=\"\[\e[0m\e[47m\e[1;30m\] :: anvi'o :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>\[\e[m\] \[\e[0m\]\"" >>> /etc/profile.d/prompt
.sh
 
CMD /bin/bash -l
