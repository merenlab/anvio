FROM ubuntu:trusty
 
ENV ANVIO_VERSION 1.2.3

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        gsl-bin \
        hmmer \
        libgsl0-dbg \
        libgsl0-dev \
        libgsl0ldbl \
        python-dev \
        python-numpy \
        python-pip \
        python-scipy \
        sqlite3 \
        wget \
        zlib1g-dev \
        libhdf5-serial-dev \
        libhdf5-dev \
    && mkdir -p /tmp/build \
    && wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.tar.gz -O /tmp/build/v2.6.2.tar.gz \
    && tar -zxvf //tmp/build/v2.6.2.tar.gz -C /tmp/build/ \
    && make -C /tmp/build/Prodigal-2.6.2/ \
    && cp /tmp/build/Prodigal-2.6.2/prodigal /usr/bin/ \
    && pip install \
        cython==0.23 \
        bottle==0.12.8 \
        hcluster==0.2.0 \
        ete2==2.3.6 \
        scikit-learn==0.16.1 \
        django==1.8.4 \
        pysam==0.8.3 \
        h5py==2.5.0 \
    && pip install anvio==$ANVIO_VERSION \
    && apt-get remove -y \
        binutils \
        build-essential \
        g++ \
        g++-4.8 \
        gcc \
        gcc-4.8 \
        make \
        patch \
    && apt-get -y autoremove \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

COPY tests /anvio-tests

RUN echo "export PS1=\"\[\e[0m\e[47m\e[1;30m\] :: anvi'o :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>>\[\e[m\] \[\e[0m\]\"" >> /etc/profile.d/prompt.sh

CMD /bin/bash -l
