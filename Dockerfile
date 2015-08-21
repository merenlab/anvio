FROM ubuntu:trusty

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

RUN mkdir -p /tmp/build && cd /tmp/build && wget https://github.com/hyattpd/Prodigal/archive/v2.6.2.tar.gz && tar -zxvf v2.6.2.tar.gz && cd Prodigal-2.6.2/ && make && cp prodigal /usr/bin/ && cd -

RUN pip install bottle==0.12.8 \
    hcluster==0.2.0 \
    ete2==2.3.6 \
    scikit-learn==0.16.1 \
    django==1.8.4 \
    pysam==0.8.3 \
    cython==0.23
RUN pip install anvio==1.0.1

COPY tests /tmp/build/tests
COPY anvio/data /tmp/build/anvio/data

# anvi-import-collection: error: unrecognized arguments: -a test-output/ANNOTATION.db
#RUN cd /tmp/build/tests && ./run_mini_test.sh

RUN rm -rf /tmp/build
