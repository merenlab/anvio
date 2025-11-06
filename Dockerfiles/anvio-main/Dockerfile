# To start a build in an absolutely pristine Docker environment you can use the following
# to clean up your docker (although please note that it will remove all existing containers
# and cached states):
#
#     docker rmi --force $(docker images -a --filter=dangling=true -q)
#     docker rm --force $(docker ps --filter=status=exited --filter=status=created -q)
#     docker system prune --force -a
#
# after that, you can start the build with the following:
#
#     BUILDKIT_PROGRESS=plain docker build -t meren/anvio:$ANVIO_VERSION .
#

FROM continuumio/miniconda3:4.11.0
ENV ANVIO_VERSION "8"

SHELL ["/bin/bash", "--login", "-c"]

RUN conda config --env --add channels bioconda
RUN conda config --env --add channels conda-forge

# Create a conda environment for anvi'o, activate it, and make sure it will
# always be activated
RUN conda create -n anvioenv python=3.10
RUN conda init bash
RUN conda activate anvioenv
RUN echo "conda activate anvioenv" >> ~/.bashrc

# Activate environment variables
ENV PATH /opt/conda/envs/anvioenv/bin:$PATH
ENV CONDA_DEFAULT_ENV anvioenv
ENV CONDA_PREFIX /opt/conda/envs/anvioenv

# fun stuff
RUN conda install -y conda-build
RUN conda install -y conda-verify

# Install nano early on so the container has
# something to edit files if debug is necessary
RUN conda install -y nano

# Install mamba
RUN conda install -y -c conda-forge mamba

# Setup the environment
RUN mamba install -y -c conda-forge -c bioconda python=3.10 \
        sqlite prodigal idba mcl muscle=3.8.1551 famsa hmmer diamond \
        blast megahit spades bowtie2 bwa graphviz "samtools>=1.9" \
        trimal iqtree trnascan-se fasttree vmatch r-base r-tidyverse \
        r-optparse r-stringi r-magrittr bioconductor-qvalue meme ghostscript

# try this, too. it may also fail to install. which is OK:
RUN mamba install -y -c bioconda fastani

RUN wget -qO- "https://cmake.org/files/v3.23/cmake-3.23.1-linux-"$(uname -m)".tar.gz" | tar --strip-components=1 -xz -C /usr/local

#GCC compiler
RUN apt-get update && \
    apt-get -y install gcc mono-mcs && \
    rm -rf /var/lib/apt/lists/*

# Install anvi'o from pip
RUN curl -L https://github.com/merenlab/anvio/releases/download/v8/anvio-8.tar.gz \
        --output anvio-8.tar.gz
RUN pip install anvio-8.tar.gz

# Install METABAT and DAS_TOOL 
#RUN mamba install metabat2 das_tool

# Install CONCOCT
RUN apt-get update && apt-get install -qq build-essential libgsl0-dev bedtools mummer samtools perl libssl-dev
RUN conda install cython
RUN pip install https://github.com/BinPro/CONCOCT/archive/1.1.0.tar.gz

# Install MAXBIN2 (installing fraggenescan will require cpanm, and we also need IDBA-UD)
RUN conda install -c bioconda perl-app-cpanminus
RUN cpanm --self-upgrade --sudo
RUN conda install -c bioconda idba
RUN cd /opt && wget https://downloads.sourceforge.net/project/fraggenescan/FragGeneScan1.31.tar.gz && tar zxf FragGeneScan1.31.tar.gz && cd FragGeneScan1.31 && make clean && make
RUN cpanm install LWP::Simple
ENV PERL5LIB /opt/conda/envs/anvioenv/lib/5.26.2/:/opt/conda/envs/anvioenv/lib/site_perl/5.26.2/:$PERL5LIB
RUN cd /opt && wget https://downloads.sourceforge.net/project/maxbin2/MaxBin-2.2.7.tar.gz && tar zxf MaxBin-2.2.7.tar.gz && cd MaxBin-2.2.7/src && make
RUN echo 'export PATH=/opt/FragGeneScan1.31:$PATH' >> ~/.bashrc
RUN echo 'export PATH=/opt/MaxBin-2.2.7:$PATH' >> ~/.bashrc

# Install some helper tools
RUN apt-get install vim util-linux -yy

# Cutify the environment
RUN echo "export PS1=\"\[\e[0m\e[47m\e[1;30m\] :: anvi'o v$ANVIO_VERSION :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>>\[\e[m\] \[\e[0m\]\"" >> /root/.bashrc

CMD /bin/bash -l

# To test the build, you can run this:
#
# docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:test-build
