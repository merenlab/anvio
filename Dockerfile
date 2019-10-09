FROM continuumio/miniconda3

RUN conda update conda
RUN conda config --env --add channels conda-forge
RUN conda config --env --add channels bioconda
RUN conda create -n anvioenv python=3.6

# Activate environment
ENV PATH /opt/conda/envs/anvioenv/bin:$PATH
ENV CONDA_DEFAULT_ENV anvioenv
ENV CONDA_PREFIX /opt/conda/envs/anvioenv

RUN echo "conda activate anvioenv" >> ~/.bashrc
SHELL ["/bin/bash", "-c"]

RUN conda install -y conda-build prodigal mcl muscle hmmer \
                     diamond blast megahit bowtie2 bwa \
                     samtools centrifuge trimal iqtree

# Build Anvi'o
RUN cd /tmp && \
    git clone --recursive --depth 1 https://github.com/merenlab/anvio.git && \
    cd anvio && \
    python setup.py sdist

RUN mkdir /tmp/recipe
RUN echo "{% set version = \"5.6\" %}" >> /tmp/recipe/meta.yaml
RUN echo "{% set sha256 = \"$(shasum -a 256 /tmp/anvio/dist/anvio-5.5-master.tar.gz | cut -d ' ' -f 1)\" %}" >> /tmp/recipe/meta.yaml
RUN echo "{% set url = \"file:///tmp/anvio/dist/anvio-5.5-master.tar.gz\" %}" >> /tmp/recipe/meta.yaml
RUN wget https://gist.githubusercontent.com/ozcan/7528b48853cf6af76bb34743036eaec0/raw/3fb48af408703c80dd8f48f744a1c41edb6adb04/meta.yaml -O /tmp/meta.yaml && \
    cat /tmp/meta.yaml >> /tmp/recipe/meta.yaml && \
    rm /tmp/meta.yaml
RUN conda-build /tmp/recipe

# Install Anvi'o
RUN conda index /opt/conda/envs/anvioenv/conda-bld/
RUN conda install -c file:///opt/conda/envs/anvioenv/conda-bld/ anvio-minimal=5.6

# Install METABAT and DAS_TOOL
RUN conda install metabat2 das_tool

# Install CONCOCT
RUN apt-get update && apt-get install -qq build-essential libgsl0-dev bedtools mummer samtools perl libssl-dev
RUN pip install https://github.com/BinPro/CONCOCT/archive/1.1.0.tar.gz

# Install BINSANITY
RUN pip install git+https://github.com/edgraham/BinSanity.git

# Install MAXBIN2 (installing fraggenescan will require cpanm, and we also need IDBA-UD)
RUN conda install -c bioconda perl-app-cpanminus
RUN cpanm --self-upgrade --sudo
RUN conda install -c bioconda idba
RUN cd /opt && wget https://downloads.sourceforge.net/project/fraggenescan/FragGeneScan1.31.tar.gz && tar zxf FragGeneScan1.31.tar.gz && cd FragGeneScan1.31 && make clean && make
RUN cpanm install LWP::Simple
ENV PERL5LIB /opt/conda/envs/anvioenv/lib/5.26.2/:/opt/conda/envs/anvioenv/lib/site_perl/5.26.2/:$PERL5LIB
RUN cd /opt && wget https://downloads.sourceforge.net/project/maxbin2/MaxBin-2.2.7.tar.gz && tar zxf MaxBin-2.2.7.tar.gz && cd MaxBin-2.2.7/src && make
ENV PATH /opt/FragGeneScan1.31:/opt/MaxBin-2.2.7:$PATH

# Install some helper tools
RUN pip install virtualenv
RUN apt-get install vim

# Setup the environment
ENV ANVIO_VERSION 6
RUN echo "export PS1=\"\[\e[0m\e[47m\e[1;30m\] :: anvi'o v$ANVIO_VERSION :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>>\[\e[m\] \[\e[0m\]\"" >> /root/.bashrc

CMD /bin/bash -l