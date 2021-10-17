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

FROM continuumio/miniconda3:4.9.2
ENV ANVIO_VERSION "7.1-rc"

RUN conda config --env --add channels bioconda
RUN conda config --env --add channels conda-forge

RUN conda create -n anvioenv python=3.6

# Install nano early on so the container has
# something to edit files if debug is necessary
RUN conda install -y nano

# Activate environment
ENV PATH /opt/conda/envs/anvioenv/bin:$PATH
ENV CONDA_DEFAULT_ENV anvioenv
ENV CONDA_PREFIX /opt/conda/envs/anvioenv

RUN echo "conda activate anvioenv" >> ~/.bashrc

RUN conda install -y conda-build
RUN conda install -y conda-verify

COPY conda-recipe /tmp/conda-recipe

# build and install anvio-minimal
RUN conda-build /tmp/conda-recipe/anvio-minimal
RUN conda index /opt/conda/envs/anvioenv/conda-bld/
RUN conda install -c file:///opt/conda/envs/anvioenv/conda-bld/ anvio-minimal=$ANVIO_VERSION
RUN conda build purge-all

# build and install anvio meta package. please note that the next line is quite a memory
# intensive step. if your docker build command is getting killed abruptly, you may want
# to increase the memory docker allowed to use (from the docker preferences or otherwise).
RUN conda-build /tmp/conda-recipe/anvio
RUN conda index /opt/conda/envs/anvioenv/conda-bld/
RUN conda install -c file:///opt/conda/envs/anvioenv/conda-bld/ anvio=$ANVIO_VERSION
RUN conda build purge-all

# Install METABAT and DAS_TOOL
RUN conda install metabat2 das_tool

# Install CONCOCT
RUN apt-get update && apt-get install -qq build-essential libgsl0-dev bedtools mummer samtools perl libssl-dev
RUN conda install cython
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
RUN echo 'export PATH=/opt/FragGeneScan1.31:$PATH' >> ~/.bashrc
RUN echo 'export PATH=/opt/MaxBin-2.2.7:$PATH' >> ~/.bashrc

# Install some helper tools
RUN pip install virtualenv
RUN apt-get install vim util-linux -yy

# this one is too painful to install through Conda
RUN Rscript -e 'install.packages(c("optparse"), repos="https://cran.rstudio.com")'

# Cutify the environment
RUN echo "export PS1=\"\[\e[0m\e[47m\e[1;30m\] :: anvi'o v$ANVIO_VERSION :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>>\[\e[m\] \[\e[0m\]\"" >> /root/.bashrc

CMD /bin/bash -l

# To test the build, you can run this:
#
# docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:test-build
