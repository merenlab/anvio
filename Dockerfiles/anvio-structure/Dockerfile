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
ENV ANVIO_VERSION "7.1_structure"

SHELL ["/bin/bash", "--login", "-c"]

RUN conda config --env --add channels bioconda
RUN conda config --env --add channels conda-forge

# Create a conda environment for anvi'o, activate it, and make sure it will
# always be activated
RUN conda create -n anvioenv python=3.6
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

# Install anvi'o 7.1 packages
# NOTE R 3.6.1 or maybe 3.6.3 is desired
RUN mamba install -y -c bioconda -c conda-forge -c salilab python=3.6 \
        "sqlite>=3.31.1" prodigal mcl muscle=3.8.1551 hmmer \
        diamond blast megahit spades bowtie2 tbb=2019.8 bwa \
        samtools=1.9 centrifuge trimal iqtree trnascan-se \
        r-base=3.6.1 r-stringi r-tidyverse r-magrittr r-optparse \
        bioconductor-qvalue fasttree vmatch fastani \
        dssp=3.0.0 paml modeller

# install pymol
RUN mamba install -c schrodinger pymol

# install qvalue
RUN Rscript -e 'install.packages("BiocManager", repos="https://cran.rstudio.com"); BiocManager::install("qvalue")'

# Setup APT
RUN apt-get update

##############################################################
# Install anvi'o from pip
##############################################################
RUN curl -L https://github.com/merenlab/anvio/releases/download/v7.1/anvio-7.1.tar.gz -o anvio-7.1.tar.gz
RUN pip install file:///anvio-7.1.tar.gz
RUN rm anvio-7.1.tar.gz

##############################################################
# Setup anvi'o databases
##############################################################
RUN anvi-setup-interacdome
RUN anvi-setup-kegg-data --mode modules --kegg-snapshot v2020-12-23
RUN anvi-setup-pfams --pfam-version 33.1
RUN anvi-setup-ncbi-cogs --cog-version COG20

##############################################################
# install dependencies for the reproducible workflow
##############################################################
RUN pip install xlrd==1.2.0 openpyxl==3.0.9 prody==2.0.1 tmscoring 'iminuit<2'

##############################################################
# Other dependencies
##############################################################

# Install METABAT and DAS_TOOL
RUN conda install metabat2 das_tool

# Install CONCOCT
RUN apt-get install -qq build-essential libgsl0-dev bedtools mummer samtools perl libssl-dev
RUN mamba install cython
RUN pip install https://github.com/BinPro/CONCOCT/archive/1.1.0.tar.gz

# Install some helper tools
RUN apt-get install vim util-linux -yy

# Cutify the environment
RUN echo "export PS1=\"\[\e[0m\e[41m\e[1;31m\] :: anvi'o v7.1 STRUCTURE :: \[\e[0m\e[0m \[\e[1;34m\]\]\w\[\e[m\] \[\e[1;32m\]>>>\[\e[m\] \[\e[0m\]\"" >> /root/.bashrc

CMD /bin/bash -l

# To test the build, you can run this:
#
# docker run --rm -it -v `pwd`:`pwd` -w `pwd` -p 8080:8080 meren/anvio:test-build
