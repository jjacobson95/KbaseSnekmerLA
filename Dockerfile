# use multi-stage build to install snekmer and avoid package dependency issues
# environment.yml needs to be updated manually if the snekmer version is updated here
# FROM condaforge/mambaforge AS mambasetup
# COPY ./environment.yml .
# RUN mamba env create -f ./environment.yml

# combine mambasetup build into kbase image
FROM kbase/sdkbase2
MAINTAINER KBase Developer

# share packages from the mambasetup stage build
# the path must also be updated in the entrypoint.sh
# COPY --from=mambasetup /opt/conda/envs/snekmer/. /opt/conda/envs/snekmer/
# ENV PATH /opt/conda/envs/snekmer/bin:$PATH


# Update the sources.list to use archived repositories
RUN echo "deb http://archive.debian.org/debian/ stretch main contrib non-free" > /etc/apt/sources.list && \
    echo "deb http://archive.debian.org/debian-security/ stretch/updates main contrib non-free" >> /etc/apt/sources.list

RUN apt-get update && apt-get install -y \
    build-essential \
    libpython3-dev \
    git \
    software-properties-common

RUN apt-get update && apt-get install -y \
    build-essential \
    libffi-dev \
    libssl-dev \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    wget \
    curl \
    llvm \
    libncurses5-dev \
    xz-utils \
    tk-dev \
    libxml2-dev \
    libxmlsec1-dev \
    libffi-dev \
    liblzma-dev

# Download and compile Python 3.9
ENV PYTHON_VERSION=3.9.1
RUN curl -O https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tar.xz \
    && tar -xf Python-${PYTHON_VERSION}.tar.xz \
    && cd Python-${PYTHON_VERSION} \
    && ./configure --enable-optimizations \
    && make -j `nproc` \
    && make altinstall \
    && cd .. \
    && rm -rf Python-${PYTHON_VERSION} \
    && rm Python-${PYTHON_VERSION}.tar.xz

# Other installations
RUN apt-get install -y git python3-pip

# Clone Snekmer repository
RUN git clone https://github.com/PNNL-CompBio/Snekmer.git /opt/Snekmer
WORKDIR /opt/Snekmer
RUN git checkout kmer-association

RUN pip3.9 install --upgrade pip
RUN pip3.9 install -r requirements.txt
RUN pip3.9 install -e git+https://github.com/PNNL-CompBio/Snekmer@2c2db4b17d68764d229920287283392c9b3dc743#egg=snekmer

# kbase sdk code
COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module
WORKDIR /kb/module
RUN make all
# ENTRYPOINT [ "./scripts/entrypoint.sh" ]
# CMD [ ]

# RUN source activate root