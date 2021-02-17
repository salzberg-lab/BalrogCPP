# Bacterial Annotation by Learned Representation Of Genes (C++ version)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/installer/conda.svg)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/balrog.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/platforms.svg)](https://anaconda.org/bioconda/balrog)
[![BioConda Install](https://anaconda.org/bioconda/balrog/badges/license.svg)](https://anaconda.org/bioconda/balrog)

## Overview
Balrog is a prokaryotic gene finder based on a Temporal Convolutional Network. We took a data-driven approach to prokaryotic gene finding, relying on the large and diverse collection of already-sequenced genomes. By training a single, universal model of bacterial genes on protein sequences from many different species, we were able to match the sensitivity of current gene finders while reducing the overall number of gene predictions. Balrog does not need to be refit on any new genome.

## Publication
Preprint available on bioRxiv [here](https://www.biorxiv.org/content/10.1101/2020.09.06.285304v1).

PLOS Computational Biology publication coming soon...

## Install Balrog via conda (relatively slow performance)
Currently, the bioconda version of Balrog does not take advantage of AVX2 instructions due to known PyTorch/LibTorch issues and is ~5-10x slower than the version built from source. We are working on a fix.

    conda create -n balrog_env python=3.7 -y
    
    conda activate balrog_env
    (alternatively: "source activate balrog_env")
    
    conda install balrog -c conda-forge -c bioconda -y
    conda install pytorch=1.7.1 -c conda-forge -y
    
    balrog --help


## Compile Balrog from source (HIGHLY RECOMMENDED for now for speed)

### Install MMseqs2
Balrog depends on MMseqs2 at runtime to help reduce false positive gene predictions. Fortunately, MMseqs2 is well supported on both Linux and MacOS. Detailed installation instructions for MMseqs2 can be found on the MMseqs2 GitHub [here](https://github.com/soedinglab/MMseqs2#installation)

    # install by brew
    brew install mmseqs2
    # install via conda
    conda install -c conda-forge -c bioconda mmseqs2
    # install docker
    docker pull soedinglab/mmseqs2
    # static build with AVX2 (fastest)
    wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
    # static build with SSE4.1
    wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz; tar xvfz mmseqs-linux-sse41.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
    # static build with SSE2 (slowest, for very old systems)
    wget https://mmseqs.com/latest/mmseqs-linux-sse2.tar.gz; tar xvfz mmseqs-linux-sse2.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH

### Build Balrog
    # Linux
    git clone https://github.com/salzberg-lab/BalrogCPP
    cd BalrogCPP
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${PREFIX} ..
    make
    export PATH=$(pwd):$PATH
    
    # MacOS
    git clone https://github.com/salzberg-lab/BalrogCPP
    cd BalrogCPP
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${PREFIX} ..
    make
    export PATH=$(pwd):$PATH

Balrog also requires zlib, which is likely already installed on your system. 

    # Most Linux systems 
    sudo apt install zlib
    
    # Ubuntu
    sudo apt-get install zlib1g
    sudo apt-get install zlib1g-dev

    # MacOS
    brew install zlib
    




