#!/bin/bash

set -e

############################
# build third party tools
############################
cd third_party_tools

# fastp
tar -xvzf fastp.tar.gz
make -C fastp -j

# bwa 
tar -xvzf bwa.tar.gz
make -C bwa -j

# samtools
tar -xvzf htslib.tar.gz
tar -xvzf samtools.tar.gz
make -C samtools -j

# megahit
tar -xvzf megahit_v1.2.9.tar.gz
if [[ -d megahit-1.2.9/build ]]; then
    rm -rf megahit-1.2.9/build
fi
mkdir -p megahit-1.2.9/build
cmake -S megahit-1.2.9 -B megahit-1.2.9/build
make -C megahit-1.2.9/build -j

# spades
tar -xzf SPAdes-3.15.5-Linux.tar.gz
if [[ -d SPAdes-3.15.5 ]]; then
    rm -rf SPAdes-3.15.5
fi
mv SPAdes-3.15.5-Linux SPAdes-3.15.5
rm -rf SPAdes-3.15.5-Linux

# DeepVirFinder
tar -xvzf DeepVirFinder.tar.gz

cd -
