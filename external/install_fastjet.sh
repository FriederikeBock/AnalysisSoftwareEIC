#!/usr/bin/env bash

# Based on: https://github.com/scikit-hep/pyjet/blob/master/install-fastjet.sh

# If you need to specify CC and CXX (for example, if there are multiple compilers),
# you can do so here.
# NOTE: This isn't required on the ORNL cluster anymore
#export CC=gcc-10
#export CXX=g++-10

# Setup
export CXXFLAGS="-std=c++17"
fastjet_version=3.3.4
fjcontrib_version=1.045
prefix=$PWD/install/fastjet

# Build in a new folder
mkdir -p fastjet
cd fastjet

if [ ! -d fastjet-${fastjet_version} ]; then
    wget http://fastjet.fr/repo/fastjet-${fastjet_version}.tar.gz
    tar xfz fastjet-${fastjet_version}.tar.gz
fi

if [ ! -d fjcontrib-${fjcontrib_version} ]; then
    wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-${fjcontrib_version}.tar.gz
    tar xfz fjcontrib-${fjcontrib_version}.tar.gz
fi

cd fastjet-${fastjet_version}
make clean
# NOTE: Need to disable autoptr because we're using c++17
./configure --prefix=$prefix --enable-allcxxplugins --enable-all-plugins --disable-auto-ptr
make -j4
make install
cd ../fjcontrib-${fjcontrib_version}
make clean
./configure --prefix=$prefix --fastjet-config=$prefix/bin/fastjet-config
make -j4
make install
make fragile-shared-install
