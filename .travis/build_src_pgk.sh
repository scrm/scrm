#!/bin/bash

# Build the manual
sudo apt-get install -y pandoc || exit 1
./doc/create-manual.sh

# Build the release
CXXFLAGS="-O3" ./bootstrap || exit 1
make clean && make dist-check || exit 1
mv scrm-*.tar.gz scrm-src.tar.gz || exit 1
