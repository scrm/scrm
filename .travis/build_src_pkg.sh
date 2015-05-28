#!/bin/bash

# Build the manual
sudo add-apt-repository -y ppa:marutter/c2d4u || exit 1
sudo apt-get update -qq && sudo apt-get install -qqy pandoc || exit 1
./doc/create-manual.sh || exit 1

# Build the release
CXXFLAGS="-O3" ./bootstrap || exit 1
make clean && make distcheck || exit 1
mv scrm-*.tar.gz scrm-src.tar.gz || exit 1
