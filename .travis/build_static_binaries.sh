#!/bin/bash

# Build statically linked binaries
CXXFLAGS="-O3" LDFLAGS='-static' ./configure || exit 1 
make clean && make && \
  tar -zcvf "scrm-x64-static.tar.gz" scrm doc/manual.html &&
  cp scrm scrm-x64-static

