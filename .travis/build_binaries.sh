#!/bin/bash

# Build statically linked binaries
CXXFLAGS="-O3" LDFLAGS='-static' ./configure || exit 1 
make clean && make && \
  tar -zcvf "scrm-$version-x64-static.tar.gz" scrm doc/manual.html

# Build Windows Binaries
CXX=i686-w64-mingw32-g++ CXXFLAGS='-O3' LDFLAGS='-static-libgcc -static -lpthread' \
  ./configure --host=i686-w64-mingw32 || exit 1
make clean && make && zip scrm-$version-win64.zip scrm.exe

# Build the 32bit version
CXX=i686-w64-mingw32-g++ CXXFLAGS='-O3 -m32' \
  LDFLAGS='-static-libgcc -static -lpthread' \
  ./configure --host=i686-w64-mingw32 || exit 1
make clean && make && zip scrm-$version-win32.zip scrm.exe
