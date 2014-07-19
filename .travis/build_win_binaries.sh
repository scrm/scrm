#!/bin/bash

# 64 bit
CXX=i686-w64-mingw32-g++ CXXFLAGS='-O3' LDFLAGS='-static-libgcc -static -lpthread' \
  ./configure --host=i686-w64-mingw32 || exit 1
make clean && make && zip scrm-win64.zip scrm.exe doc/manual.html

# 32bit 
CXX=i686-w64-mingw32-g++ CXXFLAGS='-O3 -m32' \
  LDFLAGS='-static-libgcc -static -lpthread' \
  ./configure --host=i686-w64-mingw32 || exit 1
make clean && make && zip scrm-win32.zip scrm.exe doc/manual.html
