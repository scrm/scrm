#!/bin/bash

os=$1
cxx=$2
cxxflags=$3

echo "os: $os; cxx: $cxx; cxxflags: $cxxflags"

if [ "$os" == "linux" ]; then
  if [[ "$cxxflags" == *"-m32"* ]]; then
    echo "On 32bit. Skipping coverage test."
  else
    coveralls --exclude lib --exclude tests --gcov-options '\-lp'
  fi
fi

if [ "$os" == "osx" ]; then
    echo "On OS X. Skipping coverage test."
fi

