#!/bin/bash

os=$1
cxx=$2
cxxflags=$3

echo "os: $os; cxx: $cxx; cxxflags: $cxxflags"

if [ "$os" == "linux" ]; then
  if [[ "$cxxflags" == *"-m32"* ]]; then
    sudo dpkg --add-architecture i386
    sudo apt-get update -qq
    sudo apt-get install -qqy automake valgrind g++-multilib libc6-dbg:i386 libcppunit-dev:i386
  else
    sudo apt-get update -qq
    sudo apt-get install -qqy automake libcppunit-dev valgrind;
  fi
fi

if [ "$os" == "osx" ]; then
  brew update
  brew install cppunit valgrind
fi

