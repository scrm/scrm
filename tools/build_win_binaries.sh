#!/bin/bash
# ---------------------------------------------------------------------------
# build_win_binaries.sh
# Cross-compiles scrm for Windows. Requires that mingw-w64 is installed.
# 
# usage ./build_win_binaries.sh vx.y.z
# ---------------------------------------------------------------------------


if [ "${1:0:1}" != "v" ]; then
  echo "Wrong argument; usage $0 vx.y.z"
  exit 1
fi

version=${1:1}
tempdir=$(mktemp -d)
cd $tempdir

echo "Getting scrm $1"
git clone -q https://github.com/scrm/scrm.git || exit 1
cd scrm || exit 1
git checkout -q "$1" || exit 1


./bootstrap

# Build the 64bit version
CXX=i686-w64-mingw32-g++ CXXFLAGS='-O3' LDFLAGS='-static-libgcc -static -lpthread' \
  ./configure --host=i686-w64-mingw32 || exit 1
make && zip ../scrm-$version-win64.zip scrm.exe

# Build the 32bit version
CXX=i686-w64-mingw32-g++ CXXFLAGS='-O3 -m32' LDFLAGS='-static-libgcc -static -lpthread' \
  ./configure --host=i686-w64-mingw32 || exit 1
make clean && make && zip ../scrm-$version-win32.zip scrm.exe

cd ..
rm -rf scrm
echo "Binaries in $tempdir"
