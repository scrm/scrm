#!/bin/bash

# Determine if the current commit is tagged
tag=`git describe --exact-match HEAD 2> /dev/null` || exit 0
[ "${tag:0:1}" == "v" ] || exit 0 
version=${tag:1}

# Assert that we only do this once
[ "${CXX:='g++'}" == "g++" ] || exit 0

# Build the manual
sudo apt-get install -y pandoc || exit 1
./doc/create-manual.sh

# Build the release
echo "Building release $tag..."
CXXFLAGS="-O3" ./bootstrap || exit 1
make dist || exit 1
release=`ls scrm-*.tar.gz` || exit 1

# Remove gcc-4.8
# sudo apt-get install -y mingw-w64 || exit 1

# Build Windows Binaries
# Build the 64bit version
#CXX=i686-w64-mingw32-g++ CXXFLAGS='-O3' LDFLAGS='-static-libgcc -static -lpthread' \
#  ./configure --host=i686-w64-mingw32 || exit 1
#make && zip scrm-$version-win64.zip scrm.exe

# Build the 32bit version
#CXX=i686-w64-mingw32-g++ CXXFLAGS='-O3 -m32' LDFLAGS='-static-libgcc -static -lpthread' \
#  ./configure --host=i686-w64-mingw32 || exit 1
#make clean && make && zip scrm-$version-win32.zip scrm.exe

# Get & set up the homepage repo
echo "Setting up git..."
git clone https://github.com/scrm/scrm.github.io.git
cd scrm.github.io
git config user.email "${GH_EMAIL}"
git config user.name "${GH_NAME}"
echo "https://${GH_TOKEN}:@github.com" > .git/credentials
git config credential.helper "store --file=.git/credentials"

# Copy tar.gz into releases folder and generate hashes
echo "Adding files..." 
mv ../$release releases/ || exit 1
mv ../scrm-$version-* releases/ || exit 1
rm releases/scrm-current*
cd releases
ln -s $release scrm-current.tar.gz
#ln -s scrm-$version-win64.zip scrm-current-win64.zip
#ln -s scrm-$version-win32.zip scrm-current-win32.zip
cd ..
make
git add releases/$release releases/scrm-$version-* releases/scrm-current-* releases/releases*

# And push everything to GitHub
echo "Pushing to GitHub..."
git commit -m "Add release $version"
git push -fq origin master 2>&1 > /dev/null || exit 1
rm .git/credentials
echo "Done"
