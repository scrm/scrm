#!/bin/bash

# Determine if the current commit is tagged
#tag=`git describe --exact-match HEAD 2> /dev/null` || exit 0
tag='TestRelease'

# Assert that we only do this once
[ "${CXX:='g++'}" == "g++" ] || exit 0

# Build the release
echo "\nBuilding release $tag..."
CXXFLAGS="-O3" ./bootstrap || exit 1
make distcheck || exit 1
release=`ls scrm-*.tar.gz` || exit 1

# Get & set up the homepage repo
echo "\nSetting up git..."
git clone https://github.com/scrm/scrm.github.io.git
cd scrm.github.io
git config user.email "${GH_EMAIL}"
git config user.name "${GH_NAME}"
echo "https://${GH_TOKEN}:@github.com" > .git/credentials
git config credential.helper "store --file=.git/credentials"

# Copy tar.gz into releases folder and generate hashes
mv $release releases/
rm releases/scrm-current*
ln -s releases/$release releases/scrm-current.tar.gz
make
git add releases/$release releases/scrm-current.tar.gz releases/releases*

# And push everything to GitHub
echo "\nPushing to GibHub..."
git commit -m "Add release $tag"
git push -fq origin master 2>&1 > /dev/null || exit 1
rm .git/credentials
echo "Done"
