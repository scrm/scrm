#!/bin/bash

# Determine if the current commit is tagged
tag=`git describe --exact-match HEAD 2> /dev/null` || exit 0

# If it is, build the .tar.gz
echo "Building release $tag..."
make distcheck
release=`ls scrm-*.tar.gz` || exit 1

# Copy the generated .tar.gz to the gh-pages branch
git checkout --track origin/gh-pages
mv $release releases/
git add releases/$release
git commit -m "Add release $tag"
