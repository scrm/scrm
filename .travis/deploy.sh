#!/bin/bash

# Determine if the current commit is tagged
# tag=`git describe --exact-match HEAD 2> /dev/null` || exit 0

# Assert that we only do this once
CXX=${CXX:="g++"}
[ "$CXX" == "g++" ] || exit 0

# If it is, build the .tar.gz
echo "Building release $tag..."
make distcheck || exit 1
release=`ls scrm-*.tar.gz` || exit 1

# Copy the generated .tar.gz to the gh-pages branch
git config --global user.email "${GH_EMAIL}"
git config --global user.name "${GH_NAME}"
git remote rm origin
git remote add origin "https://${GH_TOKEN}@github.com/paulstaab/scrm"
git checkout --track origin/gh-pages

mv $release releases/
cd releases
md5sum $release > $release.md5
sha512sum $release > $release.sha512
git add $release*

# Create symlinks for current release
rm scrm-current*
ln -s $release scrm-current.tar.gz
md5sum scrm-current.tar.gz > scrm-current.tar.gz.md5
sha512sum scrm-current.tar.gz > scrm-current.tar.gz.sha512
git add scrm-current.*

git commit -m "Add release $tag"

# And push everything to GitHub
git push origin gh-pages
