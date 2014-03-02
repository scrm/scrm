#!/bin/bash

# Determine if the current commit is tagged
# tag=`git describe --exact-match HEAD 2> /dev/null` || exit 0

# Assert that we only do this once
[ "${CXX:='g++'}" == "g++" ] || exit 0

# Install Latex to build the documentation
sudo apt-get install --no-install-recommends \
  texlive-base texlive-latex-base texlive-generic-recommended \
  texlive-fonts-recommended texlive-fonts-extra texlive-extra-utils \
  texlive-latex-recommended texlive-latex-extra texinfo lmodern

# If it is, build the .tar.gz
echo "Building release $tag..."
CXXFLAGS="-O3" ./bootstrap
make distcheck || exit 1
release=`ls scrm-*.tar.gz` || exit 1

# Configure git and make the repo use https-auth
git config --global user.email "${GH_EMAIL}"
git config --global user.name "${GH_NAME}"
git remote rm origin
git remote add origin "https://${GH_TOKEN}@github.com/paulstaab/scrm"
git fetch origin gh-pages
git checkout gh-pages

# Copy tar.gz into releases folder and generate hashsums
mv $release releases/
cd releases
md5sum $release > $release.md5
sha512sum $release > $release.sha512
git add $release*

# Also create a symlink point to the current release
rm scrm-current*
ln -s $release scrm-current.tar.gz
md5sum scrm-current.tar.gz > scrm-current.tar.gz.md5
sha512sum scrm-current.tar.gz > scrm-current.tar.gz.sha512
git add scrm-current.*

# And push everything to GitHub
git commit -m "Add release $tag"
git push origin gh-pages
