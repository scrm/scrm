#!/bin/bash

#-------------------------------------------------#
# Builds scrm's manual from scrm's wiki on GitHub #
#-------------------------------------------------#

pages="scrm.wiki/Installation.md scrm.wiki/Command-Line-Options.md scrm.wiki/Output.md"

scrm_dir=$PWD
if [ ! -d "$scrm_dir/doc" ]; then
  "Error: please execute this script from scrm's folder"
  exit 1
fi

tmpdir=$(mktemp -d)
cd "$tmpdir" || exit 1
git clone https://github.com/scrm/scrm.wiki.git || exit 1
cp "$scrm_dir/doc/knitr.css" ./ || exit 1

pandoc -t html5 -s -S --toc --toc-depth=2 \
  -M title='The scrm Quick Reference' \
  -M author='Paul R. Staab, Joe (Sha) Zhu, Dirk Metzler and Gerton Lunter' \
  -c knitr.css \
  --self-contained \
  -o manual.html ${pages}

cp -v manual.html "$scrm_dir/doc/"
cd "$scrm_dir"
rm -rf "$tmpdir"
