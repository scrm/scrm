#!/bin/bash
#

pars=$@
[ $# == 0 ] && pars='5 10 -r 100 1000 -l 50'
make -j2 scrm

i=0
while [ 1 ]; do
  ((i++))
  echo "Seed: $i"
  ./scrm $pars -seed $i > /dev/null 
  [ $? -eq 0 ] || { echo "Failed: ./scrm $pars -seed $i"; exit 1; }
done
