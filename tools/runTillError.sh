#!/bin/bash
#
# runTillError.sh
# Systematically calls the debug version of scrm for different seeds until an
# error occurs. Passes command line parameters to scrm. 
#
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-04-19
# Licence:  GPLv3 or later
#

pars=$@
[ $# == 0 ] && pars='5 10 -r 100 1000 -l 50'
make -j2 scrm scrm_dbg

i=0
while [ 1 ]; do
  ((i++))
  echo "Seed: $i"
  ./scrm_dbg $pars -seed $i > /dev/null 
  [ $? -eq 0 ] || { echo "Failed: ./scrm_dbg $pars -seed $i"; exit 1; }
  valgrind --error-exitcode=1 --leak-check=full -q ./scrm $pars -seed $i > /dev/null  
  [ $? -eq 0 ] || { echo "valgrind failed: ./scrm_dbg $pars -seed $i"; exit 1; }
done
