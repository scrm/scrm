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

count=0
i=0

pars=$@
[ $# == 0 ] && pars=5

while [ 1 ]; do
  ((i++))
  echo "Seed: $i"
  ./src/scrm_dbg $pars -seed $i > /dev/null 
  [ $? -eq 0 ] || exit 1
done
