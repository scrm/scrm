#!/bin/bash
#
# This script checks if the current version of scrm produces the
# same results as the system wide installed version.
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2014-07-21
# Licence:  GPLv3 or later
#

function test_scrm {
  echo -n " scrm $@ "
  for i in `seq 1 10`; do
    echo -n "."

    hash_old=$(scrm $@ -seed $i | sha512sum)
    hash_new=$(./scrm $@ -seed $i | sha512sum)

    if [ "$hash_old" != "$hash_new" ]; then
      echo
      echo "Different results for scrm $@ -seed $i" 
      exit 1
    fi
  done
  echo " done."
}

test_scrm 10 10 -r 10 100 -t 5 -T
test_scrm 2000 10 -r 10 100 -l 10 -t 2 -T
