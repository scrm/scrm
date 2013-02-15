#!/bin/bash
#
# %title%
# %description%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-02-01
# Licence:  GPLv3 or later
#

count=0
i=0

while [ 1 ]; do
  ((i++))
  ./scrm 2>&1 > tmp/run_$i.log
  count=$?
  echo "$i: return: $count"
  [ $count -eq 0 ] || exit 1
  juhu=$(grep -c "JUHU" tmp/run_$i.log)
  echo "   juhus: $juhu"
  [ $juhu -eq 0 ] || exit 1
  sleep 2
done
