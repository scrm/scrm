#!/bin/bash
#
# %title%
# %description%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-01-14
# Licence:  GPLv3 or later
#

cd src/
model="100 5 -t 5 -r 7 1001"

echo "SCRM:"
time for i in `seq 1 100` 
do
  ./scrm $model -l 100 > /dev/null 
done

echo "MS:"
time for i in `seq 1 100` 
do
  ms $model > /dev/null 
done
