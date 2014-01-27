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

model="10 1"
points=`seq 1 1 100`

echo "ARG"
echo "bases,runtime,memory" > run_times_arg.csv
for i in $points 
do
  echo $i
  for j in `seq 1 50` 
  do
    /usr/bin/time -f "$((i)),%e,%M" scrm $model -r $((i*25)) $((i*1000)) > /dev/null 2>> run_times_arg.csv 
  done
done

echo "\nAPPROX"
echo "bases,runtime,memory" > run_times_approx.csv
for i in $points 
do
  echo $i
  for j in `seq 1 50` 
  do
    /usr/bin/time -f "$((i)),%e,%M" scrm $model -r $((i*25)) $((i*1000)) -l 10000 > /dev/null 2>> run_times_approx.csv 
  done
done

echo "\nMS"
echo "bases,runtime,memory" > run_times_ms.csv
for i in $points 
do
  echo $i
  for j in `seq 1 50` 
  do
    /usr/bin/time -f "$((i)),%e,%M" ms $model -r $((i*25)) $((i*1000)) -t 5 > /dev/null 2>> run_times_ms.csv 
  done
done
grep -v "Command" run_times_ms.csv > run_times_ms2.csv
