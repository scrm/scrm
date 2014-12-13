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
  time_file_stable=$(mktemp)
  time_file_new=$(mktemp)
  for i in `seq 1 100`; do

    equal_results='Yes'
    hash_stable=$(/usr/bin/time -f "%e %M" scrm $@ -seed $i 2>> $time_file_stable | sha512sum)
    hash_new=$(/usr/bin/time -f "%e %M" ./scrm $@ -seed $i 2>> $time_file_new | sha512sum)

    if [ "$hash_stable" != "$hash_new" ]; then
      equal_results='No'
    fi
  done
  time_stable=`Rscript -e "options(digits=2); cat(mean(read.table(\"$time_file_stable\")[,1]), 's\n', sep='')"`
  ram_stable=`Rscript -e "cat(round(mean(read.table(\"$time_file_stable\")[,2])), 'kb\n', sep='')"`
  time_new=`Rscript -e "options(digits=2); cat(mean(read.table(\"$time_file_new\")[,1]), 's\n', sep='')"`
  ram_new=`Rscript -e "cat(round(mean(read.table(\"$time_file_new\")[,2])), 'kb\n', sep='')"`
  rm "$time_file_stable" "$time_file_new"

  echo "scrm $@ | $equal_results | $time_stable | $time_new | $ram_stable | $ram_new"
}

echo "Command | Equal results | time stable | time new | memory stable | memory new"
echo "------- | ------------- | ----------- | -------- | ------------- | ---------"
test_scrm 100 100 -t 5 -T
test_scrm 20 100 -r 10 100 -t 3 -T -L
test_scrm 4 10000 -r 10 100 -t 5
test_scrm 20 2 -r 500 100 -L
test_scrm 2000 1 -r 10 100 -l 50 -L
test_scrm 20 1 -r 4000 10000000 -l 300000 -T
test_scrm 10 100 -r 10 100 -I 3 3 3 4 0.5 -eN 0.1 0.05 -eN 0.2 0.5
test_scrm 50 100 -r 20 200 -G 1.5 -l 100
