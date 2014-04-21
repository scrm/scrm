#!/bin/bash

qsub submit_scrm0.sh
qsub submit_scrm50000.sh
qsub submit_scrm10000.sh
qsub submit_scrm1000.sh
qsub submit_ms.sh
qsub submit_scrm100000.sh


paramfile=Constantpopsizeparam
./ld_test.py ${paramfile}
